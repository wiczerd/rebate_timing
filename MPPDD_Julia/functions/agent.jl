"""
Return the asset grids, given the total (labor) income grid, T, and N_k
α = 1 for uniformly spaced gripoints.
"""
function assetGrid(T::Int64, income_grid::Array{Float64,3}, N_k::Int64, r::Float64,  α::Int64 = 1)

    asset_max  = zeros(T+1);
    asset_min  = zeros(T+1);
    asset_grid = zeros(T+1, N_k); # T + 1 = year after death, define for convenience in backwardSolve

    # Minimum assets (Maximum debt)
    for t = T:-1:2
        income_min   = minimum(income_grid[t,:,:]);
        asset_min[t] = asset_min[t+1]/(1 + r) - income_min + 1e-4;
    end

    # Make k_max the largest amount the consumer could ever save, even
    # if they got the highest income and saved half of income every quarter
    for t = 1:T
        income_max     = maximum(income_grid[t,:,:]);
        asset_max[t+1] = min(asset_max[t]*(1 + r) + 0.5*income_max - asset_min[t+1], 100)
    end

    for t = 1:T
        N_p = ceil(Int64, N_k*0.70);
        N_n = N_k - N_p;

        assets_p  = collect(LinRange(0.0, asset_max[t]^(1/α), N_p)).^α;  # optional convex grid for positive assets
        assets_n  = collect(LinRange(asset_min[t], 0, N_n));             # uniformly spaced points for negative assets

        # Asset grid for the given time period (linear below 0, nonlinear above 0)
        asset_grid[t, 1:N_n]     = assets_n;
        asset_grid[t, N_n+1:N_k] = assets_p;
    end

    asset_grid[T,: ] = collect(LinRange(asset_min[T], asset_max[T]^(1/α), N_k)).^α;

    return asset_grid
end

"""
Define the ϕ function
"""
@inline function ϕ(s::Float64, ϕ_1 ::Float64,ϕ_2::Float64)
    return ϕ_1 *(1 -  exp(-ϕ_2*(s-1)))/(1-exp(ϕ_2))
end

"""
CRRA utility
"""
@inline function u(c::Float64, γ::Int64)
    # Nudge to avoid log(0) - type errors during runtime
    if γ == 1
        return log(max(1e-8, c - cbar))
    else
        return ((c-cbar)^(1 - γ) - 1)/(1 - γ)
    end
end

"""
Define current utility + continuation value
"""
function objective(b::Float64, s::Float64, y::Float64, qhr::Array{Float64,1},
    asset_grid::Array{Float64,1}, a::Float64, EV_grid::Array{Float64,1}, ψ::Float64, β::Float64,
    ϕ_1::Float64, ϕ_2::Float64, γ::Int64)
    # Note: b, s, y, ψ are levels. a0idx is the index of a0.
    # qhr, asset_grid, EV are an arrays
    # Get level of capital implied by next-period consumption choice
    ap = b + (1.0 - s)*a;
    if ap>asset_grid[length(asset_grid)] || ap<asset_grid[1]
        println("ap out of bounds: (ap,b,s,a)= ($ap, $b, $s, $a)")
    end
    EV = LinearInterpolation(asset_grid, EV_grid)(ap);
    c = y + a*s - LinearInterpolation(asset_grid, qhr)(ap)*b ;
    return u(c, γ) - ϕ(s, ϕ_1, ϕ_2) + β*ψ*EV
end

"""
Solve for optimal b, given
"""
@inline function max_objective(s::Float64, y::Float64, qhr::Array{Float64,1},
     asset_grid::Array{Float64,1}, a::Float64, EV_grid::Array{Float64,1}, ψ::Float64,
     β::Float64, ϕ_1_1::Float64, ϕ_1_2::Float64, γ::Int64, bmin::Float64, bmax::Float64)

     # bstard[di] = max_gss(b_min, b_max, (s, i, y, self.phi, qhr, psi_t, asset_grid[t + 1], EV_grid))
     # default routine is Brent's Method
     res = optimize(b -> -1*objective(b, s, y, qhr, asset_grid, a, EV_grid, ψ, β, ϕ_1_1, ϕ_1_2, γ), bmin, bmax);
     return Optim.minimizer(res);
end

"""
Find s* by maximizing the interpolated v*(s),
"""
@inline function maxs(s_grid::Array{Float64,1}, vstars::Array{Float64,1},)
    # maxd = lambda dhr : -1.0*np.interp( np.subtract(np.ones(N_delinq), self.s_grid), vstard, dhr)
    # sopt, vopt  = opt.minimize(maxd, brack=(0, 1))
    #itp = LinearInterpolation(ones(N_delinq)-s_grid, vstard; extrapolation_bc=Throw());
    itp = LinearInterpolation(s_grid, vstars; extrapolation_bc=Throw());
    res = optimize(s-> -1*itp(s), 0, 1); # default routine is Brent's Method
    return Optim.minimizer(res), -1*Optim.minimum(res)
end

"""
Taking Q as given, solve for V,B,C,S,K
"""
function backwardSolve(rebate::Float64, t_rebate::Int64, β::Float64,
    γ::Int64, ϕ_1::Float64, ϕ_2::Float64, r::Float64, ε_probs::Array{Float64,2},
    z_probs::Array{Float64,3}, income_grid::Array{Float64,3}, s_grid::Array{Float64,1},
    asset_grid::Array{Float64,2}, T::Int64, N_k::Int64, N_s::Int64, N_z::Int64, N_ε::Int64,
    cbar::Float64, Q::Array{Float64,3})

    V = zeros(T+1, N_k, N_ε, N_z);
    C = zeros(T, N_k, N_ε, N_z);
    S = zeros(T, N_k, N_ε, N_z);
    B = zeros(T, N_k, N_ε, N_z);
    Ap= zeros(T, N_k, N_ε, N_z);

    for t =  T:(-1):1

        # a' grid
        at  = asset_grid[t + 1, :];
        a0  = asset_grid[t, :];

        # Rebate
        reb = rebate*(t == t_rebate);

        # Exogenous survival prob for the period
        ψ = survivalprobs[t];


        #for i = 1:N_z
        @inbounds Threads.@threads for yidx = 1:(N_z*N_ε)
        #for yidx = 1:(N_z*N_ε)
            j = floor(Int,mod(yidx-1,N_z*N_ε)/N_z )+1
            i =           mod(yidx-1,N_z)+1

            qhr    = Q[t, i,:];
            #for j = 1:N_ε
                y = income_grid[t, j, i] + reb;

                # Income transition probabilities
                tmat = ε_probs[t,:]*z_probs[t,i,:]'; # N_ε x N_z

                @inbounds for k = 1:N_k #loop over a0

                    bstars = zeros(N_s);
                    vstars = zeros(N_s);

                    # Right before certain death, the agent chooses to consume all that they can
                    if t == T
                        S[t, k, j, i] = 1.0;
                        B[t, k, j, i] = 0.0;
                        C[t, k, j, i] = a0[k] + y + reb;
                        V[t, k, j, i] = u(C[t, k, j, i], γ);
                        Ap[t, k, j, i] = 0.0;
                        continue
                    end

                    # Expected value grid next period, given income states and assets
                    EV_grid = [dot(V[t + 1, ik, :, :], tmat) for ik = 1:N_k];

                    if a0[k] >= 0 || t>=(T-1)
                        srange_a0 = N_s;
                    else
                        srange_a0 = 1;
                    end
                    # this is actually a choice -- we're going to loop over the best
                    @inbounds for si = srange_a0:N_s

                        s    = s_grid[si];

                        # b cannot put a' outside bounds and can't make c negative:
                        bmax = min(at[N_k] - (1.0 - s) * a0[k] - eps(Float32)*10, (y+a0[k]*s) *(1.0+r));
                        bmin = at[1]   - (1.0 - s) * a0[k] + eps(Float32)*10;
                        # c = y + a0[k]*s - bmin /(1.0 + r);  # check this is positve?

                        # Note: EV interpolation should take place over assets_{t+1} grid
                        bstars[si]  = max_objective(s, y, qhr, at, a0[k], EV_grid, ψ, β, ϕ_1 , ϕ_2, γ, bmin, bmax);
                        vstars[si]  = objective(bstars[si], s, y, qhr, at, a0[k], EV_grid, ψ, β, ϕ_1, ϕ_2, γ);
                    end
                    if a0[k] >= 0 || t>=(T-1)
                        vopt          = vstars[N_s];
                        sopt          = s_grid[N_s];
                        bopt          = bstars[N_s];
                    else
                        #println("bstar: $bstard vstar: $vstard")
                        sopt, vopt    = maxs(s_grid, vstars)
                        bopt          = LinearInterpolation(s_grid, bstars)(sopt);
                    end

                    S[t, k, j, i] = sopt;
                    V[t, k, j, i] = vopt;
                    B[t, k, j, i] = bopt;

                    # Budget constraint
                    # a_{t + 1} = 1 / q() * (a * s + y - ϕ(s) - c) + (1 - s) * a =>
                    # c= y + a * s - ϕ(s) - q(a(),a') * b
                    Ap[t, k, j, i] = bopt + (1.0 - sopt)*a0[k];
                    if Ap[t, k, j, i] > at[N_k] || Ap[t, k, j, i] < at[1]
                        Aphr = Ap[t, k, j, i] > at[N_k] ? at[N_k] : at[1];
                        println("Outside of A grid at $t, $k, $j, $i : $Aphr and $bopt, $sopt");
                        #return V, B, C, S, Ap;
                    else
                    #    print("..");
                        C[t, k, j, i] = y + a0[k]*sopt - LinearInterpolation(at, qhr)(Ap[t, k, j, i])*bopt - ϕ(sopt, ϕ_1, ϕ_2);
                    end
                end
            #end
        end
    end
    return V, B, C, S, Ap
end

"""
Solve for implied Q, given S
"""
function equilibriumQ(Q0::Array{Float64,3}, S::Array{Float64,4} , Ap::Array{Float64,4}, asset_grid::Array{Float64,2}, r::Float64, ϵ_probs::Array{Float64,2}, z_probs::Array{Float64,3})
    #  Dimension is T, z, ap, b
    qimplied = similar(Q0);
    T = size(Q0)[1];
    N_z=size(Q0)[2];
    N_k=size(Q0)[3];
    N_ϵ = size(ϵ_probs)[2];
    for t = T:-1:1
        if t >= T - 2
            qimplied[t, :, :] .= 1.0 / (1.0 + r)
            continue
        end

        for i = 1:N_z
            for k = 1:N_k
                qimplied[t, i, k] = 0.0;
                for ii=1:N_z
                    for j =1:N_ϵ
                        app = Ap[t+1, k, j, ii];
                        app = app < asset_grid[t+2,1]   ? asset_grid[t+2,1] : app #force in bounds (this shouldn't bind)
                        app = app > asset_grid[t+2,N_k] ? asset_grid[t+2,N_k] : app
                        shr = S[t+1, k, j, ii];
                        qp  = LinearInterpolation(asset_grid[t+2,:],Q0[t+1,ii,:])(app);
                        qimplied[t, i, k] =  (shr + (1-shr)*qp )*ϵ_probs[t+1,j]*z_probs[t+1,i,ii] +
                            qimplied[t, i, k] ;
                    end
                end
                #double check that these are defined and w/in natural bounds:
                if isfinite(qimplied[t, i, k]) && qimplied[t, i, k] > 0.0
                    qimplied[t, i, k] =  qimplied[t, i, k];
                elseif sum(S[t+1,k,:,:]) >= N_ϵ*N_z-1e-2
                    println("infinite q, all 1 at $t,$i,$k");
                    Sj1 = S[t+1,k,1,:];
                    Sj2 = S[t+1,k,2,:];
                    Sj3 = S[t+1,k,3,:];
                    println("S is $Sj1, $Sj2, $Sj3");
                    qimplied[t, i, k] = 1.0 / (1.0 + r);
                else
                    qhr = qimplied[t, i, k]
                    println("infinite/non-positive q ($qhr), not 1 at $t,$i,$k");
                    Sj1 = S[t+1,k,1,:];
                    Sj2 = S[t+1,k,2,:];
                    Sj3 = S[t+1,k,3,:];
                    println("S is $Sj1, $Sj2, $Sj3");
                    qimplied[t, i, k] = 0.0;
                end
            end
        end
        qimplied[t,:,:] .= qimplied[t,:,:]./(1+r); # needs to be w/in the t loop otherwise for periods with T>=2 will double multiply by 1/(1+r)
    end
    return qimplied
end
