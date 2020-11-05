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
    return ϕ_1 *(1 -  exp(-ϕ_2*(s-1.0)))/(1.0-exp(ϕ_2))
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
function bobjective(b::Float64, s::Float64, y::Float64, qhr::Array{Float64,1},
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
Solve for optimal b, given s 
"""
@inline function max_bobjective(s::Float64, y::Float64, qhr::Array{Float64,1},
     at::Array{Float64,1}, a::Float64, EV_grid::Array{Float64,1}, ψ::Float64,
     β::Float64, ϕ_1_1::Float64, ϕ_1_2::Float64, γ::Int64)


     bmax = min(at[N_k] - (1.0 - s) * a - eps(Float32)*10, (y+a*s) *(1.0+r));
     bmin = at[1]   - (1.0 - s) * a + eps(Float32)*10;
                        
     # bstard[di] = max_gss(b_min, b_max, (s, i, y, self.phi, qhr, psi_t, asset_grid[t + 1], EV_grid))
     # default routine is Brent's Method
     res = optimize(b -> -1*bobjective(b, s, y, qhr, at, a, EV_grid, ψ, β, ϕ_1_1, ϕ_1_2, γ), bmin, bmax);
     return -1.0*Optim.minimum(res);  # the maximum
end

"""
Find s* by maximizing the interpolated v*(s),
"""
@inline function maxs_interp(s_grid::Array{Float64,1}, vstars::Array{Float64,1})
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
    γ::Int64, ϕ_1::Float64, ϕ_2::Float64, r::Float64, εprobs::Array{Float64,2},
    zprobs::Array{Float64,3}, Ygrid::Array{Float64,3}, s_grid::Array{Float64,1},
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

        # 
        @inbounds Threads.@threads for zi = 1:N_z #
        #@inbounds Threads.@threads for idx = 1:(N_z*N_ε*N_k)
        #for idx = 1:(N_z*N_ε*N_k)
        #    ai = floor(Int,mod(idx-1,N_z*N_ε*N_k)/N_ε/N_z )+1;
        #    zi = floor(Int,mod(idx-1,N_z*N_ε)/N_ε )+1;
        #    εi =           mod(idx-1,N_ε)+1;
        #    fidx = (ai-1)*N_ε*N_k + (zi-1)*N_ε + εi 
        #    println("ai = $ai εi=$εi zi=$zi, filled in idx is $fidx")
        #end 
            qhr    = Q[t, zi,:];
            @inbounds for εi = 1:N_ε
                y = Ygrid[t, εi, zi] + reb;
                # Income transition probabilities
                tmat = εprobs[t,:]*zprobs[t,zi,:]'; # N_ε x N_z
                
                @inbounds for ai = 1:N_k #loop over a0

                    
                    a0hr = a0[ai];
                    # Right before certain death, the agent chooses to consume all that they can
                    if t == T
                        S[t, ai, εi, zi] = 1.0;
                        B[t, ai, εi, zi] = 0.0;
                        C[t, ai, εi, zi] = a0hr + y + reb;
                        V[t, ai, εi, zi] = u(C[t, ai, εi, zi], γ);
                        Ap[t, ai, εi, zi] = 0.0;
                        continue
                    end

                    # Expected value grid next period, given income states and assets
                    EV_grid = [dot(V[t + 1, aii, :, :], tmat) for aii = 1:N_k];
                    
                    apmin = at[1];
                    apmax = at[N_k];
                    #=
                    bstars = zeros(N_s);
                    vstars = zeros(N_s);

                    if a0hr >= 0 || t>=(T-1)
                        srange_a0 = N_s;
                    else
                        if ai==1
                            srange_a0 = 1;
                        else #exploit monotonicity
                            srange_a0 = slower[ai-1];
                    end
                    @inbounds for si = srange_a0:N_s

                        s    = s_grid[si];

                        # b cannot put a' outside bounds and can't make c negative:
                        bmax = min(apmax - (1.0 - s) * a0hr - eps(Float32)*10, (y+a0hr*s) *(1.0+r));
                        bmin = apmin   - (1.0 - s) * a0hr + eps(Float32)*10;
                        # c = y + a0[ai]*s - bmin /(1.0 + r);  # check this is positve?

                        # Note: EV interpolation should take place over assets_{t+1} grid
                        bstars[si],vstars[si]   = max_bobjective(s, y, qhr, at, a0hr, EV_grid, ψ, β, ϕ_1 , ϕ_2, γ);
                    
                    end
                    =#
                    if a0hr >= 0 || t>=(T-1)
                        sopt          = s_grid[N_s];
                        #bstars[si],vstars[si]   = max_bobjective(sopt, y, qhr, at, a0hr, EV_grid, ψ, β, ϕ_1 , ϕ_2, γ);
                        bmax = min(apmax - (1.0 - sopt) * a0hr - eps(Float32)*10, (y+a0hr*sopt) *(1.0+r));
                        bmin = apmin   - (1.0 - sopt) * a0hr + eps(Float32)*10;                       
                        b_solution_struct = optimize(b -> -1*bobjective(b, sopt, y, qhr, at, a0hr, EV_grid, ψ, β, ϕ_1, ϕ_2, γ), bmin, bmax);
                        vopt          = -1.0 * Optim.minimum(b_solution_struct);
                        bopt          = Optim.minimizer(b_solution_struct);
                    else
                        
                        res_s = optimize( shr-> -1.0*max_bobjective(shr, y, qhr, at, a0hr, EV_grid, ψ, β, ϕ_1 , ϕ_2, γ), 0.0,1.0, 
                            Brent(); iterations = 20  );
                        sopt    = Optim.minimizer(res_s);
                        vopt    = -1.0*Optim.minimum(res_s);
                        
                        bmax = min(apmax - (1.0 - sopt) * a0hr - eps(Float32)*10, (y+a0hr*sopt) *(1.0+r));
                        bmin = apmin   - (1.0 - sopt) * a0hr + eps(Float32)*10;                       
                        b_solution_struct = optimize(b -> -1*bobjective(b, sopt, y, qhr, at, a0hr, EV_grid, ψ, β, ϕ_1, ϕ_2, γ), bmin, bmax);
                        bopt          = Optim.minimizer(b_solution_struct);
                        
                    end
                    
                    S[t, ai, εi, zi] = sopt;
                    V[t, ai, εi, zi] = vopt;
                    B[t, ai, εi, zi] = bopt;

                    # Budget constraint
                    # a_{t + 1} = 1 / q() * (a * s + y - c) + (1 - s) * a =>
                    # c= y + a * s - q(a(),a' * b
                    Ap[t, ai, εi, zi] = bopt + (1.0 - sopt)*a0[ai];
                    if Ap[t, ai, εi, zi] > at[N_k] || Ap[t, ai, εi, zi] < at[1]
                        Aphr = Ap[t, ai, εi, zi] > at[N_k] ? at[N_k] : at[1];
                        println("Outside of A grid at $t, $ai, $εi, $zi : $Aphr and $bopt, $sopt");
                        #return V, B, C, S, Ap;
                    else
                        C[t, ai, εi, zi] = y + a0[ai]*sopt - LinearInterpolation(at, qhr)(Ap[t, ai, εi, zi])*bopt - ϕ(sopt, ϕ_1, ϕ_2);
                    end
                end
            end
        end
    end
    return V, B, C, S, Ap
end

"""
Solve for implied Q, given S
"""
function equilibriumQ(Q0::Array{Float64,3}, S::Array{Float64,4} , Ap::Array{Float64,4}, asset_grid::Array{Float64,2}, r::Float64, εprobs::Array{Float64,2}, zprobs::Array{Float64,3})
    #  Dimension is T, z, ap, b
    qimplied = similar(Q0);
    T = size(Q0)[1];
    N_z=size(Q0)[2];
    N_k=size(Q0)[3];
    N_ϵ = size(εprobs)[2];
    for t = T:-1:1
        if t >= T - 2
            qimplied[t, :, :] .= 1.0 / (1.0 + r)
            continue
        end

        for zi = 1:N_z
            for ai = 1:N_k
                qimplied[t, zi, ai] = 0.0;
                for zii=1:N_z
                    for εi =1:N_ϵ
                        app = Ap[t+1, ai, εi, zii];
                        app = app < asset_grid[t+2,1]   ? asset_grid[t+2,1] : app #force in bounds (this shouldn't bind)
                        app = app > asset_grid[t+2,N_k] ? asset_grid[t+2,N_k] : app
                        shr = S[t+1, ai, εi, zii];
                        qp  = LinearInterpolation(asset_grid[t+2,:],Q0[t+1,zii,:])(app);
                        qimplied[t, zi, ai] =  (shr + (1-shr)*qp )*εprobs[t+1,εi]*zprobs[t+1,zi,zii] +
                            qimplied[t, zi, ai] ;
                    end
                end
                #double check that these are defined and w/in natural bounds:
                if isfinite(qimplied[t, zi, ai]) && qimplied[t, zi, ai] > 0.0
                    qimplied[t, zi, ai] =  qimplied[t, zi, ai];
                elseif sum(S[t+1,ai,:,:]) >= N_ϵ*N_z-1e-2
                    println("infinite q, all 1 at $t,$i,$ai");
                    Sj1 = S[t+1,ai,1,:];
                    Sj2 = S[t+1,ai,2,:];
                    Sj3 = S[t+1,ai,3,:];
                    println("S is $Sj1, $Sj2, $Sj3");
                    qimplied[t, zi, ai] = 1.0 / (1.0 + r);
                else
                    qhr = qimplied[t, zi, ai]
                    println("infinite/non-positive q ($qhr), not 1 at $t,$zi,$ai");
                    Sj1 = S[t+1,ai,1,:];
                    Sj2 = S[t+1,ai,2,:];
                    Sj3 = S[t+1,ai,3,:];
                    println("S is $Sj1, $Sj2, $Sj3");
                    qimplied[t, zi, ai] = 0.0;
                end
            end
        end
        qimplied[t,:,:] .= qimplied[t,:,:]./(1+r); # needs to be w/in the t loop otherwise for periods with T>=2 will double multiply by 1/(1+r)
    end
    return qimplied
end


function steadystateDist!( Aεz_dist::Array{Float64,4}, Ap::Array{Float64,4}, survivalprobs::Vector{Float64}, T, N_k, N_s, N_z, N_ε, asset_grid::Array{Float64,2}, εprobs::Array{Float64,2}, zprobs::Array{Float64,3})

    survivalprobs_cml = sum(1.0 .- survivalprobs[1:(T-1)]) #fraction dead before period T+1
    survivalprobs[T] = 0.0;
    Nstate = (T-1)*N_k*N_ε*N_z;
    a0idx = findlast(asset_grid[2,:] .<= 0.0)
    # state is (t-1)*N_k*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi
    flatstatemat = spzeros(Float64,Nstate,Nstate);
    for t=1:(T-1)
        for ai=1:N_k
            for ei=1:N_ε
                for zi=1:N_z
                    apL   = Ap[t,ai,ei,zi]> asset_grid[t+1,1] ? findlast( Ap[t,ai,ei,zi] .< asset_grid[t+1,:]) : 1;
                    apLwt =  apL < N_k ? (asset_grid[t+1,apL+1] - Ap[t,ai,ei,zi] )/(asset_grid[t+1,apL+1] - asset_grid[t+1,apL]) : 1.0 ;
                    if t<(T-1)
                        for zzi=1:N_z
                            for eei=1:N_ε
                                flatstatemat[ (t-1)*N_k*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi, (t)*N_k*N_z*N_ε + (apL-1)*N_ε*N_z + (eei-1)*N_z + zzi] =
                                    apLwt*zprobs[t,zi,zzi]*εprobs[t,eei] * survivalprobs[t];
                                if apL<N_k
                                    flatstatemat[ (t-1)*N_k*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi, (t)*N_k*N_z*N_ε + apL*N_ε*N_z + (eei-1)*N_z + zzi] =
                                        (1.0 - apLwt)*zprobs[t,zi,zzi]*εprobs[t,eei] * survivalprobs[t];;
                                end
                            end
                        end
                    else
                        #everyone dies and goes to assets = 0
                        flatstatemat[ (t-1)*N_k*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi, 0*N_k*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi] = 1.0; 
                    end    
                end
            end
        end
    end

    Npara = (T-1)*N_k*N_ε*N_z;
    zerocol = Int64[];
    for ipara =1:Npara
        statemat_rowsum = sum(flatstatemat[:,ipara]);
        if statemat_rowsum <= 0.0
            push!(zerocol,ipara);
        end
    end 
    # eigen value for steady state dist:
    nzerocol = size(zerocol);
    println("Computing steady state via eigen-values, transition matrix has $nzerocol columns never visited:")
    # This bit seems not to work with the full matrix, but does with smaller ones. Size limit?
    #numeigv = 3; # just to check
    #flatstateeigs, flatstatevecs = eigs(flatstatemat; nev=numeigv, ncv=max(20,2*numeigv+1), which=:LR );
    #println( "The eigen values are $flatstateeigs" )
    flatstatevecs = (flatstatemat^(Npara) )';
    
    #unroll it into the matrix
    for t=1:(T-1)
        mass_t = 0.0;
        for ai=1:N_k
            for ei=1:N_ε
                for zi=1:N_z
                    Aεz_dist[t,ai,ei,zi] = flatstatevecs[(t-1)*N_k*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi , 1];
                    mass_t += Aεz_dist[t,ai,ei,zi];
                end
            end
        end
        Aεz_dist[t,:,:,:] .= Aεz_dist[t,:,:,:] ./ mass_t;
    end

end
