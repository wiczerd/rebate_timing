"""
Return the asset grids, given the total (labor) income grid, T, and N_k
α = 1 for uniformly spaced gripoints.
"""
function assetGrid(T::Int64, income_grid::Array{Float64,3}, N_k::Int64,Q::Array{Float64}, asset_grid_old::Array{Float64,2},  α::Int64 = 1)

    asset_max  = zeros(T+1);
    asset_min  = zeros(T+1);
    asset_grid = zeros(T+1, N_k); # T + 1 = year after death, define for convenience in backwardSolve

    # Minimum assets (Maximum debt)
    for t = T:-1:2
        qhr = Q[t,1,:];
        income_min = minimum(income_grid[t,:,:]);
        if t<T 
            
            qbar = max(LinearInterpolation(asset_grid_old[t+1,:], qhr; extrapolation_bc=Line())(asset_min[t+1] ) , 1e-3);
        else 
            qbar = 1.0/(1.0 +r );
        end
        if t>T_retire
            #spps pay back the whole thing:
            asset_min[t] = asset_min[t+1]*qbar - income_min + 1e-4;
        else
            #spps not pay back at all
            #asset_min[t] = asset_min[t+1] - (income_min + 1e-4)/qbar;
            asset_min[t] = asset_min[t+1]*qbar - income_min + 1e-4;
        end
        if t==T 
            asset_min[t]=0
        end
    end

    # Make k_max the largest amount the consumer could ever save, even
    # if they got the highest income and saved half of income every quarter
    for t = 1:T
        income_max     = maximum(income_grid[t,:,:]);
        asset_max[t+1] = min(asset_max[t]/Q[t,1,N_k] + 0.95*income_max , 100)
    end

    for t = 1:T
        N_p = ceil(Int64, N_k*0.70);
        N_n = N_k - N_p;

        assets_p  = collect(LinRange(0.0, asset_max[t]^(1/α), N_p)).^α;  # optional convex grid for positive assets
        assets_n  = collect(LinRange(asset_min[t], -1.0/N_n, N_n));             # uniformly spaced points for negative assets

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
@inline function u(c::Float64,cbar, γ::Int64)
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
function bobjective(b::Float64, s::Float64, y::Float64, q_spn::Array{Float64,1}, #Schumaker, #
    asset_grid::Array{Float64,1}, a::Float64, EV_spn:: Schumaker, #Array{Float64,1}, #
    cbar::Float64, ψ::Float64, β::Float64, ϕ_1::Float64, ϕ_2::Float64, γ::Int64)
    # Note: b, s, y, ψ are levels. a0idx is the index of a0.
    # qhr, asset_grid, EV are an arrays
    # Get level of capital implied by next-period consumption choice
    ap = b + (1.0 - s)*a;
    if ap>asset_grid[length(asset_grid)] || ap<asset_grid[1]
        println("ap out of bounds: (ap,b,s,a)= ($ap, $b, $s, $a)");
    end
    
    EV = SchumakerSpline.evaluate(EV_spn,ap) ; #LinearInterpolation(asset_grid, EV_spn)(ap); #
    #ap_step = abs(ap)*0.05;
    #if  ap + ap_step < asset_grid[N_k]  && ap - ap_step > asset_grid[1]
    #    EV_ap_noise = EV/3.0 + ( max( LinearInterpolation(asset_grid, EV_spn)(ap+ap_step) , 0.0) + max( LinearInterpolation(asset_grid, EV_spn)(ap-ap_step) , 0.0) ) / 3.0;
    #else 
    #    EV_ap_noise = EV;
    #end 
    
    q_ap = max( LinearInterpolation(asset_grid, q_spn)(ap) , 0.0);  #max( SchumakerSpline.evaluate(q_spn,ap), 0.0); #
    
    c = y + a*s - q_ap*b;
    util_penalty = 0.0;
    if c<0.0
        util_penalty = c^2;
        c=1e-4;
        instutil = u(c,cbar, γ) ;
        #println("negative c with b as $b and $instutil");
    end
    return u(c,cbar, γ) - util_penalty - ϕ(s, ϕ_1, ϕ_2) + β*ψ*EV
end

"""
Solve for optimal b, given s 
"""

function max_bobjective(s::Float64, y::Float64, q_spn::Array{Float64,1}, #Schumaker, 
    at::Array{Float64,1}, a::Float64, EV_spn::Schumaker, #Array{Float64,1}, #
    cbar::Float64,ψ::Float64, β::Float64, ϕ_1_1::Float64, ϕ_1_2::Float64, γ::Int64)

    bmin = at[1]   - (1.0 - s) * a + eps(Float32)*10;
    bmax = at[N_k] - (1.0 - s) * a - eps(Float32)*10;
    
    b_c0 = max( (y + a*s-1e-4)*(1.0+r), at[1] - (1.0 - s)*a ) #the b with c close to 0
    qmin =  b_c0 + (1.0 - s)*a <at[N_k] ? max( LinearInterpolation(at, q_spn)(b_c0 + (1.0 - s)*a) , 0.0) : 1.0/(1.0 + r);
      #SchumakerSpline.evaluate(q_spn, b_c0 + (1.0 - s)*a )
    b_c0 = (y + a*s-1e-4)/qmin  #the b with c close to 0
    
    bmax = min(bmax , b_c0);
    
    if bmax > bmin 
        res = optimize(b -> -1.0*bobjective(b, s, y, q_spn, at, a, EV_spn, cbar,ψ, β, ϕ_1_1, ϕ_1_2, γ), bmin, bmax);
        bopt = Optim.minimizer(res);
        return -1.0*Optim.minimum(res), bopt;  # the maximum
    else 
        println("bmax < bmin ($bmax > $bmin) in max_bojb  at a = $a , s = $s")
        bopt = bmin;
        return bobjective(bmin, s, y, q_spn, at, a, EV_spn, cbar, ψ, β, ϕ_1_1, ϕ_1_2, γ), bopt
    end 
end

"""
Find s* by maximizing the interpolated v*(s),
"""
function maxs_interp(s_grid::Array{Float64,1}, vstars::Array{Float64,1})
    # maxd = lambda dhr : -1.0*np.interp( np.subtract(np.ones(N_delinq), self.s_grid), vstard, dhr)
    # sopt, vopt  = opt.minimize(maxd, brack=(0, 1))
    #itp = LinearInterpolation(ones(N_delinq)-s_grid, vstard; extrapolation_bc=Throw());
    vspn = Schumaker(s_grid, vstars); #LinearInterpolation(s_grid, vstars);
    res = optimize(s-> -1.0*SchumakerSpline.evaluate(vspn,s), 0, 1); # default routine is Brent's Method
    #vspn = LinearInterpolation(s_grid, vstars);
    #res = optimize(s-> -1*vspn(s), 0, 1); # default routine is Brent's Method
    return Optim.minimizer(res), -1*Optim.minimum(res) ; 
end


"""
Taking Q as given, solve for V,B,C,S,K
"""
function backwardSolve!(ms::sol, mod::model, 
    T::Int64, N_k::Int64, N_z::Int64, N_ε::Int64, Q::Array{Float64,3})

    for t =  T:(-1):1
        
        # a' grid
        at  = mod.asset_grid[t + 1, :];
        a0  = mod.asset_grid[t, :];

        # Rebate
        reb = mod.rebate*(t == mod.t_rebate);

        # Exogenous survival prob for the period
        ψ = survivalprobs[t];

        # @inbounds Threads.@threads 
        @inbounds for zi = 1:N_z 
        #@inbounds Threads.@threads for idx = 1:(N_z*N_ε*N_k)
        #for idx = 1:(N_z*N_ε*N_k)
        #    ai = floor(Int,mod(idx-1,N_z*N_ε*N_k)/N_ε/N_z )+1;
        #    zi = floor(Int,mod(idx-1,N_z*N_ε)/N_ε )+1;
        #    εi =           mod(idx-1,N_ε)+1;
        #    fidx = (ai-1)*N_ε*N_k + (zi-1)*N_ε + εi 

        
            qhr    = Q[t, zi,:];

            # Expected value grid next period, given income states and assets
            EV_grid = zeros(N_k); #[dot(V[t + 1, aii, :, :], tmat) for aii = 1:N_k];
            for aii = 1:N_k
                for eii = 1:N_ε
                    for zii = 1:N_z 
                        EV_grid[aii] = ms.V[t+1, aii, eii, zii ]*mod.εprobs[t,eii]*mod.zprobs[t,zi,zii] + EV_grid[aii];
                    end 
                end
            end
            EV_spn = Schumaker(at,EV_grid); #Spline1D(at,EV_grid; w=ones(N_k), k=3, bc="nearest");


            @inbounds Threads.@threads for εi = 1:N_ε
                y = mod.Ygrid[t, εi, zi] + reb;
                
                slower = 0.0;
                @inbounds for ai = 1:N_k #loop over a0

                    a0hr = a0[ai];

                    # Right before certain death, the agent chooses to consume all that they can
                    if t == T
                        ms.S[t, ai, εi, zi] = 1.0;
                        ms.B[t, ai, εi, zi] = 0.0;
                        ms.C[t, ai, εi, zi] = a0hr + y + reb;
                        ms.V[t, ai, εi, zi] = u(ms.C[t, ai, εi, zi], mod.cbar, γ);
                        ms.Ap[t, ai, εi, zi] = 0.0;
                        continue
                    end
                    
                    apmin = at[1];
                    apmax = at[N_k];
                    
                    #=
                    # Coment from here for 2-dim Brent Method. Below is grid search
                    vstars = zeros(N_s);
                    
                    srange_a0 = 1;
                    srange_aN = N_s;
                    maxfeasibles= N_s;
                    #if ai==1
                        srange_a0 = 1;
                    #else #exploit monotonicity
                    #    srange_a0 = slower[ai-1];
                    #end
                    
                    if a0hr >= 0 || t>=(T-1)
                        sopt = 1.0;
                        vopt = max_bobjective(sopt, y, q_spn, at, a0hr, EV_spn, ψ, β, mod.ϕ_1 , mod.ϕ_2, γ);;
                    #elseif y - SchumakerSpline.evaluate(q_spn, at[N_k] )*(at[1] - a0hr) < 0.0 #can't get to positive consmption in this situation

                    #    sopt = 0.0;
                    #    vopt = max_bobjective(sopt, y, q_spn, at, a0hr, EV_spn, ψ, β, mod.ϕ_1 , mod.ϕ_2, γ);
                    else
                        @inbounds for si = srange_a0:srange_aN

                            shr    = s_grid[si];
    
                            # b cannot put a' outside bounds and can't make c negative:
                            # bmax = min(apmax - (1.0 - shr) * a0hr - eps(Float32)*10, (y+a0hr*shr) *(1.0+r));
                            # bmin = apmin   - (1.0 - shr) * a0hr + eps(Float32)*10;
                            # c = y + a0[ai]*s - bmin /(1.0 + r);  # check this is positve?
    
                            # Note: EV interpolation should take place over assets_{t+1} grid
                            vstars_hr   = max_bobjective(shr, y, q_spn, at, a0hr, EV_spn, ψ, β, mod.ϕ_1 , mod.ϕ_2, γ);
                            if isfinite( vstars_hr)
                                vstars[si] = vstars_hr;
                            else 
                                maxfeasibles = si-1;
                                break;
                            end
                            
                        end
                        #println("at $εi , $ai, maxfeasibles: $maxfeasibles")
                        if maxfeasibles > 0 
                            #sopt, vopt = maxs_interp(s_grid[srange_a0:maxfeasibles], vstars[srange_a0:maxfeasibles]);
                            vopt, sopti = findmax(vstars);
                            sopt = s_grid[sopti];
                        else 
                            sopt = 0.0
                            vopt = max_bobjective(sopt, y, q_spn, at, a0hr, EV_spn, ψ, β, mod.ϕ_1 , mod.ϕ_2, γ)
                            println("vopt with no options for s is $vopt at $zi, $εi, $ai ")
                        end 
                    end 
                    bmax = min(apmax - (1.0 - sopt) * a0hr - eps(Float32)*10, (y+a0hr*sopt) *(1.0+r));
                    bmin = apmin   - (1.0 - sopt) * a0hr + eps(Float32)*10;                       
                        
                    b_solution_struct = optimize(b -> -1*bobjective(b, sopt, y, q_spn, at, a0hr, EV_spn, ψ, β, mod.ϕ_1, mod.ϕ_2, γ), bmin, bmax);
                    bopt          = Optim.minimizer(b_solution_struct);
                    =#
                    # Comment from here for grid search. Below is Brent's method
                    
                    if  a0hr >= 0 || t>=T-1
                        sopt          = 1.0;
                        
                        bmax = min(apmax - (1.0 - sopt) * a0hr - eps(Float32)*10, (y+a0hr*sopt )*(1.0+r) );
                        bmin = apmin   - (1.0 - sopt) * a0hr + eps(Float32)*10  ;  
                        if bmin > bmax
                            println("bmin > bmax: $bmin > $bmax at $t, $ai, and $a0hr")
                            bopt = bmin;
                            vopt = bobjective(bopt, sopt, y, qhr, at, a0hr, EV_spn, mod.cbar,ψ, β, mod.ϕ_1, mod.ϕ_2, γ);
                            println("bmin > bmax: $bmin > $bmax with $vopt, $bopt")
                            
                        else 
                            b_solution_struct = optimize(b -> -1*bobjective(b, sopt, y, qhr, at, a0hr, EV_spn, mod.cbar,ψ, β, mod.ϕ_1, mod.ϕ_2, γ), bmin, bmax);
                            vopt          = -1.0 * Optim.minimum(b_solution_struct);
                            bopt          = Optim.minimizer(b_solution_struct);
                        end 
                    else
                        smax = 1.0;

                        if y+ a0hr - qhr[1] *( at[1]-a0hr ) <= 0 #consumption in this case
                            #have to default a little bit:
                            smax = min(1.0,  (1e-6 - y + qhr[1]*at[1]-qhr[1]*a0hr)/(a0hr - a0hr*qhr[1])  )
                        end
                        
                        res_s = optimize( shr-> -1.0*(max_bobjective(shr, y, qhr, at, a0hr, EV_spn, mod.cbar,ψ, β, mod.ϕ_1 , mod.ϕ_2, γ)[1]), 0.0,smax, 
                            Brent(); iterations = 100,abs_tol=1e-5  );
                        sopt    = Optim.minimizer(res_s);
                        
                        vopt,bopt = max_bobjective(sopt, y, qhr, at, a0hr, EV_spn, mod.cbar,ψ, β, mod.ϕ_1 , mod.ϕ_2, γ);
                        slower = sopt;
                        
                    end
                    

                    ms.S[t, ai, εi, zi] = sopt;
                    ms.V[t, ai, εi, zi] = vopt;
                    ms.B[t, ai, εi, zi] = bopt;

                    # Budget constraint
                    # a_{t + 1} = 1 / q() * (a * s + y - c) + (1 - s) * a =>
                    # c= y + a * s - q(a(),a' * b
                    ms.Ap[t, ai, εi, zi] = bopt + (1.0 - sopt)*a0[ai];
                    if ms.Ap[t, ai, εi, zi] > at[N_k] || ms.Ap[t, ai, εi, zi] < at[1]
                        Aphr = ms.Ap[t, ai, εi, zi] #> at[N_k] ? at[N_k] : at[1];
                        at1 = at[1];
                        println("Outside of A grid at $t, $ai, $εi, $zi : $Aphr < $at1 and $bopt, $sopt");
                        #return V, B, C, S, Ap;
                    else
                        ms.C[t, ai, εi, zi] = y + a0[ai]*sopt - LinearInterpolation(at, qhr)(ms.Ap[t, ai, εi, zi])*bopt ;
                    end
                end
            end
        end
    end
    #return V, B, C, S, Ap
    return 1
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
                for εi =1:N_ϵ
                    for zii=1:N_z
                
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
                    qimplied[t, zi, ai] = 1.0 / (1.0 + r);
                else
                    qhr = qimplied[t, zi, ai]
                    println("infinite/non-positive q ($qhr), not 1 at $t,$zi,$ai");
                    qimplied[t, zi, ai] = 0.0;
                end
            end
        end
        qimplied[t,:,:] .= qimplied[t,:,:]./(1+r); # needs to be w/in the t loop otherwise for periods with T>=2 will double multiply by 1/(1+r)
    end
    return qimplied
end

function draw_shocks!(ht::hists,mod::model, seed::Int64=12281951)
    Random.seed!(seed)

    ht.a0qtlhist .= rand(Nsim);
    
    ageproto = rand(Nsim);
    for i=1:Nsim
        ht.agehist[i,1] = ceil(Int64,ageproto[i]*T);
        for t=2:(Tsim+Tburnin)
            agp1 = ht.agehist[i,t-1] + 1
            ht.agehist[i,t] = agp1<T ? agp1 : 1
        end 
    end
    
    εproto =  rand(Nsim,Tsim+Tburnin);
    zproto =  rand(Nsim,Tsim+Tburnin);
    cumεprobs = zeros(T,N_ε+1);
    cumzprobs = zeros(T,N_z,N_z+1);
    cumzergod = zeros(T,N_z+1);
    for t=1:T 
        cumεprobs[t,2:(N_ε+1)] .= cumsum(mod.εprobs[t,:]);
        for zi=1:N_z
            cumzprobs[t,zi,2:(N_z+1)] .= cumsum(mod.zprobs[t,zi,:]);
        end 
        zergod_t = mod.zprobs[t,:,:]^40;
        cumzergod[t,2:(N_z+1)] .= cumsum(zergod_t[1,:]);
    end 
    
    for i=1:Nsim
        ht.zidxhist[i,1] = sum( zproto[i,1].>cumzergod[ ht.agehist[i,1],: ] );
        for t=1:(Tsim+Tburnin)
            if t>1
                izlast = ht.zidxhist[i,t-1];
                ht.zidxhist[i,t] = sum( zproto[i,t] .>= cumzprobs[ht.agehist[i,t], izlast,: ]  );
            end
            ht.εidxhist[i,t] = sum( εproto[i,t] .>= cumεprobs[ht.agehist[i,t], : ]  );
        end 
    end
end 

function sim_hists!(ht::hists, ms::sol, mod::model, T, N_k,N_z, N_ε)

    
    asset_alpha = ones(T,N_ε, N_z); asset_beta=ones(T,N_ε, N_z);  #will iteratively fit a beta distribution 
    for dist_iter = 1:5
        
        for i=1:Nsim 

            for it=1:(Tsim+Tburnin)
                ei_hr = ht.εidxhist[i,it];
                zi_hr = ht.zidxhist[i,it];
                age_hr=  ht.agehist[i,it];

                if it==1 # Initialize

                    if age_hr >1 && asset_alpha[age_hr,ei_hr,zi_hr] > 0
                        ht.ahist[i,1] = quantile(Beta(asset_alpha[age_hr,ei_hr,zi_hr],asset_beta[age_hr,ei_hr,zi_hr]),ht.a0qtlhist[i] )*(mod.asset_grid[age_hr,N_k] - mod.asset_grid[age_hr,1]) + mod.asset_grid[age_hr,1];
                    elseif age_hr==1 #young all go to the same place
                        ht.ahist[i,1] = 0.0;
                    elseif asset_alpha[age_hr,ei_hr,zi_hr] <=0 #if there was a degenerate distribution I stored the average as -asset_alpha
                        ht.ahist[i,1] = -asset_alpha[age_hr,ei_hr,zi_hr]*(mod.asset_grid[age_hr,N_k] - mod.asset_grid[age_hr,1]) + mod.asset_grid[age_hr,1];
                    end 
                end     
                aL = 1; aLwt = 0.0;
                if it<(Tsim+Tburnin) 
                    if age_hr > 1
                        aL = sum( ht.ahist[i,it] .>=  mod.asset_grid[age_hr,:] );
                        aL = aL>= N_k ? N_k-1 : aL;
                        aLwt =  (mod.asset_grid[age_hr,aL+1] - ht.ahist[i,it])/(mod.asset_grid[age_hr,aL+1] - mod.asset_grid[age_hr,aL])
                    else #for first period, doesn't matter
                        aL = 1; aLwt = 0.0;
                    end 
                    ht.ahist[i,it+1] = aLwt*ms.Ap[age_hr,aL,ei_hr,zi_hr ] +
                        (1.0 - aLwt)*ms.Ap[ age_hr,aL+1,ei_hr,zi_hr ];
                end
                ht.chist[i,it] = aLwt*ms.C[age_hr,aL,ei_hr,zi_hr] + 
                    (1.0 - aLwt)*ms.C[age_hr,aL+1,ei_hr,zi_hr];
                ht.shist[i,it] = aLwt*ms.S[age_hr,aL,ei_hr,zi_hr] + 
                    (1.0 - aLwt)*ms.S[age_hr,aL+1,ei_hr,zi_hr];
                ht.bhist[i,it] = aLwt*ms.B[age_hr,aL,ei_hr,zi_hr] + 
                    (1.0 - aLwt)*ms.B[age_hr,aL+1,ei_hr,zi_hr];
            end
        end

        # update asset_alpha, asset_beta
        for t=1:T
            for ei=1:N_ε
                for zi=1:N_z
                    mask_hr = (ht.εidxhist.==ei) .& (ht.zidxhist.==zi) .& (ht.agehist .== t) ;
                    ahist_hr = (ht.ahist[mask_hr] .- mod.asset_grid[t,1])./(mod.asset_grid[t,N_k] - mod.asset_grid[t,1]);
                    Ea_hr = 0.5; Va_hr = 0.0;
                    if(sum(mask_hr)>0)
                        Ea_hr = mean(ahist_hr);
                        Va_hr = var(ahist_hr);
                    end
                    #if Ea_hr<=0 || Va_hr<0 || isnan(Ea_hr) || isnan(Va_hr)
                    #    println(" Ea_hr,Va_hr, t = $Ea_hr, $Va_hr, $t, $ei, $zi ");
                    #end 
                    if Va_hr>0 && !isnan(Va_hr)
                        asset_alpha[t,ei,zi] =      Ea_hr *(Ea_hr*(1.0-Ea_hr)/Va_hr-1.0);
                        asset_beta[t,ei,zi]  = (1.0-Ea_hr)*(Ea_hr*(1.0-Ea_hr)/Va_hr-1.0);
                    else #in the case where everyone goes to the same place:
                        asset_alpha[t,ei,zi] = -Ea_hr;
                        asset_beta[t,ei,zi]  = -Ea_hr;
                    end

                
                end
            end
        end
    end

    Ea_tot = mean(ht.ahist);
    Es_tot = mean(ht.shist);

    println("The average asset position is $Ea_tot and pay-back is $Es_tot");


end


function steadystateDist!( ms::sol, mod::model, survivalprobs::Vector{Float64}, T, N_k, N_z, N_ε)

    survivalprobs_cml = sum(1.0 .- survivalprobs[1:(T-1)]) #fraction dead before period T+1
    survivalprobs[T] = 0.0;
    N_krnd = Int64(N_k/5);
    asteps = Int64(N_k/N_krnd);
    asset_grid_rnd = zeros(T,N_krnd);
    asset_space    = zeros(T,N_krnd);
    for t=1:(T-1)
        for airnd=1:N_krnd
            for ai=1:asteps
                if (airnd-1)*asteps + ai+1 < N_k 
                    asset_space[t,airnd]+= mod.asset_grid[t,(airnd-1)*asteps + ai+1] - mod.asset_grid[t,(airnd-1)*asteps + ai];
                end
                asset_grid_rnd[t,airnd] += mod.asset_grid[t,(airnd-1)*asteps + ai]/asteps;
            end
        end
        # need to make sure the ends are the same
        asset_grid_rnd[t,1] = mod.asset_grid[t,1]
        asset_grid_rnd[t,N_krnd] = mod.asset_grid[t,N_k]
    end
    Nstate = (T-1)*N_krnd*N_ε*N_z;
    a0idx = findlast(  asset_grid_rnd[2,:] .<= 0.0 );  #mod.asset_grid[2,:] .<= 0.0)
    # state is (t-1)*N_k*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi
    flatstatemat = spzeros(Float64,Nstate,Nstate);
    for t=1:(T-1)
        for ai=1:N_krnd
            for ei=1:N_ε
                for zi=1:N_z
                    #apL   = ms.Ap[t,ai,ei,zi]>  mod.asset_grid[t+1,1] ? findlast( ms.Ap[t,ai,ei,zi] .>  mod.asset_grid[t+1,:]) : 1;
                    #apL   = isnan(apL) | isnothing(apL) ? N_k : apL;
                    #apLwt =  apL < N_k ? ( mod.asset_grid[t+1,apL+1] - ms.Ap[t,ai,ei,zi] )/( mod.asset_grid[t+1,apL+1] -  mod.asset_grid[t+1,apL]) : 1.0 ;
                    Ap_rnd = 0.0;
                    for aihr=1:asteps
                        Ap_rnd = ms.Ap[t,(ai-1)*asteps+aihr,ei,zi]/asteps + Ap_rnd;
                    end 
                    apL   = Ap_rnd >  asset_grid_rnd[t+1,1] ? sum( Ap_rnd .>=  asset_grid_rnd[t+1,:]) : 1;
                    apL   = isnan(apL) | isnothing(apL) ? N_krnd : apL;
                    apLwt =  apL < N_krnd ? ( asset_grid_rnd[t+1,apL+1] - Ap_rnd )/( asset_grid_rnd[t+1,apL+1] -  asset_grid_rnd[t+1,apL]) : 1.0 ;
                    if t<(T-1)
                        for zzi=1:N_z
                            for eei=1:N_ε
                                flatstatemat[ (t-1)*N_krnd*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi, (t)*N_krnd*N_z*N_ε + (apL-1)*N_ε*N_z + (eei-1)*N_z + zzi] =
                                    apLwt* mod.zprobs[t,zi,zzi]* mod.εprobs[t,eei] * survivalprobs[t] +
                                    flatstatemat[ (t-1)*N_krnd*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi, (t)*N_krnd*N_z*N_ε + (apL-1)*N_ε*N_z + (eei-1)*N_z + zzi];
                                if apL<N_krnd
                                    flatstatemat[ (t-1)*N_krnd*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi, (t)*N_krnd*N_z*N_ε + apL*N_ε*N_z + (eei-1)*N_z + zzi] =
                                        (1.0 - apLwt)*mod.zprobs[t,zi,zzi]* mod.εprobs[t,eei] * survivalprobs[t] +
                                        flatstatemat[ (t-1)*N_krnd*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi, (t)*N_krnd*N_z*N_ε + apL*N_ε*N_z + (eei-1)*N_z + zzi];
                                end
                            end
                        end
                    else
                        #everyone dies and goes to same assets in ag 0
                        flatstatemat[ (t-1)*N_krnd*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi, 0*N_krnd*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi] = 1.0; 
                    end    
                end
            end
        end
    end

    Npara = (T-1)*N_krnd*N_ε*N_z;
    zerocol = Int64[];
    for ipara =1:Npara
        statemat_rowsum = sum(flatstatemat[:,ipara]);
        if statemat_rowsum <= 0.0
            push!(zerocol,ipara);
        end
    end 

    # eigen value for steady state dist:
    nzerocol = size(zerocol);
    println("Computing steady state via eigen-value method, transition matrix has $nzerocol columns never visited")
    # This bit seems not to work with the full matrix, but does with smaller ones. Size limit?
    numeigv = 3; # just to check
    flatstateeigs, flatstatevecs = eigs(flatstatemat; nev=numeigv, ncv=max(20,2*numeigv+1), which=:LR );
    println( "The eigen values are $flatstateeigs" )
    
    
    #unroll it into the matrix
    for t=1:(T-1)
        mass_t = 0.0;
        for ai=1:N_krnd
            for ei=1:N_ε
                for zi=1:N_z
                    dist_hr = real(flatstatevecs[ (t-1)*N_krnd*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi ,1]);
                    mass_t += dist_hr*asset_space[t,ai];
                    for ai_step = 1:asteps 
                        ai_hr = (ai-1)*asteps + ai_step;
                        ms.Aεz_dist[t,ai_hr,ei,zi] = dist_hr;
                    end
                end
            end
        end
        ms.Aεz_dist[t,:,:,:] .= ms.Aεz_dist[t,:,:,:] ./ mass_t;
    end

end
