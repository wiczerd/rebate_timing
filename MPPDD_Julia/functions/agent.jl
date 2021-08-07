"""
Return the asset grids, given the total (labor) income grid, N_ages, and N_k
α = 1 for uniformly spaced gripoints.
"""
function assetGrid(N_ages::Int64, income_grid::Array{Float64,3}, N_k::Int64,Q::Array{Float64}, asset_grid_old::Array{Float64,2}, r::Float64, α::Int64 = 1)

    asset_max  = zeros(N_ages+1);
    asset_min  = zeros(N_ages+1);
    asset_grid = zeros(N_ages+1, N_k); # N_ages + 1 = year after death, define for convenience in backwardSolve

    # Minimum assets (Maximum debt)
    for t = N_ages:-1:1
        qhr = Q[t,1,:];
        income_min = minimum(income_grid[t,:,:]);
        if t<N_ages 
            qbar = max(LinearInterpolation(asset_grid_old[t+1,:], qhr; extrapolation_bc=Line())(asset_min[t+1] ) , 1e-2);
            if isnan(qbar)
                qbar = qhr[1];
            end
        else 
            qbar = 1.0/(1.0 +r );
        end
        if t>=age_retire
            #spps pay back the whole thing:
            asset_min[t] = asset_min[t+1]*qbar - income_min/2.0 + 1e-4;
        else
            #spps not pay back at all
            #asset_min[t] = asset_min[t+1] - (income_min + 1e-4)/qbar;
            asset_min[t] = asset_min[t+1]*qbar - income_min + 1e-4;
            if isnan( asset_min[t])
                asset_min[t]= asset_min[t+1]
            end 
        end
        #if t==N_ages 
        #    asset_min[t]=0
        #end
    end

    # Make k_max the largest amount the consumer could ever save, even
    # if they got the highest income and saved half of income every quarter
    for t = 1:N_ages
        income_max     = maximum(income_grid[t,:,:]);
        if t==1
            asset_max[t] = 0.95*income_max
        end
        asset_max[t+1] = min(asset_max[t]/Q[t,1,N_k] + 0.95*income_max , 100)
        if t==N_ages 
            asset_max[t+1] = asset_max[t];
        end 
    end

    # temporarily have everyone with the same grid
    max_asset_max = maximum(asset_max);
    min_asset_min = minimum(asset_min);

    for t = 1:N_ages
        N_p = ceil(Int64, N_k*0.70);
        N_n = N_k - N_p;

        #assets_p  = collect(LinRange(0.0, asset_max[t]^(1/α), N_p)).^α; 
        #assets_n  = collect(LinRange(asset_min[t], -1.0/N_n, N_n));   
        assets_p  = collect(LinRange(0.0, asset_max[t]^(1/α), N_p)).^α;  # optional convex grid for positive assets
        assets_n  = collect(LinRange(asset_min[t], 0.0, N_n+1));      # uniformly spaced points for negative assets

        # Asset grid for the given time period (linear below 0, nonlinear above 0)
        asset_grid[t, 1:N_n]     = assets_n[1:N_n];
        asset_grid[t, N_n+1:N_k] = assets_p;

        #double check that it's monotone:
        sort!(asset_grid[t,:]);
    end
    


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
    asset_grid::Array{Float64,1}, a::Float64, t::Int64, EVp_spn:: Schumaker, #Array{Float64,1}, #
    cbar::Float64, ψ::Float64, β::Float64, ϕ_1::Float64, ϕ_2::Float64, γ::Int64)::Float64
    # Note: b, s, y, ψ are levels. 
    
    # Get level of capital implied by next-period consumption choice
    ap = b + (1.0 - s)*a;
    if ap>asset_grid[length(asset_grid)] || ap<asset_grid[1]
        println("ap out of bounds: (ap,b,s,a)= ($ap, $b, $s, $a)");
    end
    
    EVp= SchumakerSpline.evaluate(EVp_spn,ap) ; #LinearInterpolation(asset_grid, EV_spn)(ap); #
    
    q_ap = max( LinearInterpolation(asset_grid, q_spn)(ap) , 0.0);  
    
    c = y + a*s - q_ap*b;
    util_penalty = 0.0;
    if c<1e-4
        util_penalty = c^2;
        c=1e-4;
        
        #println("negative c with b as $b and $instutil");
    end
    return u(c,cbar, γ) - util_penalty - ϕ(s, ϕ_1, ϕ_2) + β*ψ*EVp
end

"""
Solve for optimal b, given s 
"""

function max_bobjective(s::Float64, y::Float64, q_spn::Array{Float64,1}, #Schumaker, 
    atp::Array{Float64,1}, a::Float64, t::Int64, EVp_spn::Schumaker,
    cbar::Float64,r::Float64,ψ::Float64, β::Float64, ϕ_1_1::Float64, ϕ_1_2::Float64, γ::Int64)::Tuple{Float64,Float64}

    bmin = atp[1]   - (1.0 - s) * a ;
    bmax = atp[N_k] - (1.0 - s) * a ;
    
    #this portion tried to ensure c>0 when saving a lot:
    # 0 < y + a*s - q_(b + (1-s)*a)*bmax2; sometimes: (y+a*s)*(1+r)
    bmax_c0 = (y+a*s)*(1+r);
    if bmax_c0 > 0 # can save enough to get into positive assets, so q = 1/(1+r)
        bmax = min(bmax,bmax_c0);
    else 
        qmin = bmax_c0 + (1.0 - s)*a > atp[1] ?  LinearInterpolation(atp, q_spn)( bmax_c0 + (1.0 - s)*a ) : 
        LinearInterpolation(atp, q_spn)( bmax_c0 + (1.0 - s)*a )(atp[1] );
        bmax = min( (y+a*s)/qmin , bmax )
    end
       
    if bmax > bmin 
        res::Optim.OptimizationResults = optimize(b -> -1.0*bobjective(b, s, y, q_spn, atp, a,t, EVp_spn, cbar,ψ, β, ϕ_1_1, ϕ_1_2, γ), bmin, bmax);
        bopt = Optim.minimizer(res);
        return -1.0*Optim.minimum(res), bopt;  # the maximum
    else 
        println("age is $t and atp is $atp ");
        println("bmax < bmin ($bmax > $bmin) in max_bojb  at a = $a , s = $s")
        bopt = bmin;
        exit();
        return bobjective(bmin, s, y, q_spn, atp, a, t, EVp_spn, cbar, ψ, β, ϕ_1_1, ϕ_1_2, γ), bopt
    end 
end


function solve_ai( mod::model, y::Float64,atp::Array{Float64,1},a0::Array{Float64,1},
    EVp_spn::Schumaker, N_k::Int64, zi::Int64, εi::Int64,ai::Int64,t::Int64,ψ::Float64, qhr::Array{Float64,1})::Tuple{Float64,Float64,Float64,Float64,Float64}
    
    a0hr = a0[ai];

    apmin = atp[1];
    apmax = atp[N_k];
    # Brent's method for solving both s and b
    sopt::Float64 =  1.0;
    vopt::Float64 =  1.0;
    bopt::Float64 =  1.0; #just to set the scope
    
    if  a0hr >= 0 
        sopt          = 1.0;
        
        bmax::Float64 = min(apmax - (1.0 - sopt) * a0hr , (y+a0hr*sopt )*(1.0+ mod.r) );
        bmin::Float64 = apmin   - (1.0 - sopt) * a0hr  ;  
        if bmin > bmax
            #println("bmin > bmax: $bmin > $bmax at $t, $ai, and $a0hr")
            bopt = bmin;
            vopt = bobjective(bopt, sopt, y, qhr, atp, a0hr, t, EVp_spn, mod.cbar,ψ, mod.β, mod.ϕ_1, mod.ϕ_2, γ);
             println("bmin > bmax: $bmin > $bmax with $vopt, $bopt")
            
        else 
            b_solution_struct::Optim.OptimizationResults = optimize(b -> -1*bobjective(b, sopt, y, qhr, atp, a0hr, t,EVp_spn, mod.cbar,ψ, mod.β, mod.ϕ_1, mod.ϕ_2, γ), bmin, bmax);
            vopt          = -1.0 * Optim.minimum(b_solution_struct);
            bopt          = Optim.minimizer(b_solution_struct);
        end 
    else
        smax = 1.0;

        if y+ a0hr - qhr[1] *( atp[1]-a0hr ) <= 0 #consumption in this case
            #have to default a little bit:
            smax = min(1.0,  (1e-6 - y + qhr[1]*atp[1]-qhr[1]*a0hr)/(a0hr - a0hr*qhr[1])  )
        end
        if smax > 0
            res_s::Optim.OptimizationResults = optimize( shr-> -1.0*(max_bobjective(shr, y, qhr, atp, a0hr, t, EVp_spn, mod.cbar,mod.r,ψ, mod.β, mod.ϕ_1 , mod.ϕ_2, γ)[1]), 0.0,smax, 
                Brent(); iterations = 100,abs_tol=1e-5  );
            minout    = Optim.minimizer(res_s);
            sopt = minout[1];
        else 
            sopt = 0.0;
        end
        
        vopt,bopt = max_bobjective(sopt, y, qhr, atp, a0hr, t, EVp_spn, mod.cbar,mod.r,ψ,mod.β, mod.ϕ_1 , mod.ϕ_2, γ);
        
        
    end
    
    # Budget constraint
    # a_{t + 1} = 1 / q() * (a * s + y - c) + (1 - s) * a =>
    # c= y + a * s - q(a(),a' * b
    Apopt = bopt + (1.0 - sopt)*a0[ai];
    if Apopt > atp[N_k] || Apopt < atp[1]
        at1 = atp[1];
        println("Outside of A grid at $t, $ai, $εi, $zi : $Apopt < $at1 and $bopt, $sopt");
        #return V, B, C, S, Ap;
    else
        copt = y + a0[ai]*sopt - LinearInterpolation(atp, qhr)(Apopt)*bopt ;
    end
    
    return vopt,sopt,bopt,copt, Apopt
end


"""
Taking Q as given, solve for V,B,C,S,K
"""
function backwardSolve!(ms::sol, mod::model, 
    N_ages::Int64, N_k::Int64, N_z::Int64, N_ε::Int64, Q::Array{Float64,3})



    for t =  N_ages:(-1):1
        
        # a' grid
        atp = mod.asset_grid[t + 1, :];
        a0  = mod.asset_grid[t, :];

        # println("age $t atp: $atp a0: $a0")
        
        # Exogenous survival prob for the period
        ψ = survivalprobs[t];

        triloop = N_t; #lastiter==0 ? 1 : N_t

        for tri=1:triloop
            # Rebate
            transf = mod.transfer_grid[tri];


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
                EVp_grid = zeros(N_k); 
                for aii = 1:N_k
                    for eii = 1:N_ε
                        for zii = 1:N_z 
                            if t<N_ages
                                EVp_grid[aii] = ms.V[t+1, 1, aii, eii, zii ]*mod.εprobs[t,eii]*mod.zprobs[t,zi,zii] + EVp_grid[aii];
                            else 
                                c = mod.Ygrid[t, eii, zii] + transf + mod.asset_grid[t];
                                uc = c > 0.0 ? u(c,mod.cbar, γ) : u(1e-6,mod.cbar,γ);
                                EVp_grid[aii] = uc *mod.εprobs[t,eii]*mod.zprobs[t,zi,zii] + EVp_grid[aii];
                            end
                        end 
                    end
                end
                
                EVp_spn = Schumaker(atp,EVp_grid); 


                @inbounds for εi = 1:N_ε
                    y = mod.Ygrid[t, εi, zi] + transf;

                    @inbounds Threads.@threads for ai = 1:N_k #loop over a0

                        # With stochastic aging, not a thing:
                        # Right before certain death, the agent chooses to consume all that they can
                        if t == N_ages
                            ms.S[t, tri, ai, εi, zi] = 1.0;
                            ms.B[t, tri, ai, εi, zi] = 0.0;
                            ms.C[t, tri, ai, εi, zi] = a0[ai] + y + transf;
                            ms.V[t, tri, ai, εi, zi] = u(ms.C[t, tri, ai, εi, zi], mod.cbar, γ);
                            ms.Ap[t,tri, ai, εi, zi] = 0.0;
                            continue
                        end
                        
                        (vopt,sopt,bopt,copt, Apopt) = solve_ai(mod,y,atp,a0,EVp_spn,N_k,zi,εi,ai,t,ψ,qhr);
                        
                        ms.V[t, tri, ai, εi, zi] = vopt;
                        ms.S[t, tri, ai, εi, zi] = sopt;
                        ms.B[t, tri, ai, εi, zi] = bopt;
                        ms.C[t, tri, ai, εi, zi] = copt;
                        ms.Ap[t, tri, ai, εi, zi] = Apopt;
                        
                        
                    end
                    
                    #ms.Ap[t, tri, :, εi, zi] .= msAp_ai;
                    #ms.C[t, tri, :, εi, zi] .= msC_ai;
                    #ms.S[t, tri, :, εi, zi] .= msS_ai;
                    #ms.V[t, tri, :, εi, zi] .= msV_ai;
                    #ms.B[t, tri, :, εi, zi] .= msB_ai;

                end
            end
        end #tri loop over transfer
    end


    #return V, B, C, S, Ap
    return 1
end

"""
Solve for implied Q, given S
"""
function equilibriumQ(Q0::Array{Float64,3}, ms::sol, asset_grid::Array{Float64,2}, r::Float64, λ::Float64, εprobs::Array{Float64,2}, zprobs::Array{Float64,3})
    #  Dimension is N_ages, z, ap, b
    qimplied = similar(Q0);
    #N_ages = size(Q0)[1];
    #N_z=size(Q0)[2];
    #N_k=size(Q0)[3];
    #N_ϵ = size(εprobs)[2];
    for t = N_ages:-1:1
        #if t >= T - 2
        #    qimplied[t, :, :] .= 1.0 / (1.0 + r + I_{b>0} )
        #    continue
        #end

        for zi = 1:N_z
            for ai = 1:N_k
                qimplied[t, zi, ai] = 0.0;
                for εi =1:N_ε
                    for zii=1:N_z
                        #=if age:
                        # what to do with old dead guys? Replaced by young people exactly the same as them
                        PrAge_tp1 = ageprobs[t]

                        app_tp1 = ms.Ap[tp1, 1, ai, εi, zii];
                        app_tp1 = app_tp1 < asset_grid[tp1,1]   ? asset_grid[tp1,1] : app_tp1
                        app_tp1 = app_tp1 > asset_grid[tp1,N_k] ? asset_grid[tp1,N_k] : app_tp1
                        
                        shr_tp1 = t < N_ages ? ms.S[tp1, 1, ai, εi, zii] : ms.S[t, 1, ai, εi, zii];
                        qp_tp1  = LinearInterpolation(asset_grid[tp1,:],Q0[tp1,zii,:])(app_tp1);
                        =#
                        tp1 = t < N_ages ? t+1 : t;
                        
                        # introduce trembles:
                        #for aii=1:(1+2*N_trmbl)
                        app = ms.Ap[tp1, 1, ai, εi, zii]; #ms.Aptrmbl[t, 1, ai, εi, zii,aii] 
                        app = app < asset_grid[tp1,1]   ? asset_grid[tp1,1] : app #force in bounds (this shouldn't bind)
                        app = app > asset_grid[tp1,N_k] ? asset_grid[tp1,N_k] : app
                        shr = ms.S[tp1, 1, ai, εi, zii];
                    
                        qp      = LinearInterpolation(asset_grid[tp1  ,:],Q0[tp1  ,zii,:])(app);
                        qimplied[t, zi, ai] =  ( shr + (1-shr)*qp
                            )*εprobs[t,εi]*zprobs[t,zi,zii] +
                            qimplied[t, zi, ai] ;
                        #end
                        

                    end
                end
                #double check that these are defined and w/in natural bounds:
                if isfinite(qimplied[t, zi, ai]) && qimplied[t, zi, ai] > 0.0
                    qimplied[t, zi, ai] =  qimplied[t, zi, ai];
                elseif sum(ms.S[t,1,ai,:,:]) >= N_ϵ*N_z-1e-2
                    println("infinite q, all 1 at $t,$i,$ai");
                    qimplied[t, zi, ai] = 1.0 / (1.0 + r);
                else
                    qhr = qimplied[t, zi, ai]
                    println("infinite/non-positive q ($qhr), not 1 at $t,$zi,$ai");
                    qimplied[t, zi, ai] = 0.0;
                end
                if asset_grid[t,ai] >= 0.0 
                    qimplied[t,zi,ai] =   qimplied[t,zi,ai]./(1.0 + r); # needs to be w/in the t loop otherwise for periods with T>=2 will double multiply by 1/(1+r)
                else 
                    qimplied[t,zi,ai] =   qimplied[t,zi,ai]./(1.0 + r + λ); 
                end
            end #looping on ai 
        end #looping on zi
        
    end
    return qimplied
end

function draw_shocks!(mod::model, ht::hists,seed::Int64=12281951)
    Random.seed!(seed)

    ht.a0qtlhist .= rand(Nsim);
    
    ageproto = rand(Nsim,Thist*freq);
    @inbounds for i=1:Nsim
        ht.agehist[i,1] = ceil(Int64,ageproto[i,1]*N_ages);
        @inbounds for t=2:(Thist*freq)
            agp1 = ageproto[i,t] < ageprobs[ht.agehist[i,t-1]] ?
                ht.agehist[i,t-1] + 1 : ht.agehist[i,t-1];
            ht.agehist[i,t] = agp1<=N_ages ? agp1 : 1
        end 
    end
    
    εproto =  rand(Nsim,Thist*freq);
    zproto =  rand(Nsim,Thist*freq);
    cumεprobs = zeros(N_ages,N_ε+1);
    cumzprobs = zeros(N_ages,N_z,N_z+1);
    cumzergod = zeros(N_ages,N_z+1);
    for t=1:N_ages 
        cumεprobs[t,2:(N_ε+1)] .= cumsum(mod.εprobs[t,:]);
        for zi=1:N_z
            cumzprobs[t,zi,2:(N_z+1)] .= cumsum(mod.zprobs[t,zi,:]);
        end 
        if t>1
            zergod_t = mod.zprobs[t,:,:]^(1000*freq);
            cumzergod[t,2:(N_z+1)] .= cumsum(zergod_t[1,:]);
        else 
            cumzergod[1,2:(N_z+1)] .= cumsum(ones(N_z)./N_z );
        end
    end 
    
    @inbounds for i=1:Nsim
        ht.zidxhist[i,1] = sum( zproto[i,1].>cumzergod[ ht.agehist[i,1],: ] );
        @inbounds for t=1:(Thist*freq)
            if t>1
                izlast = ht.zidxhist[i,t-1];
                ht.zidxhist[i,t] = sum( zproto[i,t] .>= cumzprobs[ht.agehist[i,t], izlast,: ]  );
            end
            ht.εidxhist[i,t] = sum( εproto[i,t] .>= cumεprobs[ht.agehist[i,t], : ]  );
            income_hr = mod.Ygrid[ht.agehist[i,t],ht.εidxhist[i,t],ht.zidxhist[i,t]];
            ht.incomehist[i,t] = income_hr;
        end 
    end
end 

function sim_hists!(mod::model, ht::hists, ms::sol, mmt::moments,N_ages::Int64, N_k::Int64,N_z::Int64, N_ε::Int64)

    ithist = zeros(Nsim,Thist*freq);
    asset_alpha = ones(N_ages,N_ε, N_z); asset_beta=ones(N_ages,N_ε, N_z);  #will iteratively fit a beta distribution 
    a95         = ones(N_ages,N_ε, N_z);
    for dist_iter = 1:1
        
        @inbounds for i=1:Nsim 

            @inbounds for it=1:(Thist*freq)
                ei_hr = ht.εidxhist[i,it];
                zi_hr = ht.zidxhist[i,it];
                age_hr=  ht.agehist[i,it];
                ithist[i,it] = it;
                if it==1 # Initialize
                    a95_hr = a95[it,ei_hr,zi_hr];
                    if age_hr >1 && asset_alpha[age_hr,ei_hr,zi_hr] > 0 && asset_beta[age_hr,ei_hr,zi_hr] > 0
                        ht.ahist[i,1] = quantile(Beta(asset_alpha[age_hr,ei_hr,zi_hr],asset_beta[age_hr,ei_hr,zi_hr]),ht.a0qtlhist[i] )*(a95_hr - mod.asset_grid[age_hr,1]) + mod.asset_grid[age_hr,1];
                    elseif age_hr==1 #young all go to the same place
                        ht.ahist[i,1] = 0.0;
                    elseif asset_alpha[age_hr,ei_hr,zi_hr] <=0 #if there was a degenerate distribution I stored the average as -asset_alpha
                        ht.ahist[i,1] = -asset_alpha[age_hr,ei_hr,zi_hr]*(a95_hr - mod.asset_grid[age_hr,1]) + mod.asset_grid[age_hr,1];
                    end 
                end  # end if it==1   
                aL = 1; aLwt = 0.0;
                if age_hr > 1
                    aL = sum( ht.ahist[i,it] .>=  mod.asset_grid[age_hr,:] );
                    aL = aL>= N_k ? N_k-1 : (aL <1 ? 1 : aL);
                    aLwt =  (mod.asset_grid[age_hr,aL+1] - ht.ahist[i,it])/(mod.asset_grid[age_hr,aL+1] - mod.asset_grid[age_hr,aL])
                else #for new-borns, start everyone with no assets or debt
                    ht.ahist[i,it] = 0.0;
                    aL = sum( ht.ahist[i,it] .>=  mod.asset_grid[age_hr,:] );
                    aL = aL>= N_k ? N_k-1 : (aL <1 ? 1 : aL);
                    aLwt =  (mod.asset_grid[age_hr,aL+1] - ht.ahist[i,it])/(mod.asset_grid[age_hr,aL+1] - mod.asset_grid[age_hr,aL])
                end 
                if it<=(Thist*freq - 1)
                    ap = aLwt*ms.Ap[age_hr,1,aL,ei_hr,zi_hr ] +
                        (1.0 - aLwt)*ms.Ap[ age_hr,1,aL+1,ei_hr,zi_hr ];
                    ht.ahist[i,it+1] = ap;
                    apL = sum( ap .>=  mod.asset_grid[age_hr,:] );
                    apL = apL>= N_k ? N_k-1 : (apL <1 ? 1 : apL);
                    apLwt =  (mod.asset_grid[age_hr,aL+1] - ap)/(mod.asset_grid[age_hr,aL+1] - mod.asset_grid[age_hr,aL])
                    ht.qhist[i,it] = apLwt*ms.Q0[age_hr,zi_hr,apL] + 
                        (1.0 - apLwt)*ms.Q0[age_hr,zi_hr,apL+1];
                end
                ht.chist[i,it] = aLwt*ms.C[age_hr,1,aL,ei_hr,zi_hr] + 
                    (1.0 - aLwt)*ms.C[age_hr,1,aL+1,ei_hr,zi_hr];
                
                ht.shist[i,it] = aLwt*ms.S[age_hr,1,aL,ei_hr,zi_hr] + 
                    (1.0 - aLwt)*ms.S[age_hr,1,aL+1,ei_hr,zi_hr];
                if ht.shist[i,it] <.98 && ht.ahist[i,it] >0
                    shr = ht.shist[i,it];
                    ahr = ht.ahist[i,it];
                    #println("With age $age_hr, wrong s,a: $shr,$ahr with aL,aLwt: $aL,$aLwt")
                    ht.shist[i,it] =1.0;
                end
                ht.bhist[i,it] = aLwt*ms.B[age_hr,1,aL,ei_hr,zi_hr] + 
                    (1.0 - aLwt)*ms.B[age_hr,1,aL+1,ei_hr,zi_hr];

                #compute mpc, mpd,mpb:
                ht.mpbhist[1,i,it] = (aLwt*ms.B[age_hr,2,aL,ei_hr,zi_hr] + 
                    (1.0 - aLwt)*ms.B[age_hr,2,aL+1,ei_hr,zi_hr] -
                    ht.bhist[i,it])/mod.transfer_grid[2]
                ht.mpbhist[2,i,it] = (aLwt*ms.B[age_hr,3,aL,ei_hr,zi_hr] + 
                    (1.0 - aLwt)*ms.B[age_hr,3,aL+1,ei_hr,zi_hr] -
                    ht.bhist[i,it])/mod.transfer_grid[3]
                                
                ht.mpchist[1,i,it] = (aLwt*ms.C[age_hr,2,aL,ei_hr,zi_hr] + 
                    (1.0 - aLwt)*ms.C[age_hr,2,aL+1,ei_hr,zi_hr] -
                    ht.chist[i,it])/mod.transfer_grid[2]
                ht.mpchist[2,i,it] = (aLwt*ms.C[age_hr,3,aL,ei_hr,zi_hr] + 
                    (1.0 - aLwt)*ms.C[age_hr,3,aL+1,ei_hr,zi_hr] -
                    ht.chist[i,it])/mod.transfer_grid[3]
                    
                ht.mpshist[1,i,it] = (aLwt*ms.S[age_hr,2,aL,ei_hr,zi_hr] + 
                    (1.0 - aLwt)*ms.S[age_hr,2,aL+1,ei_hr,zi_hr] -
                    ht.shist[i,it])/mod.transfer_grid[2]
                ht.mpshist[2,i,it] = (aLwt*ms.S[age_hr,3,aL,ei_hr,zi_hr] + 
                    (1.0 - aLwt)*ms.S[age_hr,3,aL+1,ei_hr,zi_hr] -
                    ht.shist[i,it])/mod.transfer_grid[3]

                if it<=(Thist*freq - 1)
                    #a' = b + (1.0 - s) a 
                    ap_tr2 = aLwt*ms.B[age_hr,2,aL,ei_hr,zi_hr] + 
                        (1.0 - aLwt)*ms.B[age_hr,2,aL+1,ei_hr,zi_hr] +
                        (1.0 - aLwt*ms.S[age_hr,2,aL,ei_hr,zi_hr] - 
                        (1.0 - aLwt)*ms.S[age_hr,2,aL+1,ei_hr,zi_hr])*ht.ahist[i,it];
                    apL = sum( ap_tr2 .>=  mod.asset_grid[age_hr,:] );
                    apL = apL>= N_k ? N_k-1 : (apL <1 ? 1 : apL);
                    apLwt =  (mod.asset_grid[age_hr,aL+1] - ap_tr2)/(mod.asset_grid[age_hr,aL+1] - mod.asset_grid[age_hr,aL])
                    q_tr2 = apLwt*ms.Q0[age_hr,zi_hr,apL] + 
                        (1.0 - apLwt)*ms.Q0[age_hr,zi_hr,apL+1];
                    ht.mpqhist[1,i,it] = (q_tr2 - ht.qhist[i,it])/mod.transfer_grid[2];

                    ap_tr3 = aLwt*ms.B[age_hr,3,aL,ei_hr,zi_hr] + 
                        (1.0 - aLwt)*ms.B[age_hr,3,aL+1,ei_hr,zi_hr] +
                        (1.0 - aLwt*ms.S[age_hr,3,aL,ei_hr,zi_hr] - 
                        (1.0 - aLwt)*ms.S[age_hr,3,aL+1,ei_hr,zi_hr])*ht.ahist[i,it];
                    apL = sum( ap_tr3 .>=  mod.asset_grid[age_hr,:] );
                    apL = apL>= N_k ? N_k-1 : (apL <1 ? 1 : apL);
                    apLwt =  (mod.asset_grid[age_hr,aL+1] - ap_tr3)/(mod.asset_grid[age_hr,aL+1] - mod.asset_grid[age_hr,aL])
                    q_tr3 = apLwt*ms.Q0[age_hr,zi_hr,apL] + 
                        (1.0 - apLwt)*ms.Q0[age_hr,zi_hr,apL+1];
                    ht.mpqhist[2,i,it] = (q_tr3 - ht.qhist[i,it])/mod.transfer_grid[3]
                end 
            end #end it=1:Tsim+burnin
        end #end i =1:Nsim

        # update asset_alpha, asset_beta
        for t=1:N_ages
            for ei=1:N_ε
                for zi=1:N_z
                    mask_hr = (ht.εidxhist.==ei) .& (ht.zidxhist.==zi) .& (ht.agehist .== t) .& (ithist .>= (Tburnin*freq));
                    if(sum(mask_hr)>0)
                        
                        a95[t,ei,zi] = quantile( ht.ahist[mask_hr],.95 );
                        ahist_hr = (ht.ahist[mask_hr] .- mod.asset_grid[t,1])./(a95[t,ei,zi] - mod.asset_grid[t,1]);
                        ahist_hr[ahist_hr .> 1.] .= 1;
                        Ea_hr = 0.5; Va_hr = 0.0;
                        #masktot = sum(mask_hr)
                        #println("masktot $masktot")
                        Ea_hr = mean(ahist_hr);
                        Va_hr = var(ahist_hr);
                    else 
                        Va_hr = 0.
                    end
                    #if Ea_hr<=0 || Va_hr<0 || isnan(Ea_hr) || isnan(Va_hr)
                    #    println(" Ea_hr,Va_hr, t = $Ea_hr, $Va_hr, $t, $ei, $zi ");
                    #end 
                    if Va_hr>0 && !isnan(Va_hr)
                        asset_alpha[t,ei,zi] =      Ea_hr *(Ea_hr*(1.0-Ea_hr)/Va_hr-1.0);
                        asset_beta[t,ei,zi]  = (1.0-Ea_hr)*(Ea_hr*(1.0-Ea_hr)/Va_hr-1.0);
                    else #in the case where everyone goes to the same place:
                        Ea_hr = (mean(ht.ahist[mask_hr])- mod.asset_grid[t,1])/(mod.asset_grid[t,N_k] - mod.asset_grid[t,1]);
                        if(isnan(Ea_hr))
                            asset_alpha[t,ei,zi] = 1.0;
                        else 
                            asset_alpha[t,ei,zi] = Ea_hr/(1.0-Ea_hr);
                        end
                        asset_beta[t,ei,zi]  = 1.0;
                    end                
                end
            end
        end
        asset_alpha_t1 = mean(asset_alpha[:,2,2]);
        asset_beta_t1 = mean(asset_beta[:,2,2]);
        println("Asset α = $asset_alpha_t1, β=$asset_beta_t1");
    end #iterated simulations

    Ea_tot = mean(ht.ahist);
    Es_tot = mean(ht.shist);

    mmt.fr_neg_asset = 0.0
    mmt.avg_indiv_debt_income=0.0
    mmt.avg_s=0.0
    mmt.fr_s1=0.0
    denom_assets = 0.0;
    mmt.extent_delinq = 0.0;
    mmt.avg_debt = 0.0;
    mmt.avg_income_neg_asset = 0.0;
    mmt.mpc_ac0 = 0.0;

    Tsim0 = (Tburnin+1)*freq;
    TsimT = (Tsim+Tburnin)*freq;

    MPCs_negassets = ones(Nsim*(TsimT-Tsim0+1),5)
    MPCs_posassets = ones(Nsim*(TsimT-Tsim0+1),3)
    idx_neg::Int64 = 0
    idx_pos::Int64 = 0
    @inbounds for i=1:Nsim 
        @inbounds for it=Tsim0:TsimT
            if ht.ahist[i,it] < 0.0 
            #    println("negative assets! $i, $it")
                mmt.fr_neg_asset =  1.0 + mmt.fr_neg_asset ;
                mmt.avg_income_neg_asset = ht.incomehist[i,it] + mmt.avg_income_neg_asset 
                mmt.avg_debt = -ht.ahist[i,it] + mmt.avg_debt;
                mmt.avg_indiv_debt_income =  -ht.ahist[i,it] /ht.incomehist[i,it] + mmt.avg_indiv_debt_income;
                mmt.avg_s = ht.shist[i,it] + mmt.avg_s
                mmt.extent_delinq = (-ht.ahist[i,it])*(1.0-ht.shist[i,it]) + mmt.extent_delinq
                mmt.fr_s1 = ht.shist[i,it] > 0.99 ? 1.0 + mmt.fr_s1 : mmt.fr_s1;

                idx_neg += 1;
                
                MPCs_negassets[idx_neg,1] = -ht.ahist[i,it] ;
                MPCs_negassets[idx_neg,2] =  ht.mpchist[i,it];
                MPCs_negassets[idx_neg,3] =  ht.mpbhist[i,it];
                MPCs_negassets[idx_neg,4] =  ht.mpshist[i,it];
                MPCs_negassets[idx_neg,5] =  ht.mpqhist[i,it];

            else 
                denom_assets = 1.0 + denom_assets;
                mmt.mpc_ac0 = ht.mpchist[i,it] + mmt.mpc_ac0;
                MPCs_posassets[idx_pos,1] = ht.ahist[i,it] ;
                MPCs_posassets[idx_pos,2] = ht.mpchist[i,it] ;
                MPCs_posassets[idx_pos,3] = ht.mpbhist[i,it] ;

                idx_pos += 1;

            end
        end
    end
    mmt.avg_income = mean(ht.incomehist);
    mmt.avg_income_neg_asset = mmt.avg_income_neg_asset / mmt.fr_neg_asset;
    mmt.avg_indiv_debt_income = mmt.avg_indiv_debt_income/mmt.fr_neg_asset;
    mmt.avg_s = mmt.avg_s/mmt.fr_neg_asset;
    mmt.fr_s1 = mmt.fr_s1/mmt.fr_neg_asset;
    mmt.avg_debt = mmt.avg_debt/mmt.fr_neg_asset;
    mmt.extent_delinq = mmt.extent_delinq/mmt.fr_neg_asset;
    
    mmt.mpc_ac0 = mmt.mpc_ac0 / denom_assets
    mmt.mpc_cor_neg_asset = cor(MPCs_negassets[1:idx_neg,:]);
    mmt.mpc_cor_neg_asset = cor(MPCs_posassets[1:idx_pos,:]);
    
    mmt.fr_neg_asset = mmt.fr_neg_asset/(mmt.fr_neg_asset+ denom_assets);

    println("The average asset position is $Ea_tot and pay-back is $Es_tot");

end


function steadystateDist!( ms::sol, mod::model, survivalprobs::Vector{Float64}, N_ages, N_k, N_z, N_ε)

    survivalprobs_cml = sum(1.0 .- survivalprobs[1:(N_ages-1)]) #fraction dead before period N_ages+1
    survivalprobs[N_ages] = 0.0;
    N_krnd = Int64(N_k/5);
    asteps = Int64(N_k/N_krnd);
    asset_grid_rnd = zeros(N_ages,N_krnd);
    asset_space    = zeros(N_ages,N_krnd);
    for t=1:(N_ages-1)
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
    Nstate = N_ages*N_krnd*N_ε*N_z;
    a0idx = findlast(  asset_grid_rnd[2,:] .<= 0.0 );  #mod.asset_grid[2,:] .<= 0.0)
    # state is (t-1)*N_k*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi
    flatstatemat = spzeros(Float64,Nstate,Nstate);
    for t=1:N_ages
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
                    if t<(N_ages-1)
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

    Npara = N_ages*N_krnd*N_ε*N_z;
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
    for t=1:N_ages
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
