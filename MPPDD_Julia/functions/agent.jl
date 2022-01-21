"""
Return the asset grids, given the total (labor) income grid, N_ages, and N_a
α = 1 for uniformly spaced gripoints.
"""
function assetGrid(N_ages::Int64, income_grid::Array{Float64,3}, N_a::Int64,Q::Array{Float64}, asset_grid_old::Array{Float64,2}, r::Float64, α::Float64 = 2.0)

    asset_max  = zeros(N_ages+1);
    asset_min  = zeros(N_ages+1);
    asset_grid = zeros(N_ages+1, N_a); # N_ages + 1 = year after death, define for convenience in backwardSolve

    # Minimum assets (Maximum debt)
    if natborrowlimit == true
		for t = N_ages:-1:1
		    qhr = Q[t,1,:];
		    income_min = minimum(income_grid[t,:,:]);
		    if t<N_ages
		        qbar = max(LinearInterpolation(asset_grid_old[t+1,:], qhr; extrapolation_bc=Line())(asset_min[t+1] ) , 1e-2);
                qbar = min(qbar, 1.0);
		        if isnan(qbar)
		            qbar = 1/(1.0+ r );
		        end
		    else
		        qbar = 1.0/(1.0 +r );
		    end
		    if t>=age_retire
		        #spps pay back the whole thing:
		        asset_min[t] = qbar*asset_min[t+1] - income_min + 1e-4;
		    else
		        #spps not pay back at all
		        #asset_min[t] = asset_min[t+1] - (income_min + 1e-4)/qbar;
		        asset_min[t] = qbar*asset_min[t+1] - income_min + 1e-4;
		        if isnan( asset_min[t])
		            asset_min[t]= asset_min[t+1]
		        end
		    end
            #double-check that it's never increasing
            asset_min[t] = asset_min[t] > asset_min[t+1] ? asset_min[t+1] : asset_min[t] ;
		end
    else
		for t = N_ages:-1:1
		    #where does QB stop increasing?
		    if t<N_ages
		        if minimum(asset_grid_old[t+1,:])<0
		            asset_min_tz = Inf;
		            for zi=1:N_z
		                qa_aip1 = 0.0;
		                for ai = N_a:-1:1

		                    if asset_grid_old[t+1,ai]<0
		                        qhr = Q[t,zi,:];
		                        qbar = LinearInterpolation(asset_grid_old[t+1,:], qhr; extrapolation_bc=Line())(asset_grid_old[t+1,ai]);
		                        if isnan(qbar)
		                            qbar = qhr[1];
		                        end
		                        qa = -qbar *asset_grid_old[t+1,ai]
		                        if qa < qa_aip1
		                            asset_min_tz = asset_grid_old[t+1,ai] < asset_min_tz ? asset_grid_old[t+1,ai] : asset_min_tz;
		                            #println(" The ai is $ai ");
		                            break;
		                        else
		                            qa_aip1 = qa;
		                        end
		                    end
		                end
		            end
		            asset_min[t+1] = asset_min_tz <Inf ? asset_min_tz : asset_grid_old[t+1,1];
		            if asset_min[t+1] == asset_grid_old[t+1,1] # we went too high, and never get to the hump in q
		                #println("Went to high with min asset at $t")
		                qbar = LinearInterpolation(asset_grid_old[t+1,:], Q[t,1,:]; extrapolation_bc=Line())(asset_grid_old[t+1,1]);
		                income_min = minimum(income_grid[t,:,:]);
		                asset_min[t] = asset_min[t+1]*qbar - income_min + 1e-4;
		                if isnan( asset_min[t])
		                    asset_min[t]= asset_min[t+1]
		                end
		            end
		        else #never goes negative
		            qbar = LinearInterpolation(asset_grid_old[t+1,:], Q[t,zi,:]; extrapolation_bc=Line())(asset_grid_old[t+1,1]);
		            income_min = minimum(income_grid[t,:,:]);
		            asset_min[t] = asset_min[t+1]*qbar - income_min + 1e-4;
		        end
		    end
		end
	end #end switch for natural or endog borrowing limit
    # Make k_max the largest amount the consumer could ever save, even
    # if they got the highest income and saved half of income every quarter
    for t = 1:N_ages
        income_max     = maximum(income_grid[t,:,:]);
        income_mean    = mean(income_grid[t,:,N_z]);
        if t==1
            asset_max[t] = income_max
        end
        asset_max[t+1] =  income_max*5; #min(asset_max[t]/Q[t,1,N_a] + income_max , income_max*5)
        asset_max[t+1] = max(asset_max[t+1] , asset_max[t])
        if t==N_ages
            asset_max[t+1] = asset_max[t];
        end
    end


	if unevengrid == false
		for t = 1:N_ages
		    #asset_grid[t, :] = LinRange(asset_min[t],asset_max[t],N_a);
            N_p = floor(Int64, N_a*0.34);
		    N_n = N_a - N_p;
            asset_grid[t, 1:N_n] = LinRange(asset_min[t],0.0,N_n);
            asset_grid[t, N_n:N_a] = LinRange(0.0,asset_max[t],N_p + 1);
		end
	else
		for t = 1:N_ages
		    N_p = floor(Int64, N_a*0.50);
		    N_n = N_a - N_p;

		    #assets_p  = collect(LinRange(0.0, asset_max[t]^(1/α), N_p)).^α;
		    #assets_n  = collect(LinRange(asset_min[t], -1.0/N_n, N_n));
		    assets_p  = collect(LinRange(0.0, asset_max[t]^(1/α), N_p+1)).^α;  # optional convex grid for positive assets
		    assets_n  = collect(LinRange(0.0,(-asset_min[t])^(1/α), N_n)).^α;      # optional convex grid for negative assets

		    # Asset grid for the given time period (linear below 0, nonlinear above 0)
		    asset_grid[t, 1:N_n]   = -assets_n[N_n:-1:1];
		    asset_grid[t, N_n:N_a] = assets_p;

		    #double check that it's monotone:
		    sort!(asset_grid[t,:]);
		end
    end

    return asset_grid
end

"""
Define the ϕ function
"""
@inline function ϕ(s::Float64, ϕ_1 ::Float64,ϕ_2::Float64,ϕ_3::Float64)
    return s<1 ? ϕ_1 *(1 -  exp(-ϕ_2*(s-1.0)))/(1.0-exp(ϕ_2)) - ϕ_3 : 0.0;
end

"""
CRRA utility
"""
@inline function u(c::Float64,s::Float64, γ::Int64,ϕ_1 ::Float64,ϕ_2::Float64,ϕ_3::Float64)
    # Nudge to avoid log(0) - type errors during runtime
    if γ == 1
        #return log(max(1e-8, c - cbar))
        return log( max(1e-8,c)) -  ϕ(s, ϕ_1,ϕ_2,ϕ_3)
        #return log( max(1e-8,c)* (ϕ_1-  ϕ(s, ϕ_1,ϕ_2,ϕ_3)) )
    else
        return ((c )^(1 - γ) - 1)/(1 - γ) -  ϕ(s, ϕ_1,ϕ_2,ϕ_3)
        #return ((c *(ϕ_1-  ϕ(s, ϕ_1,ϕ_2,ϕ_3)))^(1 - γ) - 1)/(1 - γ)
    end
end

"""
Define current utility + continuation value
"""
function bobjective(b::Float64, s::Float64, y::Float64, q_spn::Schumaker, #Array{Float64,1}, #
    atp_grid::Array{Float64,1}, a::Float64, t::Int64, EVp_spn::Array{Float64,1}, #Schumaker, #
    reclaimrate::Float64, ψ::Float64, β::Float64, ϕ_1::Float64, ϕ_2::Float64, ϕ_3::Float64, γ::Int64)::Float64
    # Note: b, s, y, ψ are levels.

    # Get level of capital implied by next-period consumption choice
    ap = b + (1.0 - s)*a*reclaimrate;

    EVp= LinearInterpolation(atp_grid, EVp_spn,extrapolation_bc=Line())(ap); #SchumakerSpline.evaluate(EVp_spn,ap) ; #
    q_ap = max( SchumakerSpline.evaluate(q_spn,ap) , 0.0); #max( CubicSplineInterpolation(atp_grid, q_spn)(ap) , 0.0);

    if ap>atp_grid[length(atp_grid)] || ap<atp_grid[1]
        if print_lev>0
            println("ap out of bounds: (ap,b,s,a)= ($ap, $b, $s, $a)");
        end 
        ap0 = ap<atp_grid[1] ? atp_grid[1] : atp_grid[length(atp_grid)]
        EVp=  LinearInterpolation(atp_grid, EVp_spn,extrapolation_bc=Line())(ap); # SchumakerSpline.evaluate(EVp_spn,ap0) - (ap-ap0)^2;
        q_ap = max( SchumakerSpline.evaluate(q_spn,ap0) , 0.0);

    end

    c = y + a*s - q_ap*b;
    util_penalty = 0.0;
    if c<1e-7
        util_penalty = (c-1e-7)^2;
        if print_lev >0
            println("negative c with b as $b and c as $c");
        end
        c=1e-7;
    end
    #return u(c,cbar, γ, ϕ_1, ϕ_2, ϕ_3) - util_penalty - ϕ(s, ϕ_1, ϕ_2, ϕ_3) + β*ψ*EVp
    return u(c,s, γ, ϕ_1, ϕ_2, ϕ_3) - util_penalty  + β*ψ*EVp
end

"""
Solve for optimal b, given s
"""

function max_bobjective(s::Float64, y::Float64, q_spn::Schumaker, #Array{Float64,1}, #
    atp::Array{Float64,1}, a::Float64, t::Int64, EVp_spn::Array{Float64,1}, #Schumaker,
    reclaimrate::Float64,r::Float64,ψ::Float64, β::Float64, ϕ_1::Float64, ϕ_2::Float64, ϕ_3::Float64, γ::Int64)::Tuple{Float64,Float64}

    bmin = atp[1]   - (1.0 - s) * a ;
    bmax = atp[N_a] - (1.0 - s) * a ;

    #this portion tried to ensure c>0 when saving a lot:
    # 0 < y + a*s - q_(b + (1-s)*a)*bmax2; sometimes: (y+a*s)*(1+r)
    bmax_c0 = (y+a*s)*(1+r);
    # can save enough to get into positive assets, so q = 1/(1+r)
    # check if I'm consuming positive and
    if (y - bmax_c0/(1.0 + r) + a*s) > 0 && bmax_c0 + (1-s)*a >0
        bmax = min(bmax,bmax_c0);

    else #can't save enough to get positive assets
        #qmin = bmax_c0 + (1.0 - s)*a > atp[1] ?  CubicSplineInterpolation(atp, q_spn)( bmax_c0 + (1.0 - s)*a ) :
        #CubicSplineInterpolation(atp, q_spn)( bmax_c0 + (1.0 - s)*a )(atp[1] );
        qmin = bmax_c0 + (1.0 - s)*a > atp[1] ? SchumakerSpline.evaluate(q_spn, bmax_c0 + (1.0 - s)*a ) :
            SchumakerSpline.evaluate(q_spn, atp[1] );
        bmax = qmin > 0 ? min( (y+a*s)/qmin , bmax ) : bmin;
    end

    if bmax > bmin
        res::Optim.OptimizationResults = optimize(b -> -1.0*bobjective(b, s, y, q_spn, atp, a,t, EVp_spn, reclaimrate,ψ, β, ϕ_1, ϕ_2, ϕ_3, γ), bmin, bmax);
        bopt = Optim.minimizer(res);
        return -1.0*Optim.minimum(res), bopt;  # the maximum
    else
        # println("age is $t and atp is $atp ");
        # println("bmax < bmin ($bmax < $bmin) in max_bojb  at a = $a , s = $s")
        bopt = bmin;

        return bobjective(bmin, s, y, q_spn, atp, a, t, EVp_spn, reclaimrate, ψ, β, ϕ_1, ϕ_2,ϕ_3, γ), bopt
    end
end


function solve_ai( mod::model, y::Float64,atp::Array{Float64,1},a0::Array{Float64,1},EVp_spn::Array{Float64,1}, #Schumaker,
    N_a::Int64, zi::Int64, εi::Int64,ai::Int64,t::Int64,ψ::Float64,q_spn::Schumaker #Array{Float64,1}
    )::Tuple{Float64,Float64,Float64,Float64,Float64}  # qhr::Array{Float64,1})::Tuple{Float64,Float64,Float64,Float64,Float64}

    a0hr = a0[ai];

    apmin = atp[1];
    apmax = atp[N_a];
    # Brent's method for solving both s and b
    sopt::Float64 =  1.0;
    vopt::Float64 =  1.0;
    bopt::Float64 =  1.0; #just to set the scope


    #q_spn = Schumaker(atp,qhr);

    if  a0hr >= 0 || fullcommit == true
        sopt          = 1.0;

        bmax::Float64 = min(apmax  , (y+a0hr*sopt )*(1.0+ mod.r) );
        bmin::Float64 = apmin    ;
        if bmin > bmax
            #println("bmin > bmax: $bmin > $bmax at $t, $ai, and $a0hr")
            bopt = bmin;
            vopt = bobjective(bopt, sopt, y, q_spn, atp, a0hr, t, EVp_spn, mod.reclaimrate,ψ, mod.β, mod.ϕ_1, mod.ϕ_2, mod.ϕ_3, γ);
            if print_lev>0
                println("bmin > bmax: $bmin > $bmax with $vopt, $bopt")
            end

        else
            b_solution_struct::Optim.OptimizationResults = optimize(b -> -1*bobjective(b, sopt, y, q_spn, atp, a0hr, t,EVp_spn, mod.reclaimrate,ψ, mod.β, mod.ϕ_1, mod.ϕ_2, mod.ϕ_3, γ), bmin, bmax);
            vopt          = -1.0 * Optim.minimum(b_solution_struct);
            bopt          = Optim.minimizer(b_solution_struct);
        end
    else
        smax = 1.0;
		#if constQ == true
		#    qhr1 = q_spn[1];
		#else
		    qhr1 = SchumakerSpline.evaluate(q_spn,atp[1]);
		#end
        if y+ a0hr - qhr1 *( atp[1]-a0hr ) <= 0 #consumption in this case
            #have to default a little bit:
            smax = min(1.0,  (1e-6 - y + qhr1*atp[1]-qhr1*a0hr)/(a0hr - a0hr*qhr1)  )
        end
        if smax > 0
            res_s::Optim.OptimizationResults = optimize( shr-> -1.0*(max_bobjective(shr, y, q_spn, atp, a0hr, t, EVp_spn, mod.reclaimrate,mod.r,ψ, mod.β, mod.ϕ_1 , mod.ϕ_2, mod.ϕ_3, γ)[1]), 0.0,smax,
                Brent(); iterations = 100,abs_tol=1e-5  );
            minout    = Optim.minimizer(res_s);
            sopt = minout[1];
        else
            sopt = 0.0;
        end

        vopt,bopt = max_bobjective(sopt, y, q_spn, atp, a0hr, t, EVp_spn, mod.reclaimrate,mod.r,ψ,mod.β, mod.ϕ_1 , mod.ϕ_2, mod.ϕ_3, γ);


    end

    # Budget constraint
    # a_{t + 1} = 1 / q() * (a * s + y - c) + (1 - s) * a =>
    # c= y + a * s - q(a(),a' * b
    Apopt = bopt + (1.0 - sopt)*a0[ai]*mod.reclaimrate;
    if Apopt > atp[N_a] || Apopt < atp[1]
        at1 = atp[1];
        if print_lev>0
            println("Outside of A grid at $t, $ai, $εi, $zi : $Apopt < $at1 and $bopt, $sopt");
        end 
        #return V, B, C, S, Ap;
        Apopt = Apopt > atp[N_a] ? atp[N_a] : Apopt > atp[1]
    end

    copt = y + a0[ai]*sopt - SchumakerSpline.evaluate(q_spn,Apopt)*bopt ;



    return vopt,sopt,bopt,copt, Apopt
end


"""
Taking Q as given, solve for V,B,C,S,K
"""
function backwardSolve!(ms::sol, mod::model,
    N_ages::Int64, N_a::Int64, N_z::Int64, N_ε::Int64, Q::Array{Float64,3})

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
            #@inbounds Threads.@threads for idx = 1:(N_z*N_ε*N_a)
            #for idx = 1:(N_z*N_ε*N_a)
            #    ai = floor(Int,mod(idx-1,N_z*N_ε*N_a)/N_ε/N_z )+1;
            #    zi = floor(Int,mod(idx-1,N_z*N_ε)/N_ε )+1;
            #    εi =           mod(idx-1,N_ε)+1;
            #    fidx = (ai-1)*N_ε*N_a + (zi-1)*N_ε + εi
                qhr    = Q[t, zi,:];

                # Expected value grid next period, given income states and assets
                EVp_grid = zeros(N_a);
                for aii = 1:N_a
                    for eii = 1:N_ε
                        for zii = 1:N_z
                            if t<N_ages
                                EVp_grid[aii] = ms.V[t+1, 1, aii, eii, zii ]*mod.εprobs[t,eii]*mod.zprobs[t,zi,zii] + EVp_grid[aii];
                            else
                                c = mod.Ygrid[t, eii, zii] + transf + mod.asset_grid[t];
                                uc = c > 0.0 ? u(c,1.0, γ, mod.ϕ_1,mod.ϕ_2,mod.ϕ_3 ) : u(1e-6,1.0,γ, mod.ϕ_1,mod.ϕ_2,mod.ϕ_3);
                                EVp_grid[aii] = uc *mod.εprobs[t,eii]*mod.zprobs[t,zi,zii] + EVp_grid[aii];
                            end
                        end
                    end
                end

                #EVp_spn = Schumaker(atp,EVp_grid);
				#if constQ == true
				#	  q_spn = qhr;
				#else
					  q_spn = Schumaker(atp,qhr);
			    #end

                @inbounds for εi = 1:N_ε
                    y = mod.Ygrid[t, εi, zi] + transf;

                    #@inbounds Threads.@threads for  ai = 1:N_a #loop over a0
                    for  ai = 1:N_a #loop over a0

                        # With stochastic aging, not a thing:
                        # Right before certain death, the agent chooses to consume all that they can
                        if t == N_ages
                            ms.S[t, tri, ai, εi, zi] = 1.0;
                            ms.B[t, tri, ai, εi, zi] = 0.0;
                            ms.C[t, tri, ai, εi, zi] = a0[ai] + y + transf;
                            ms.V[t, tri, ai, εi, zi] = u(ms.C[t, tri, ai, εi, zi], 1.0, γ,  mod.ϕ_1,mod.ϕ_2,mod.ϕ_3);
                            ms.Ap[t,tri, ai, εi, zi] = 0.0;
                            continue
                        end

                        (vopt,sopt,bopt,copt, Apopt) = solve_ai(mod,y,atp,a0,EVp_grid, #EVp_spn,
                                    N_a,zi,εi,ai,t,ψ, q_spn #qhr #
                                    );

                                    # save the whole frontier of b choices
                        if debug_saves==1 && tri==1
                            for api=1:N_a 
                                bhr =  atp[api] - (1.0 - sopt)*a0[ai]*mod.reclaimrate;
                                ms.B_frontier[t, ai, εi, zi,1,api] = bhr
                                ms.V_frontier[t, ai, εi, zi,1,api] =  bobjective(bhr, sopt, y, q_spn,atp, a0[ai], t, EVp_grid, mod.reclaimrate, ψ, mod.β, mod.ϕ_1, mod.ϕ_2, mod.ϕ_3, γ)
                                # full repay
                                ms.B_frontier[t, ai, εi, zi,2,api] =  atp[api];
                                ms.V_frontier[t, ai, εi, zi,2,api] =  bobjective(atp[api], sopt, y, q_spn,atp, a0[ai], t, EVp_grid, mod.reclaimrate, ψ, mod.β, mod.ϕ_1, mod.ϕ_2, mod.ϕ_3, γ)
                                #full default 
                                if a0[ai] <0
                                    bhr =  atp[api] - a0[ai]*mod.reclaimrate;
                                    ms.B_frontier[t, ai, εi, zi,3, api] = bhr
                                    ms.V_frontier[t, ai, εi, zi,3, api] =  bobjective(bhr, sopt, y, q_spn,atp, a0[ai], t, EVp_grid, mod.reclaimrate, ψ, mod.β, mod.ϕ_1, mod.ϕ_2, mod.ϕ_3, γ)
                                end
                            end
                        end

                        ms.V[t, tri, ai, εi, zi] = vopt;
                        ms.S[t, tri, ai, εi, zi] = sopt;
                        ms.B[t, tri, ai, εi, zi] = bopt;
                        ms.C[t, tri, ai, εi, zi] = copt;
                        ms.Ap[t, tri, ai, εi, zi] = Apopt;


                    end
                    nonmono = 0

                end

            end
        end #tri loop over transfer
    end


    #return V, B, C, S, Ap
    return 1
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
    for agei=1:N_ages
        cumεprobs[agei,2:(N_ε+1)] .= cumsum(mod.εprobs[agei,:]);
        for zi=1:N_z
            cumzprobs[agei,zi,2:(N_z+1)] .= cumsum(mod.zprobs[agei,zi,:]);
        end
        if agei>1
            zergod_agei = sum(mod.zprobs[agei,:,:]^(10000*freq) ,dims=1)/N_z;
            cumzergod[agei,2:(N_z+1)] .= cumsum(zergod_agei[1,:]);
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

function sim_agent!(i::Int,t0::Int,age0::Int, a0::Float64, ei0::Int, zi0::Int, 
    mod::model, ht::hists, ms::sol, mmt::moments,N_ages::Int64, N_a::Int64,N_z::Int64, N_ε::Int64)


    it= t0;
    ei_hr = ei0;
    zi_hr = zi0;
    if age0 > 0 && age0 < N_ages
        age_hr=  age0-1;
    else 
        age_hr=  ht.agehist[i,it];
    end 
    # Initialize
    ht.ahist[i,it] = a0
    
    @inbounds for it=t0:(Thist*freq)
        
        if age0 >= 0 && age0 < N_ages
            age_hr = age_hr + 1 <=N_ages ? age_hr + 1 : 1
        else 
            age_hr=  ht.agehist[i,it];
        end 
        
        aL = 1; aLwt = 0.0;
        #if age_hr > 1
            aL = sum( ht.ahist[i,it] .>=  mod.asset_grid[age_hr,:] );
            aL = aL>= N_a ? N_a-1 : (aL <1 ? 1 : aL);
            aLwt =  (mod.asset_grid[age_hr,aL+1] - ht.ahist[i,it])/(mod.asset_grid[age_hr,aL+1] - mod.asset_grid[age_hr,aL])
        
        if it<=(Thist*freq - 1)
            ap = aLwt*ms.Ap[age_hr,1,aL,ei_hr,zi_hr ] +
                (1.0 - aLwt)*ms.Ap[ age_hr,1,aL+1,ei_hr,zi_hr ];
            ht.ahist[i,it+1] = ap;
            apL = sum( ap .>=  mod.asset_grid[age_hr,:] );
            apL = apL>= N_a ? N_a-1 : (apL <1 ? 1 : apL);
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

        ht.mpshist[1,i,it] = mod.asset_grid[age_hr,aL+1] <= 0 ? (aLwt*ms.S[age_hr,2,aL,ei_hr,zi_hr] +
            (1.0 - aLwt)*ms.S[age_hr,2,aL+1,ei_hr,zi_hr] -
            ht.shist[i,it])/mod.transfer_grid[2] : 0.;
        ht.mpshist[2,i,it] = mod.asset_grid[age_hr,aL+1] <= 0 ? (aLwt*ms.S[age_hr,3,aL,ei_hr,zi_hr] +
            (1.0 - aLwt)*ms.S[age_hr,3,aL+1,ei_hr,zi_hr] -
            ht.shist[i,it])/mod.transfer_grid[3] : 0.0

        if it<=(Thist*freq - 1)
            #a' = b + (1.0 - s) a
            ap_tr2 = aLwt*ms.B[age_hr,2,aL,ei_hr,zi_hr] +
                (1.0 - aLwt)*ms.B[age_hr,2,aL+1,ei_hr,zi_hr] +
                (1.0 - aLwt*ms.S[age_hr,2,aL,ei_hr,zi_hr] -
                (1.0 - aLwt)*ms.S[age_hr,2,aL+1,ei_hr,zi_hr])*ht.ahist[i,it];
            apL = sum( ap_tr2 .>=  mod.asset_grid[age_hr,:] );
            apL = apL>= N_a ? N_a-1 : (apL <1 ? 1 : apL);
            apLwt =  (mod.asset_grid[age_hr,aL+1] - ap_tr2)/(mod.asset_grid[age_hr,aL+1] - mod.asset_grid[age_hr,aL])
            q_tr2 = apLwt*ms.Q0[age_hr,zi_hr,apL] +
                (1.0 - apLwt)*ms.Q0[age_hr,zi_hr,apL+1];
            ht.mpqhist[1,i,it] = (q_tr2 - ht.qhist[i,it])/mod.transfer_grid[2];

            ap_tr3 = aLwt*ms.B[age_hr,3,aL,ei_hr,zi_hr] +
                (1.0 - aLwt)*ms.B[age_hr,3,aL+1,ei_hr,zi_hr] +
                (1.0 - aLwt*ms.S[age_hr,3,aL,ei_hr,zi_hr] -
                (1.0 - aLwt)*ms.S[age_hr,3,aL+1,ei_hr,zi_hr])*ht.ahist[i,it];
            apL = sum( ap_tr3 .>=  mod.asset_grid[age_hr,:] );
            apL = apL>= N_a ? N_a-1 : (apL <1 ? 1 : apL);
            apLwt =  (mod.asset_grid[age_hr,aL+1] - ap_tr3)/(mod.asset_grid[age_hr,aL+1] - mod.asset_grid[age_hr,aL])
            q_tr3 = apLwt*ms.Q0[age_hr,zi_hr,apL] +
                (1.0 - apLwt)*ms.Q0[age_hr,zi_hr,apL+1];
            ht.mpqhist[2,i,it] = (q_tr3 - ht.qhist[i,it])/mod.transfer_grid[3]
        end
        
    end #end it=1:Tsim+burnin

end


function sim_hists!(mod::model, ht::hists, ms::sol, mmt::moments,N_ages::Int64, N_a::Int64,N_z::Int64, N_ε::Int64)

    ithist = repeat( collect(1:Thist*freq)',Nsim,1 ) ;
    Tsim0 = (Tburnin+1)*freq;
    TsimT = (Tsim+Tburnin)*freq;
    age_z_mean = zeros(N_ages,N_z);
    Nhere = zeros(N_ages,N_z);
    for dist_iter = 1:4
        @inbounds for i=1:Nsim
            age_hr = ht.agehist[i,1];
            ei_hr = ht.εidxhist[i,1];
            zi_hr = ht.zidxhist[i,1];
            #a95_hr = a95[it,ei_hr,zi_hr];
            if age_hr >1 && dist_iter>1
                idraw = max(min(floor(Int,ht.a0qtlhist[i]*(Nsim+1))+1,Nsim),1);
                tdraw = max(min(floor(Int,i*Tsim*freq + Tburnin*freq )+ 1,Tsim*freq),1) ;
                age_draw = ht.agehist[idraw,tdraw]
                zi_draw   = ht.zidxhist[idraw,tdraw]
                a0 = age_z_mean[age_hr,zi_hr] + (ht.ahist[age_draw,zi_draw] - age_z_mean[age_draw,zi_draw])
            elseif age_hr==1 #young all go to the same place
                a0 = 0.0;
            else
                rand_i = max(min(floor(Int,ht.a0qtlhist[i]*(N_a+1)),N_a ), 1)
                a0 = mod.asset_grid[age_hr,rand_i];
                #a0 = 0.0;
            end
            sim_agent!(i,1,-1, a0, ei_hr, zi_hr, mod, ht, ms, mmt,N_ages, N_a,N_z, N_ε)
            
        end #end i =1:Nsim

        # update asset_alpha, asset_beta
        for t=2:N_ages
            for zi=1:N_z
                mask_hr = (ht.zidxhist.==zi) .& (ht.agehist .>= t-2) .& (ht.agehist .<= t+2) .& (ithist .> (Tburnin*freq));
                Nhere[t,zi] = sum(mask_hr) 
                if(sum(mask_hr)>0)
                    age_z_mean[t,zi] = mean(ht.ahist[mask_hr]);
                else
                    age_z_mean[t,zi]= age_z_mean[t-1,zi] ;
                end
                
            end
        end
        avg_avg = mean(age_z_mean);
        println("average asset is $avg_avg")
    end #iterated simulations

    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #after the simulation, compute moments based on the simulation
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    #+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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


    MPCs_negassets = ones(Nsim*(TsimT-Tsim0+1),6)
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
                MPCs_negassets[idx_neg,2] =  (ht.mpchist[1,i,it] + ht.mpbhist[2,i,it])/2;
                MPCs_negassets[idx_neg,3] =  (ht.mpbhist[1,i,it] + ht.mpbhist[2,i,it])/2;
                MPCs_negassets[idx_neg,4] =  (ht.mpshist[1,i,it] + ht.mpshist[2,i,it])/2;
                MPCs_negassets[idx_neg,5] =  (ht.mpqhist[1,i,it] + ht.mpqhist[2,i,it])/2;
                MPCs_negassets[idx_neg,6] =  (ht.mpshist[1,i,it]+ht.mpshist[2,i,it])/2*(-ht.ahist[i,it])+ (ht.mpbhist[1,i,it]+ht.mpbhist[2,i,it])/2 + (ht.mpqhist[1,i,it]+ht.mpqhist[2,i,it])/2;
            else
                denom_assets = 1.0 + denom_assets;
                mmt.mpc_ac0 = (ht.mpchist[1,i,it] + ht.mpchist[2,i,it]) /2 + mmt.mpc_ac0;
                idx_pos += 1;
                MPCs_posassets[idx_pos,1] = ht.ahist[i,it] ;
                MPCs_posassets[idx_pos,2] = (ht.mpchist[1,i,it]+ht.mpchist[2,i,it])/2 ;
                MPCs_posassets[idx_pos,3] = (ht.mpbhist[1,i,it]+ht.mpbhist[2,i,it])/2 ;


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
    if idx_neg>0
        mmt.mpc_ac0 = mmt.mpc_ac0 / denom_assets;
        mmt.mpc_cor_neg_asset = cor(MPCs_negassets[1:idx_neg,:]);
        mmt.mpc_cor_pos_asset = cor(MPCs_posassets[1:idx_pos,:]);
    else
        mmt.mpc_ac0
    end
    mmt.fr_neg_asset = mmt.fr_neg_asset/(mmt.fr_neg_asset+ denom_assets);

    println("The average asset position is $Ea_tot and pay-back is $Es_tot");

end


function steadystateDist!( ms::sol, mod::model, survivalprobs::Vector{Float64}, N_ages, N_a, N_z, N_ε)

    survivalprobs_cml = sum(1.0 .- survivalprobs[1:(N_ages-1)]) #fraction dead before period N_ages+1
    survivalprobs[N_ages] = 0.0;
    N_arnd = Int64(N_a/5);
    asteps = Int64(N_a/N_arnd);
    asset_grid_rnd = zeros(N_ages,N_arnd);
    asset_space    = zeros(N_ages,N_arnd);
    for t=1:(N_ages-1)
        for airnd=1:N_arnd
            for ai=1:asteps
                if (airnd-1)*asteps + ai+1 < N_a
                    asset_space[t,airnd]+= mod.asset_grid[t,(airnd-1)*asteps + ai+1] - mod.asset_grid[t,(airnd-1)*asteps + ai];
                end
                asset_grid_rnd[t,airnd] += mod.asset_grid[t,(airnd-1)*asteps + ai]/asteps;
            end
        end
        # need to make sure the ends are the same
        asset_grid_rnd[t,1] = mod.asset_grid[t,1]
        asset_grid_rnd[t,N_arnd] = mod.asset_grid[t,N_a]
    end
    Nstate = N_ages*N_arnd*N_ε*N_z;
    a0idx = findlast(  asset_grid_rnd[2,:] .<= 0.0 );  #mod.asset_grid[2,:] .<= 0.0)
    # state is (t-1)*N_a*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi
    flatstatemat = spzeros(Float64,Nstate,Nstate);
    for t=1:N_ages
        for ai=1:N_arnd
            for ei=1:N_ε
                for zi=1:N_z
                    #apL   = ms.Ap[t,ai,ei,zi]>  mod.asset_grid[t+1,1] ? findlast( ms.Ap[t,ai,ei,zi] .>  mod.asset_grid[t+1,:]) : 1;
                    #apL   = isnan(apL) | isnothing(apL) ? N_a : apL;
                    #apLwt =  apL < N_a ? ( mod.asset_grid[t+1,apL+1] - ms.Ap[t,ai,ei,zi] )/( mod.asset_grid[t+1,apL+1] -  mod.asset_grid[t+1,apL]) : 1.0 ;
                    Ap_rnd = 0.0;
                    for aihr=1:asteps
                        Ap_rnd = ms.Ap[t,(ai-1)*asteps+aihr,ei,zi]/asteps + Ap_rnd;
                    end
                    apL   = Ap_rnd >  asset_grid_rnd[t+1,1] ? sum( Ap_rnd .>=  asset_grid_rnd[t+1,:]) : 1;
                    apL   = isnan(apL) | isnothing(apL) ? N_arnd : apL;
                    apLwt =  apL < N_arnd ? ( asset_grid_rnd[t+1,apL+1] - Ap_rnd )/( asset_grid_rnd[t+1,apL+1] -  asset_grid_rnd[t+1,apL]) : 1.0 ;
                    if t<(N_ages-1)
                        for zzi=1:N_z
                            for eei=1:N_ε
                                flatstatemat[ (t-1)*N_arnd*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi, (t)*N_arnd*N_z*N_ε + (apL-1)*N_ε*N_z + (eei-1)*N_z + zzi] =
                                    apLwt* mod.zprobs[t,zi,zzi]* mod.εprobs[t,eei] * survivalprobs[t] +
                                    flatstatemat[ (t-1)*N_arnd*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi, (t)*N_arnd*N_z*N_ε + (apL-1)*N_ε*N_z + (eei-1)*N_z + zzi];
                                if apL<N_arnd
                                    flatstatemat[ (t-1)*N_arnd*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi, (t)*N_arnd*N_z*N_ε + apL*N_ε*N_z + (eei-1)*N_z + zzi] =
                                        (1.0 - apLwt)*mod.zprobs[t,zi,zzi]* mod.εprobs[t,eei] * survivalprobs[t] +
                                        flatstatemat[ (t-1)*N_arnd*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi, (t)*N_arnd*N_z*N_ε + apL*N_ε*N_z + (eei-1)*N_z + zzi];
                                end
                            end
                        end
                    else
                        #everyone dies and goes to same assets in ag 0
                        flatstatemat[ (t-1)*N_arnd*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi, 0*N_arnd*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi] = 1.0;
                    end
                end
            end
        end
    end

    Npara = N_ages*N_arnd*N_ε*N_z;
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
        for ai=1:N_arnd
            for ei=1:N_ε
                for zi=1:N_z
                    dist_hr = real(flatstatevecs[ (t-1)*N_arnd*N_z*N_ε + (ai-1)*N_ε*N_z + (ei-1)*N_z + zi ,1]);
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
