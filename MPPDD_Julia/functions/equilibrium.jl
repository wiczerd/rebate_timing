

## Main function to solve model
function solveQ!(mod::model,ms::sol,ht::hists,mmt::moments)

    #mod communicates parameters, etc. ht needs to be passed in to keep the random draws across parameters. Also contains simulation histories of endogenous variables.

    updq = 0.5;  #speed to update the price function

    
    ## Initialize Q0
    ms.Q0 = ones(N_ages,N_z,N_a).*(1.0/(1.0+mod.r));
    ## Initialize asset grid
    for t=1:N_ages
        mod.asset_grid[t,:] .= collect(LinRange(1,N_a,N_a));
    end
    mod.asset_grid = assetGrid(N_ages,mod.Ygrid,N_a,ms.Q0, mod.asset_grid,mod.r, grid_curvature );

    #=  Another Q initialization? 
    Tried to initialize with some delinquency built in
    for t=1:N_ages
        #a0 = sum( asset_grid[t,:].<0.0 )+1
        for i = 1:N_z
            for api=1:N_a
                if t<age_retire
                    ms.Q0[t, i, api] = mod.asset_grid[t,api]<0.0 ? exp(.1*mod.asset_grid[t,api]) / (1.0 + mod.r) : 1.0 / (1.0 + mod.r);
                else 
                    ms.Q0[t, i, api] = 1.0/(1.0 + mod.r);
                end
            end
        end
    end
    =#
    mod.transfer_grid = zeros(N_t);
    #median quarterly income: ~ 3.0, in data it's 68,703/4 = 17175.75
    # => 1200/17175.75
    mod.transfer_grid[2] = 0.21; mod.transfer_grid[3] = 0.419;




    #++++++++++++++++++++++++++++++++++++++++++
    # use last Q as a guess or go back to our initial one?
    rlow  =-0.025;
    rhigh = 2*mod.r - rlow;


    for qiter = 1:maxiterq
       
        @time backwardSolve!(ms,mod, N_ages, N_a, N_z, N_ε, ms.Q0);
        excessdemand =0.0;
        
        if fullcommit == false && constQ == false
            Q = equilibriumQ(ms.Q0, ms, mod);

            Qdist = maximum( abs.(Q .- ms.Q0) );
            loc=CartesianIndex{2}(1,1);
            maxageQ =0;maxage=1;
            for t=1:N_ages
                if maxageQ < maximum( abs.(Q[t,:,:]  .- ms.Q0[t,:,:]) );
                    maxageQ,loc = findmax( abs.(Q[t,:,:]  .- ms.Q0[t,:,:]) );
                    maxage = t;

                end
            end
            println("sup norm of distance is $Qdist at max age: $maxage and z,a $loc");
            
            # little diagnostic:
            sinterior = 0
            sbad0 = 0
            sbad0T = 0
            for t=1:N_ages
                for k=1:N_a
                    for εi =1:N_ε
                        for zi =1:N_z
                            sinterior = ms.S[t, 1, k, εi, zi]< 0.99 && ms.S[t, 1, k, εi, zi]> 0.01 ? sinterior + 1 : sinterior
                            sbad0     = ms.S[t, 1, k, εi, zi]< 0.99 && mod.asset_grid[t, k]>0.0 ? sbad0 + 1     : sbad0
                            sbad0T    = ms.S[t, 1, k, εi, zi]< 0.99 && mod.asset_grid[t, k]>0.0 && t< N_ages-1 ? sbad0T + 1     : sbad0T
                        end
                    end
                end
            end
            println("The number of bad0 payment shares: $sbad0, $sbad0T were young . Number of interior payment shares: $sinterior")
        else 
            mmt_i = moments();
            sim_hists!(mod,ht,ms,mmt_i,N_ages,N_a,N_z,N_ε);
            excessdemand = mean(ht.ahist);
            println("Excess demand is $excessdemand")
            if excessdemand >0
                rhigh = mod.r;
            else
                rlow = mod.r;
            end 
            mod.r = 0.5*(rhigh + rlow);
        end
      
        #this seems to be an allocation/creation operation rather than element-wise assignment
        # Q0 = deepcopy(Q);
        if fullcommit == false && constQ==false
            for t=1:N_ages
                for zi = 1:N_z
                    for api = 1:N_a
                        ms.Q0[t,zi,api] = updq * Q[t,zi,api] + (1.0 - updq) * ms.Q0[t,zi,api];
                    end
                end
            end
        else 
            ms.Q0 .= 1.0/(1.0 + mod.r);
            Qdist  = rhigh-rlow;
            println("New interest rate is $(mod.r) and bounds: $rlow,$rhigh")
            
        end 
        
        
        if Qdist < qtol
            break;
        # asset grid depends on Q, so update as I go
        # update the grid a bit more slowly than the Q 
        # can't keep updating it 
        elseif  fixedassetgrid==false && fullcommit==false && constQ==false
            
            asset_grid_new = assetGrid(N_ages, mod.Ygrid, N_a, ms.Q0, mod.asset_grid, mod.r, grid_curvature);
            #asset_grid_new = (1.0 - updq/2) .* mod.asset_grid .+ updq/2 .* asset_grid_hr;
            for t=1:(N_ages-1)
                #min_agrid_old[t] = minimum(mod.asset_grid[t,:]);
                #min_agrid_new2[t] = minimum(asset_grid_new[t,:]);

                for zi = 1:N_z
                    for api = 1:N_a
                        Q[t,zi,api] = max( LinearInterpolation(mod.asset_grid[t+1,:], ms.Q0[t,zi,:]; extrapolation_bc=Line())(asset_grid_new[t+1,api]) ,1e-4 ) ;
                    end
                end
            end
            
            ms.Q0 .= Q;
            mod.asset_grid .= asset_grid_new;
            mina = minimum(mod.asset_grid);
            
            println("On iter $qiter, most possible borrowing: $mina")
        elseif fullcommit ==true || constQ==true
            println("Making new asset grid on iteration $qiter")
            asset_grid_new = assetGrid(N_ages, mod.Ygrid, N_a, ms.Q0, mod.asset_grid, mod.r, grid_curvature);
            mod.asset_grid .= asset_grid_new;
            mina = minimum(mod.asset_grid);
            
            println("On iter $qiter, most possible borrowing: $mina")
        else 
            println("Q iteration $qiter")
        end 
        
    end

    return ms.Q0;

end #end solveQ!


"""
Solve for implied Q, given S
"""
function equilibriumQ(Q0::Array{Float64,3}, ms::sol, mod::model)
    #  Dimension is N_ages, z, ap, b
    qimplied = similar(Q0);
    r = mod.r;
    asset_grid = mod.asset_grid;
    λ = mod.λ;
    εprobs = mod.εprobs;
    zprobs = mod.zprobs;

    
    for t = N_ages:-1:1
        #if t >= T - 2
        #    qimplied[t, :, :] .= 1.0 / (1.0 + r + I_{b>0} )
        #    continue
        #end

        for zi = 1:N_z
            for ai = 1:N_a
                qimplied[t, zi, ai] = 0.0;
                for εi =1:N_ε
                    for zii=1:N_z
                        #=if age:
                        # what to do with old dead guys? Replaced by young people exactly the same as them
                        PrAge_tp1 = ageprobs[t]

                        app_tp1 = ms.Ap[tp1, 1, ai, εi, zii];
                        app_tp1 = app_tp1 < asset_grid[tp1,1]   ? asset_grid[tp1,1] : app_tp1
                        app_tp1 = app_tp1 > asset_grid[tp1,N_a] ? asset_grid[tp1,N_a] : app_tp1
                        
                        shr_tp1 = t < N_ages ? ms.S[tp1, 1, ai, εi, zii] : ms.S[t, 1, ai, εi, zii];
                        qp_tp1  = LinearInterpolation(asset_grid[tp1,:],Q0[tp1,zii,:])(app_tp1);
                        =#
                        tp1 = t < N_ages ? t+1 : t;
                        tp2 = t+1 < N_ages ? t+2 : tp1;
                        
                        apL = sum(ms.Ap[tp1, 1, ai, εi, zii] .> asset_grid[tp1,:]) ;
                        apL = apL+1 > N_a ? N_a-1 : apL;
                        apLwt = (asset_grid[tp2,apL+1] - ms.Ap[tp1, 1, ai, εi, zii] )/(asset_grid[tp2,apL+1] - asset_grid[tp2,apL]) ;

                        app = apLwt *ms.Ap[tp1, 1, apL, εi, zii] + (1.0 - apLwt)* ms.Ap[tp1, 1, apL+1, εi, zii]; 
                        app = app < asset_grid[tp2,1]   ? asset_grid[tp2,1] : app #force in bounds (this shouldn't bind)
                        app = app > asset_grid[tp2,N_a] ? asset_grid[tp2,N_a] : app
                        shr = ms.S[tp1, 1, ai, εi, zii];
                    
                        qp      = LinearInterpolation(asset_grid[tp2 ,:],Q0[tp1  ,zii,:])(app);
                        qimplied[t, zi, ai] =  ( shr + (1-shr)*qp*mod.reclaimrate
                            )*εprobs[t,εi]*zprobs[t,zi,zii] +
                            qimplied[t, zi, ai] ;  #risk premium on loans. Will discount by risk-free below:

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
                # discount by the risk-free rate
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
