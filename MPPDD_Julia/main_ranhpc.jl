# run line: JULIA_NUM_THREADS=32 julia main.jl
# on RAN run with:  julia150-batch 5 1 main_ranhpc.jl

## Packages

#using Pkg
using Distributed
#when running on the RAN cluster
using ClusterManagers

const nworkers = 128;
#when running on the RAN cluster
mem      = 5;  # request memory for each worker
ENV["frbnyjuliamemory"]=string(mem)*"G"
addprocs_frbny(nworkers);


@everywhere begin 
    using Distributed
    using LinearAlgebra
    #using SparseArrays, ParSpMatVec
    #using Arpack
    using SchumakerSpline, Interpolations
    using Distributions, Random
    using Optim
    using JLD
    #using Plots
    using DataFrames
    using DataStructures
    using SharedArrays
    using DelimitedFiles
    using Parameters

end

@everywhere begin #: bcast to all workers, split from the 'using' packages b/c Julia is quirky
    #import Base.Threads.@spawn
    #cd("/Users/meghanagaur/DropBox/Consumption_Model_Codes/julia/")     # personal
    #cd("C:\\Users\\rcemxg23\\DropBox\\Consumption_Model_Codes\\MPPDD_julia") # ran virtual desktop
    #cd("/home/david/workspace/rebate_timing/MPPDD_Julia/")
    #cd("/Users/rcedxm19/Dropbox/MPPDD/codes/MPPDD_Julia")

    #set root directory for saving output here:
    saveroot = pwd()*"/"#"/home/david/workspace/rebate_timing/MPPDD_Julia/"

    mutable struct sol
        V::Array{Float64,5}    # defined for (N_ages+1,N_t,N_k,N_ε, N_z)
        C::Array{Float64,5}    # defined for (N_ages,N_t,N_k,N_ε, N_z)
        S::Array{Float64,5}    # defined for (N_ages,N_t,N_k,N_ε, N_z)
        B::Array{Float64,5}    # defined for (N_ages,N_t,N_k,N_ε, N_z)
        Ap::Array{Float64,5}   # defined for (N_ages,N_t,N_k,N_ε, N_z)
        #Vtrmbl::Array{Float64,6}# defined for (N_ages,N_t,N_k,N_ε, N_z)
        #Aptrmbl::Array{Float64,6}# defined for (N_ages,N_t,N_k,N_ε, N_z)
        Q0::Array{Float64,3}   # defined for (N_ages,N_k,N_z) known states
        
        Aεz_dist::Array{Float64,4} #Equilibrium distribution across state (N_ages,a,z,ε)
        sol() = new()
    end


    mutable struct model
        zprobs::Array{Float64,3}
        zgrid ::Array{Float64,2}
        εprobs::Array{Float64,2}
        εgrid ::Array{Float64,2}
        Ygrid::Array{Float64,3}
        asset_grid::Array{Float64,2}
        cbar::Float64
        γ::Int64
        ϕ_1::Float64
        ϕ_2::Float64
        r::Float64
        β::Float64
        λ::Float64
        transfer_grid::Array{Float64,1}

        model() = (modd = new(); modd.cbar = 0.0; modd.γ =1; return modd)
    end

    mutable struct hists
        #histories for:
        chist::Array{Float64,2}
        bhist::Array{Float64,2}
        shist::Array{Float64,2}
        ahist::Array{Float64,2}
        qhist::Array{Float64,2}
        #should be drawn once:
        εidxhist::Array{Int64,2}
        zidxhist::Array{Int64,2}
        agehist::Array{Int64,2}
        incomehist::Array{Float64,2}
        mpchist::Array{Float64,3}
        mpshist::Array{Float64,3}
        mpbhist::Array{Float64,3}
        mpqhist::Array{Float64,3}

        a0qtlhist::Array{Float64,1} #quantile of a in the asset distribution for first period.

        hists() = new()
    end
    
    mutable struct moments
        #moments from the simulation to return
        fr_neg_asset::Float64
        avg_indiv_debt_income::Float64
        avg_debt::Float64
        avg_income::Float64
        avg_income_neg_asset::Float64
        extent_delinq::Float64
        avg_s::Float64
        fr_s1::Float64
        mpc_cor_neg_asset::Array{Float64,2}
        mpc_cor_pos_asset::Array{Float64,2}
        mpc_ac0 :: Float64
    
    end 

    function moments()
        #constructor for moments
        fr_neg_asset = 0.0
        avg_indiv_debt_income =0.0
        avg_debt =0.0
        avg_income =0.0
        avg_income_neg_asset =0.0
        extent_delinq =0.0
        avg_s =0.0
        fr_s1 =0.0
        mpc_cor_neg_asset = zeros(5,5)
        mpc_cor_pos_asset = zeros(3,3)
        mpc_ac0 = 0.0

        return moments(fr_neg_asset,avg_indiv_debt_income,avg_debt,avg_income,avg_income,extent_delinq,avg_s,fr_s1,mpc_cor_neg_asset,mpc_cor_pos_asset,mpc_ac0)

    end


    include("functions/incomeProcess.jl")
    include("functions/agent.jl")

    const qtol         = 1e-3; # tolerance for convergence in Q
    const maxiterq     = 40;   # how many iterations?
    const maxiterV     = 30;
    const Vtol         = 1e-3; #can be reasonably loose b/c iterating on Q too
    ## Parameters
    const freq         = 4; # 1=annual, 4=quarterly
    const γ            = 1;


    ## State space

    const N_k          = 50;    # Number of asset states, 40-50
    const N_z          = 7;      # Number of permanent income states, 7-10
    const N_ε          = 5;         # Number of transitory income states, 5-21
    const N_t          = 3;      # Number of transfer levels

    ## Income Process  - from Kaplan & Violante (2010)
    ρ          =  1.0;
    σ_η       = (0.01/freq)^(1/2);  # Variance of permanent income shocks = 0.01
    σ_z0      = (0.15/freq)^(1/2);  # Variance of permanent income shocks = 0.15
    σ_ε       = (0.05/freq)^(1/2);  # Variance of transitory income shocks = 0.05

    ## Simulation:
    const Nsim      = 20000;  #number of Agents
    const Tsim      = 12;     #number of periods in the simulation
    const Tburnin   = 32;     #number of burnin periods
    const Thist     = Tsim+Tburnin;

    const Nparams   = 5; #calibration parameters
    const Nmoments  = 8; # calibration moments

    ## Ages:

    const YrsT = 85;
    const Yrs_retire = 65; #corresponds to 65
    const Yrs0 = 25;
    const N_ages = (YrsT-Yrs0+1)*freq;
    const age_retire = (Yrs_retire-Yrs0+1) * freq;
    const ageprobs = ones(N_ages+1); #vcat(ones(N_ages-1)*(N_ages-1.0)/(T_retire - T0)/freq, 1.0/(T-T_retire)/freq );

    ## Agents begin at  age 25, die at 95, retire at 65

    # Survival Probabilities
    const survivalprobs = ones(N_ages); # N_ages set to 0 in agent.jl


    const ϕ_2grid = [2.0 3.0];
    const ϕ_1grid = [0.2 0.8];
    const βgrid   = [0.9 0.95];


    function alloc_hists!(ht::hists)
        ht.a0qtlhist = zeros(Nsim);
        ht.zidxhist  = Array{Int64,2}(undef, Nsim, Thist*freq);
        ht.εidxhist  = Array{Int64,2}(undef, Nsim, Thist*freq);
        ht.ahist     = zeros(Nsim,Thist*freq);
        ht.shist     = zeros(Nsim,Thist*freq);
        ht.chist     = zeros(Nsim,Thist*freq);
        ht.bhist     = zeros(Nsim,Thist*freq);
        ht.qhist     = zeros(Nsim,Thist*freq);
        ht.incomehist= zeros(Nsim,Thist*freq);
        ht.agehist   = Array{Int64,2}(undef, Nsim, Thist*freq);
        ht.mpchist   = Array{Float64,3}(undef,2,Nsim,Thist*freq);
        ht.mpshist   = Array{Float64,3}(undef,2,Nsim,Thist*freq);
        ht.mpbhist   = Array{Float64,3}(undef,2,Nsim,Thist*freq);
        ht.mpqhist   = Array{Float64,3}(undef,2,Nsim,Thist*freq);

    end 

    function alloc_sol!(ms::sol)
        ms.V  = zeros(N_ages+1,N_t,N_k,N_ε, N_z);
        ms.C  = zeros(N_ages,N_t,N_k,N_ε, N_z);
        ms.S  = zeros(N_ages,N_t,N_k,N_ε, N_z);
        ms.B  = zeros(N_ages,N_t,N_k,N_ε, N_z);
        ms.Ap = zeros(N_ages,N_t,N_k,N_ε, N_z);
        #ms.Vtrmbl  = zeros(N_ages,N_t,N_k,N_ε, N_z, 2*N_trmbl+1);
        #ms.Aptrmbl = zeros(N_ages,N_t,N_k,N_ε, N_z, 2*N_trmbl+1);
    end 


    ## Main function to solve model
    function solveQ!(mod::model,ms::sol,ht::hists,mmt::moments)

        #mod communicates parameters, etc. ht needs to be passed in to keep the random draws across parameters. Also contains simulation histories of endogenous variables.

        updq = 0.5;  #speed to update the price function

        ## Initialize asset and repayment grid
        asset_max  = zeros(N_ages+1);
        asset_min  = zeros(N_ages+1);
        mod.asset_grid = zeros(N_ages+1, N_k);
         # Minimum assets (Maximum debt)
        for t = N_ages:-1:1
            income_min   = minimum(mod.Ygrid[t,:,:]);
            asset_min[t] = t>= age_retire ? income_min : asset_min[t+1]/(1.0+mod.r) - income_min+ 1e-4;
        end

        for t = 1:N_ages
            income_max     = maximum(mod.Ygrid[t,:,:]);
            asset_max[t+1] = min(asset_max[t] + 0.95*income_max, 100)
        end
        asset_max[1] = asset_max[2];

        for t = 1:N_ages
            mod.asset_grid[t, :] = LinRange(asset_min[t],asset_max[t],N_k);
        end

        mod.transfer_grid = zeros(N_t);
        #median quarterly income: ~ 3.0, in data it's 68,703/4 = 17175.75
        # => 1200/17175.75
        mod.transfer_grid[2] = 0.21; mod.transfer_grid[3] = 0.419;

        ## Initialize Q0
        ms.Q0 = ones((N_ages,N_z,N_k)).*(1/(1+mod.r));
        for t=1:N_ages
            #a0 = sum( asset_grid[t,:].<0.0 )+1
            for i = 1:N_z
                for api=1:N_k
                    if t<age_retire
                        ms.Q0[t, i, api] = mod.asset_grid[t,api]<0.0 ? exp(.1*mod.asset_grid[t,api]) / (1.0 + mod.r) : 1.0 / (1.0 + mod.r);
                    else 
                        ms.Q0[t, i, api] = 1.0/(1.0 + mod.r);
                    end
                end
            end
        end


         ϕ_2 = mod.ϕ_2
        println("Starting to compute eqm for ϕ_2 = $ϕ_2")
        println("")

        #++++++++++++++++++++++++++++++++++++++++++
        # use last Q as a guess or go back to our initial one?
        for qiter = 1:40

            ## Asset and repayment grid : note, this depends on Q
            # stop updating this after ~5 iterations:
            if qiter>1
                asset_grid_hr = assetGrid(N_ages, mod.Ygrid, N_k, ms.Q0, mod.asset_grid, mod.r);
                mod.asset_grid .= (1.0 - updq) .* mod.asset_grid .+ updq .* asset_grid_hr;
                mina = minimum(mod.asset_grid);
                for t=1:N_ages
                    sort!(mod.asset_grid[t,:]);
                end 
                println("On iter $qiter, most possible borrowing: $mina")
            else 
               println("Iter 1")
            end 
            #println("$mod.asset_grid[N_ages,:]");
            
            @time backwardSolve!(ms,mod, N_ages, N_k, N_z, N_ε, ms.Q0);

            Q = equilibriumQ(ms.Q0, ms, mod.asset_grid, mod.r, mod.λ, mod.εprobs, mod.zprobs);
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
                for k=1:N_k
                    for εi =1:N_ε
                        for zi =1:N_z
                            sinterior = ms.S[t, 1, k, εi, zi]< 0.99 && ms.S[t, 1, k, εi, zi]> 0.01 ? sinterior + 1 : sinterior
                            sbad0     = ms.S[t, 1, k, εi, zi]< 0.99 && mod.asset_grid[t, k]>0.0 ? sbad0 + 1     : sbad0
                            sbad0T    = ms.S[t, 1, k, εi, zi]< 0.99 && mod.asset_grid[t, k]>0.0 && t< N_ages-1 ? sbad0T + 1     : sbad0T
                        end
                    end
                end
            end
            #println("The number of bad0 payment shares: $sbad0, $sbad0T were young . Number of interior payment shares: $sinterior")


            #this seems to be an allocation/creation operation rather than element-wise assignment
            # Q0 = deepcopy(Q);
            for t=1:N_ages
                for zi = 1:N_z
                    for api = 1:N_k
                        ms.Q0[t,zi,api] = updq * Q[t,zi,api] + (1.0 - updq) * ms.Q0[t,zi,api];
                    end
                end
            end
            if Qdist < qtol
                break;
            end

        end

        return ms.Q0;

    end #end solveQ!

    """
    Takes in parameters, outputs moments / parameters to CSV 
    """
    function params_moments(iofilename, vals, j)
    

        # FILL IN MAIN STUFF HERE
        params   = vals[j];
        age_peak       = (Yrs_retire-10)*freq; # corresponds to 50-55, at which income peaks

        ht = hists();
        alloc_hists!(ht);

        mod      = model();
        mod.r    = params[:r];
        mod.β    = params[:β];
        mod.ϕ_1  = params[:ϕ_1];
        mod.ϕ_2  = params[:ϕ_2];
        mod.λ  = params[:λ];

        mmt = moments();
        mod.Ygrid, mod.zgrid, mod.εgrid = get_Ygrid(age_retire, age_peak, N_ages, N_z, N_ε, σ_ε, σ_z0, σ_η, ρ); # level grids of income

        mod.zprobs    = get_zprobs(N_ages, N_z, σ_z0, σ_η, ρ); # transition matrices for permanent income
        mod.εprobs    = get_εprobs(N_ages, N_ε, σ_ε);          # transition matrices for transitory income

        ms = sol();
        alloc_sol!(ms);
        pidnum = getpid();
        workernum = myid() - 1 ;

        println("Drawing shocks on PID $pidnum, $workernum at j $j !");
        draw_shocks!(mod,ht, 12281951);

        println("Solving on PID $pidnum, $workernum at j $j !");
        # optionally, test the timing of this to see how things output
        #sleep(ceil(Int,rand()*5))
        ms.Q0 = solveQ!(mod,ms,ht,mmt) 
        
        println("Simulating!  $pidnum, $workernum at j $j");
        sim_hists!(mod, ht, ms, mmt, N_ages, N_k,N_z, N_ε);
        
        saveloc = string(saveroot,"modEnvr_p",j,".jld")
        @save saveloc mod
        #saveloc = string(saveroot,"simHists_p",j,".jld")
        #@save saveloc ht
        saveloc = string(saveroot,"eqmQ_p",j,".jld")
        @save saveloc ms.Q0
        saveloc = string(saveroot,"solMats_p",j,".jld")
        @save saveloc ms
        

        outvec = reshape([j, mod.r, mod.β, mod.ϕ_1, mod.ϕ_2, mod.λ , mmt.fr_neg_asset, mmt.avg_debt,mmt.avg_s, mmt.extent_delinq,mmt.fr_s1, mmt.avg_income,mmt.avg_income_neg_asset,mmt.mpc_ac0],1,Nparams+Nmoments+1 );
        iohr = open(iofilename,"a")
        writedlm(iohr, round.(outvec,digits=5), ',')
        close(iohr);


        iohr = open("corMat_negasset_p", j,"a");
        writedlm(iohr, round.(corMat,digits=5), ',');
        close(iohr);
        
        
        iohr = open("corMat_posasset_p", j,"a");
        writedlm(iohr, round.(corMat,digits=5), ',');
        close(iohr);
    
        #sharedmom[j,1] = mmt.avg_indiv_debt_income;
        #sharedmom[j,2] = mmt.avg_s;
        #sharedmom[j,3] = mmt.fr_neg_asset;
        #sharedmom[j,4] = mmt.fr_s1;
    end
    
    
end # end parallel everywhere block




function estimate_elasticities(param_grids::OrderedDict{Symbol,Array{Float64,1}}, file::String = string(saveroot,"out.csv"))

    file =string(saveroot,"explore_params.csv");

    params_header = reshape(["idx_num","r","beta","phi_1", "phi_2", "borrow_wdg","fr_neg_asset", "avg_debt", "avg_s","extent_delinq", "fr_s1", "avg_income","avg_income_neg_asset","MCP_a0"],1,Nparams+Nmoments+1 );
    
    # note: this will 'delete/overwrite' the file if it already exists
    writedlm(file,params_header,',') # note: last argument is the delimiter which should be the same as below

    parvals = OrderedDict{Int64, OrderedDict{Symbol,Float64}}()
    i =1;
    for ω in values(param_grids[:ω]), β in values(param_grids[:β]), ϕ_1 in values(param_grids[:ϕ_1]), ϕ_2 in values(param_grids[:ϕ_2]), λ in values(param_grids[:λ])
        
        β_f = β^(1/freq); # correct for frequency of model
        temp = OrderedDict{Symbol, Float64}(
                    :β => β_f,
                    :ϕ_1 => ϕ_1/freq,
                    :ϕ_2 => ϕ_2,
                    :r =>   (1/β_f-1)*(1-ω),
                    :λ =>   λ)
        push!(parvals, i => temp)
        i+=1
    end

    
    # initialize some shared arrays in case we want to export
    #momentsout  = SharedArray{Float64}(length(parvals),Nmoments);

    # better for heavy lifting, allow julia to distribute
    pmap(j->params_moments(file, parvals, j), 1:length(parvals))

    # Save a final output (note: jld2 seems to work better for shared arrays)
    # with additional objects
    #saveloc = string(saveroot,"MPPDD_test.jld2")
    #@save saveloc momentsout
end

# parameter grids, annual values (where relevant)
param_grids  = OrderedDict{Symbol,Array{Float64,1}}([
        (:ω, collect(LinRange(0.01, 0.99, 4)) ),
        (:β, [0.90, 0.95, 0.99, 0.995]),
        (:ϕ_1,[0.01, 0.1, 0.2, 0.3]),
        (:ϕ_2, [1.25, 2.0, 4, 5.5]),
        (:λ , collect(LinRange(0.0, 0.06, 4))) ]);

estimate_elasticities(param_grids)

rmprocs(workers())
