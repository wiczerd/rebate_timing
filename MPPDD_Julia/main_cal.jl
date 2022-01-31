# run line: JULIA_NUM_THREADS=32 julia main.jl


#request memory: probably just "add procs"
#nworkers = 48;  # num workers -- REQUEST WORKERS INSIDE CODE, NOT IN TERMINAL
#mem      = 5;  # request memory for each worker
#ENV["frbnyjuliamemory"]=string(mem)*"G"
#addprocs_frbny(nworkers)




## Packages

#using Pkg
using Distributed
#when running on the RAN cluster
using ClusterManagers

const nworkers = 40*4;
addprocs_slurm(nworkers,partition="long-40core",nodes=4, time="16:50:00",mem="150GB");
#when running on the RAN cluster
#addprocs_frbny(nworkers);

#addprocs(8);

@everywhere begin 
    using Distributed
    using LinearAlgebra
    #using SparseArrays, ParSpMatVec
    #using Arpack
    using SchumakerSpline, Interpolations
    using Distributions, Random
    using Optim
    using JLD
    using Plots
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
        V::Array{Float64,5}    # defined for (N_ages+1,N_t,N_a,N_ε, N_z)
        C::Array{Float64,5}    # defined for (N_ages,N_t,N_a,N_ε, N_z)
        S::Array{Float64,5}    # defined for (N_ages,N_t,N_a,N_ε, N_z)
        B::Array{Float64,5}    # defined for (N_ages,N_t,N_a,N_ε, N_z)
        Ap::Array{Float64,5}   # defined for (N_ages,N_t,N_a,N_ε, N_z)
        Q0::Array{Float64,3}   # defined for (N_ages,N_z,N_a) known states
        V_frontier::Array{Float64,6}   # defined for (N_ages,N_a,N_ε, N_z, sidx,N_a)
        B_frontier::Array{Float64,6}   # defined for (N_ages,N_a,N_ε, N_z, sidx, N_a)
        Aεz_dist::Array{Float64,4} #Equilibrium distribution across state (N_ages,a,z,ε)
    end

    function sol()
        V  = zeros(N_ages+1,N_t,N_a,N_ε, N_z); #value function
        C  = zeros(N_ages,N_t,N_a,N_ε, N_z);   #consumption
        S  = zeros(N_ages,N_t,N_a,N_ε, N_z);   # pay back rate
        B  = zeros(N_ages,N_t,N_a,N_ε, N_z);   # change in asset position
        Ap = zeros(N_ages,N_t,N_a,N_ε, N_z);   # asset position tomorrow
        Q0 = zeros(N_ages, N_z,N_a);           # interest rate function
        V_frontier= zeros(N_ages,N_a,N_ε, N_z, 3,N_a);   
        B_frontier= zeros(N_ages,N_a,N_ε, N_z, 3,N_a);   
        Aεz_dist = zeros(N_ages,N_a,N_ε,N_z);  # distribution over states


        return sol(V,C,S,B,Ap,Q0,V_frontier,B_frontier,Aεz_dist)
    end

    mutable struct model
        zprobs::Array{Float64,3}
        zgrid ::Array{Float64,2}
        εprobs::Array{Float64,2}
        εgrid ::Array{Float64,2}
        Ygrid::Array{Float64,3}
        asset_grid::Array{Float64,2}
        cbar::Float64
        exog_asset_demand::Float64
        γ::Int64
        ϕ_1::Float64 #average cost level
        ϕ_2::Float64 # curvature
        ϕ_3::Float64 # fixed costs
        reclaimrate::Float64 #h in the model... the fraction not taken by haircut
        r::Float64
        β::Float64
        λ::Float64
        transfer_grid::Array{Float64,1}

        model() = (modd = new();modd.asset_grid=zeros(N_ages,N_a); modd.cbar = 0.0; modd.exog_asset_demand = 0.0; modd.γ =1; return modd)
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
    end

    function hists()
        a0qtlhist = zeros(Nsim);
        zidxhist  = Array{Int64,2}(undef, Nsim, Thist*freq);
        εidxhist  = Array{Int64,2}(undef, Nsim, Thist*freq);
        ahist     = zeros(Nsim,Thist*freq);
        shist     = zeros(Nsim,Thist*freq);
        chist     = zeros(Nsim,Thist*freq);
        bhist     = zeros(Nsim,Thist*freq);
        qhist     = zeros(Nsim,Thist*freq);
        incomehist= zeros(Nsim,Thist*freq);
        agehist   = Array{Int64,2}(undef, Nsim, Thist*freq);
        mpchist   = Array{Float64,3}(undef,2,Nsim,Thist*freq);
        mpshist   = Array{Float64,3}(undef,2,Nsim,Thist*freq);
        mpbhist   = Array{Float64,3}(undef,2,Nsim,Thist*freq);
        mpqhist   = Array{Float64,3}(undef,2,Nsim,Thist*freq);

        return hists(chist,bhist,
            shist,ahist,qhist,εidxhist,
            zidxhist,agehist,incomehist,
            mpchist,mpshist,mpbhist,mpqhist,a0qtlhist)
    end

    mutable struct moments
        #moments from the simulation to return
        net_asset_supply::Float64
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
        net_asset_supply = 0.0;
        fr_neg_asset = 0.0
        avg_indiv_debt_income =0.0
        avg_debt =0.0
        avg_income =0.0
        avg_income_neg_asset =0.0
        extent_delinq =0.0
        avg_s =0.0
        fr_s1 =0.0
        mpc_cor_neg_asset = zeros(6,6)
        mpc_cor_pos_asset = zeros(3,3)
        mpc_ac0 = 0.0

        return moments(net_asset_supply,fr_neg_asset,avg_indiv_debt_income,avg_debt,avg_income,avg_income,extent_delinq,avg_s,fr_s1,mpc_cor_neg_asset,mpc_cor_pos_asset,mpc_ac0)
    end


    include("functions/incomeProcess.jl")
    include("functions/agent.jl")
    include("functions/equilibrium.jl")


    const qtol         = 1e-3; # tolerance for convergence in Q
    const maxiterq     = 30;   # how many iterations?
    const maxiterV     = 30;
    const Vtol         = 1e-3; #can be reasonably loose b/c iterating on Q too
    ## Parameters
    const freq         = 4; # 1=annual, 4=quarterly
    const γ            = 1;

    const debug_saves  = 0; # should we save stuff for debugging? If 1, then yes.
    print_lev = 0;

    ## State space

    const N_a          = 60;    # Number of asset states, 40-50
    const N_z          = 5;      # Number of permanent income states, 7-10
    const N_ε          = 4;         # Number of transitory income states, 5-21
    const N_t          = 3;      # Number of transfer levels

    const grid_curvature = 2.0;  #unevenly spaced grid points. >0 means more close to 0

    ## Income Process  - from Kaplan & Violante (2010)
    ρ          =  0.999; #should've been 1, but I worry about the ergod distribution.
    σ_η       = (0.01/freq)^(1/2);  # Variance of permanent income shocks = 0.01
    σ_z0      = (0.15/freq)^(1/2);  # Variance of permanent income shocks = 0.15
    σ_ε       = (0.05/freq)^(1/2);  # Variance of transitory income shocks = 0.05

    ## Simulation:
    const Nsim      = 10000;  #number of Agents
    const Tsim      = 12;     #number of periods in the simulation
    const Tburnin   = 32;     #number of burnin periods
    const Thist     = Tsim+Tburnin;

    const Nparams   = 6; #calibration parameters
    const Nmoments  = 8; # calibration moments

    ## Ages:

    const YrsT = 85;
    const Yrs_retire = 65; #corresponds to 65
    const Yrs0 = 25;
    const N_ages = (YrsT-Yrs0+1)*freq;
    const age_retire = (Yrs_retire-Yrs0+1) * freq;
    const ageprobs = ones(N_ages+1); #vcat(ones(N_ages-1)*(N_ages-1.0)/(T_retire - T0)/freq, 1.0/(T-T_retire)/freq );
    const agegrid = collect(LinRange(1, N_ages, N_ages));
    ## Agents begin at  age 25, die at 95, retire at 65

    # Survival Probabilities
    const survivalprobs = ones(N_ages); # N_ages set to 0 in agent.jl


	const fixedassetgrid = false;
	const unevengrid = true;
	const natborrowlimit = true;

    constQ = false; #set q = 1/(1+r) or every debt level
    fullcommit = false; # do not allow s<1, which also implies constQ



    print("Done defining constants.")
    """
    Takes in parameters, outputs moments / parameters to CSV
    """
    function params_moments(iofilename, parvals, j)
        parval = parvals[j];
        
        ht = hists();
        
        mod      = model();
        mod.β    = parval[:β];
        mod.r    = parval[:r];
        mod.ϕ_1  = parval[:ϕ_1];
        mod.ϕ_2  = parval[:ϕ_2];
        mod.ϕ_3  = parval[:ϕ_3];
        mod.λ  = parval[:λ];
        mod.reclaimrate = parval[:reclaimrate];

        mmt = moments();
        age_peak       = (Yrs_retire-10)*freq; # corresponds to 50-55, at which income peaks
        mod.Ygrid, mod.zgrid, mod.εgrid = get_Ygrid(age_retire, age_peak, N_ages, N_z, N_ε, σ_ε, σ_z0, σ_η, ρ); # level grids of income

        mod.zprobs    = get_zprobs(N_ages, N_z, σ_z0, σ_η, ρ); # transition matrices for permanent income
        mod.εprobs    = get_εprobs(N_ages, N_ε, σ_ε);          # transition matrices for transitory income

        ms = sol();
        pidnum = getpid();
        workernum = myid() - 1 ;

        println("Drawing shocks on PID $pidnum, $workernum at j $j !");
        draw_shocks!(mod,ht, 12281951);
        
        ms.Q0 = solveQ!(mod,ms,ht,mmt) ;
        if fullcommit == true
            saveloc = string(saveroot,"modEnvr_commitment",j,".jld")
            @save saveloc mod
            saveloc = string(saveroot,"solMats_commitment",j,".jld")
            @save saveloc ms
        elseif constQ == true
            saveloc = string(saveroot,"modEnvr_constQ",j,".jld")
            @save saveloc mod
            saveloc = string(saveroot,"solMats_constQ",j,".jld")
            @save saveloc ms
        else  
            saveloc = string(saveroot,"modEnvr_p",j,".jld")
            @save saveloc mod
            saveloc = string(saveroot,"solMats_p",j,".jld")
            @save saveloc ms
        end 
        
        println("About to sim_hists")
        sim_hists!(mod, ht, ms, mmt, N_ages, N_a,N_z, N_ε);
        println("Done sim_hists")
        
        saveloc = string(saveroot,"mmts_p",j,".jld")
        @save saveloc mmt
        println("Done saving..")
        outvec = reshape([j, mod.r, mod.β, mod.ϕ_1, mod.ϕ_2, mod.ϕ_3, mod.λ , mmt.fr_neg_asset, mmt.avg_debt,mmt.avg_s, mmt.extent_delinq,mmt.fr_s1, mmt.avg_income,mmt.avg_income_neg_asset,mmt.mpc_ac0],1,Nparams+Nmoments+1 );

        iohr = open(iofilename,"a")
        writedlm(iohr, round.(outvec,digits=5), ',')
        close(iohr);

        #sharedmom[j,1] = mmt.avg_indiv_debt_income;
        #sharedmom[j,2] = mmt.avg_s;
        #sharedmom[j,3] = mmt.fr_neg_asset;
        #sharedmom[j,4] = mmt.fr_s1;

        
    end
    
    
end # end parallel everywhere block


function cal_fullcommit(param_grids::OrderedDict{Symbol,Array{Float64,1}}, file::String = string(saveroot,"out.csv"))

    file =string(saveroot,"fullcommit_params.csv");

    params_header = reshape(["idx_num","r","beta","phi_1", "phi_2", "phi_3", "borrow_wdg","fr_neg_asset", "avg_debt", "avg_s","extent_delinq", "fr_s1", "avg_income","avg_income_neg_asset","MPC_a0"],1,Nparams+Nmoments+1 );

    # note: this will 'delete/overwrite' the file if it already exists
    writedlm(file,params_header,',') # note: last argument is the delimiter which should be the same as below

    #first do it with full commitment
    fullcommit_old = fullcommit;
    global fullcommit = true;
    

    β = 0.9;ω = 0.10; ϕ_1 = 0.1; ϕ_2 = 2.0; ϕ_3 = 0.0; λ = 0.0;
    recrate =1.00;
    println("set stuff");
    parvals = OrderedDict{Int64, OrderedDict{Symbol,Float64}}()
    i =1;
    for β in values(param_grids[:β])
        
        β_f = β^(1/freq); # correct for frequency of model
        temp = OrderedDict{Symbol, Float64}(
                    :β => β_f,
                    :ϕ_1 => ϕ_1/freq,
                    :ϕ_2 => ϕ_2,
                    :ϕ_3 => ϕ_3,
                    :r =>   (1/β_f-1)*(1-ω),
                    :λ =>   λ,
                    :reclaimrate => recrate);
        println("created temp struct, $(temp[:r])")            
        push!(parvals, i => temp)
        i+=1
    end

    # better for heavy lifting, allow julia to distribute
    pmap(j->params_moments(file, parvals, j), 1:length(parvals))

    global fullcommit = fullcommit_old;
end


function cal_constQ(param_grids::OrderedDict{Symbol,Array{Float64,1}}, file::String = string(saveroot,"out.csv"))

    file =string(saveroot,"constQ_params.csv");

    params_header = reshape(["idx_num","r","beta","phi_1", "phi_2", "phi_3", "borrow_wdg","fr_neg_asset", "avg_debt", "avg_s","extent_delinq", "fr_s1", "avg_income","avg_income_neg_asset","MPC_a0"],1,Nparams+Nmoments+1 );

    # note: this will 'delete/overwrite' the file if it already exists
    writedlm(file,params_header,',') # note: last argument is the delimiter which should be the same as below

    #first do it with full commitment
    constQ_old = constQ;
    global constQ = true;
    

    β = 0.9;ω = 0.10; ϕ_1 = 0.1; ϕ_2 = 2.0; ϕ_3 = 0.0; λ = 0.0;
    recrate =1.00;
    println("set stuff");
    parvals = OrderedDict{Int64, OrderedDict{Symbol,Float64}}()
    i =1;
    for β in values(param_grids[:β]), ϕ_1 in values(param_grids[:ϕ_1]), ϕ_2 in values(param_grids[:ϕ_2])
        
        β_f = β^(1/freq); # correct for frequency of model
        temp = OrderedDict{Symbol, Float64}(
                    :β => β_f,
                    :ϕ_1 => ϕ_1/freq,
                    :ϕ_2 => ϕ_2,
                    :ϕ_3 => ϕ_3,
                    :r =>   (1/β_f-1)*(1-ω),
                    :λ =>   λ,
                    :reclaimrate => recrate);
        push!(parvals, i => temp)
        i+=1
    end

    # better for heavy lifting, allow julia to distribute
    pmap(j->params_moments(file, parvals, j), 1:length(parvals))
    global constQ = constQ_old;
end

# parameter grids, annual values (where relevant)
#param_grids  = OrderedDict{Symbol,Array{Float64,1}}([
#        (:ω, collect(LinRange(0.00, 0.50, 3)) ),
#        (:β, [0.95, 0.99, 0.995]),
#        (:ϕ_1,[0.01, 0.1, 0.3]),
#        (:ϕ_2, [1.25, 2.0, 4, 5.5]),
#        (:ϕ_3, [0.0, 0.1, 0.2]),
#        (:λ , collect(LinRange(0.0, 0.06, 3))) ]);

#estimate_elasticities(param_grids)
beta_grid = OrderedDict{Symbol,Array{Float64,1}}([(:β, collect(LinRange(0.83,0.97,30)) ) ]) ;
cal_fullcommit( beta_grid )

beta_phi_grid = OrderedDict{Symbol,Array{Float64,1}}([
    (:β,collect(LinRange(0.83,0.97,15))  ), (:ϕ_1, collect(LinRange(0.05,0.25,5)) ), (:ϕ_2, collect(LinRange(1.5,4,6)) ) ]) ;
cal_constQ( beta_phi_grid )


rmprocs(workers())

#=+++++++++++++++++++++++++++++++
#run it once with some pre-set parameters
β_f = 0.99^(1/freq); # correct for frequency of model
ϕ_1 = 0.5;
ϕ_2 = 4;
ω = .5;
parvals = OrderedDict{Int64, OrderedDict{Symbol,Float64}}()
parhr = OrderedDict{Symbol, Float64}(
            :β => β_f,
            :ϕ_1 => ϕ_1/freq,
            :ϕ_2 => ϕ_2,
            :r =>   (1/β_f-1)*(1-ω))
push!(parvals, 0=> parhr);
params_moments("run1.csv", parvals, 0)

#++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++
## diagnostic plots:
#++++++++++++++++++++++++++++

for j = 1:length(Nparams+1)
    
    j=607
    saveloc = string(saveroot,"/matlab_conversion/Aug4/","eqmQ_p",j,".jld")
    @load saveloc Q0
    saveloc = string(saveroot,"/matlab_conversion/Aug4/","solMats_p",j,".jld")
    @load saveloc ms
    saveloc = string(saveroot,"/matlab_conversion/Aug4/","modEnvr_p",j,".jld")
    @load saveloc mod
    saveloc = string(saveroot,"SavedOutput/param_search_May13/","simHists_p",j,".jld")
    @load saveloc ht

    pltages = [1 10 (age_retire-1)  (age_retire+1)]
    asset_grid_ages = hcat(mod.asset_grid[pltages[1],:],mod.asset_grid[pltages[2],:],mod.asset_grid[pltages[3],:],mod.asset_grid[pltages[4],:]);

    for zi=1:N_z

        
        
        Qages = hcat( Q0[pltages[1],zi,:],Q0[pltages[2],zi,:],Q0[pltages[3],zi,:], Q0[pltages[4],zi,:] )
        plot(asset_grid_ages,Qages, label=["Age 1" "Age 2" "Age 3" "Age 4"],legend=:bottomright, lw=3)
        plot!(title = "Debt price schedule", ylabel="Equilibrium debt price", xlabel = "a'", size=(1500,600))
        savefig(string("Q_ages_phi",ϕ_2,"_zi",zi,".png"))

        #=N_a_neg = sum(mod.asset_grid[2,:].<0)
        Qassets = zeros(age_retire-1, N_a_neg)
        for ai=1:N_a_neg
            for t=2:age_retire
                natborrowing = minimum(mod.asset_grid[t,:]);
                if mod.asset_grid[2,ai] > natborrowing
                    Qassets[t-1,ai] = LinearInterpolation(mod.asset_grid[t,:], Q0[t-1,zi,:])(mod.asset_grid[2,ai]);
                else
                    Qassets[t-1,ai] = NaN
                end
            end
        end
        labs = round.(mod.asset_grid[2,1:N_a_neg]*10.0)/10.0
        plot!(1:(age_retire-1),Qassets,title = "Price as a function of age", ylabel="Equilibrium debt price", xlabel = "Age",
            label= labs', legend=:bottomright, lw=3, size=(1500,600)) #label=["Age 1" "Age 15" "Age 25" "Age 35"],
        savefig(string("Q_assets_phi",ϕ_2,"_zi",zi,".png"))
        =#
        for ei=1:N_ε

            Vages = hcat( ms.V[pltages[1],1,:,ei,zi],ms.V[pltages[2],1,:,ei,zi],ms.V[pltages[3],1,:,ei,zi], ms.V[pltages[4],1,:,ei,zi] );
            plot(asset_grid_ages,Vages,label=["Age 1" "Age 2" "Age 3" "Age 4"],legend=:bottomright ,lw=3)
            plot!(title = "Value Function", ylabel="", xlabel = "a", size=(1500,600))
            savefig(string("V_ages_phi",ϕ_2,"_zi",zi,"_ei",ei,".png"))

            Cages = hcat( ms.C[pltages[1],1,:,ei,zi],ms.C[pltages[2],1,:,ei,zi],ms.C[pltages[3],1,:,ei,zi], ms.C[pltages[4],1,:,ei,zi] );
            plot(asset_grid_ages,Cages,label=["Age 1" "Age 2" "Age 3" "Age 4"],legend=:bottomright, lw=3)
            plot!(title = "Consumption", ylabel="", xlabel = "a", size=(1500,600))

            savefig(string("C_ages_phi",ϕ_2,"_zi",zi,"_ei",ei,".png"))

            Sages = hcat( ms.S[pltages[1],1,:,ei,zi],ms.S[pltages[2],1,:,ei,zi],ms.S[pltages[3],1,:,ei,zi], ms.S[pltages[4],1,:,ei,zi] );
            plot(asset_grid_ages,Sages,label=["Age 1" "Age 2" "Age 3" "Age 4"],legend=:bottomright, lw=3)
            plot!(title = "Delinquency policy", ylabel="", xlabel = "a", size=(1500,600))
            savefig(string("S_ages_phi",ϕ_2,"_zi",zi,"_ei",ei,".png"))

            Apages = hcat( ms.Ap[pltages[1],1,:,ei,zi]-mod.asset_grid[pltages[1],:],ms.Ap[pltages[2],1,:,ei,zi]-mod.asset_grid[pltages[2],:],ms.Ap[pltages[3],1,:,ei,zi]-mod.asset_grid[pltages[3],:], ms.Ap[pltages[4],1,:,ei,zi]-mod.asset_grid[pltages[4],:] );
            plot(asset_grid_ages,Apages,label=["Age 1" "Age 2" "Age 3" "Age 4"],legend=:bottomright, lw=3)
            plot!(title = "A' policy", ylabel="a'-a", xlabel = "a", size=(1500,600))
            savefig(string("Ap_ages_phi",ϕ_2,"_zi",zi,"_ei",ei,".png"))

        end



        # experiment giving transitory income shock:
        Apages = hcat((ms.Ap[pltages[1],2, :,3,zi]- ms.Ap[pltages[1],1, :,3,zi])./ (mod.transfer_grid[2] - mod.transfer_grid[1]),
                      (ms.Ap[pltages[2],2,:,3,zi]- ms.Ap[pltages[2],1,:,3,zi])./ (mod.transfer_grid[2] - mod.transfer_grid[1]),
                      (ms.Ap[pltages[3],2,:,3,zi]- ms.Ap[pltages[3],1,:,3,zi])./ (mod.transfer_grid[2] - mod.transfer_grid[1]),
                      (ms.Ap[pltages[4],2,:,3,zi]- ms.Ap[pltages[4],1,:,3,zi])./ (mod.transfer_grid[2] - mod.transfer_grid[1]) )
        epschng = round((mod.transfer_grid[2] - mod.transfer_grid[1])*100.0)/100.0
        plot(asset_grid_ages,Apages)
        plot!(title = "A' change with rebate shock of $epschng", ylabel="a' change / rebate ", xlabel = "a",
            label=["Age 1" "Age 2" "Age 3" "Age 4"],legend=:bottomright, lw=3, size=(1500,600))
        savefig(string("ApDif_ages_phi",ϕ_2,"_zi",zi,".png"))

        delQages = zeros(N_a,4 )
        for ai=1:N_a
            delQages[ai,1] = (LinearInterpolation(mod.asset_grid[pltages[1]+1,:], Q0[pltages[1],zi,:])(ms.Ap[pltages[1],2, ai,3,zi]) - LinearInterpolation(mod.asset_grid[pltages[1]+1,:], Q0[pltages[1],zi,:])(ms.Ap[pltages[1],1, ai,3,zi]))/
                (ms.Ap[pltages[1],2, ai,3,zi]- ms.Ap[pltages[1],1, ai,3,zi]);
            delQages[ai,2] = (LinearInterpolation(mod.asset_grid[pltages[2]+1,:], Q0[pltages[2],zi,:])(ms.Ap[pltages[2],2, ai,3,zi]) - LinearInterpolation(mod.asset_grid[pltages[2]+1,:], Q0[pltages[2],zi,:])(ms.Ap[pltages[2],1, ai,3,zi]))/
                (ms.Ap[pltages[2],2, ai,3,zi]- ms.Ap[pltages[2],1, ai,3,zi]);
            delQages[ai,3] = (LinearInterpolation(mod.asset_grid[pltages[3]+1,:], Q0[pltages[3],zi,:])(ms.Ap[pltages[3],2, ai,3,zi]) - LinearInterpolation(mod.asset_grid[pltages[3]+1,:], Q0[pltages[3],zi,:])(ms.Ap[pltages[3],1, ai,3,zi]))/
                (ms.Ap[pltages[3],2, ai,3,zi]- ms.Ap[pltages[3],1, ai,3,zi]);
            delQages[ai,4] = (LinearInterpolation(mod.asset_grid[pltages[4]+1,:], Q0[pltages[4],zi,:])(ms.Ap[pltages[4],2, ai,3,zi]) - LinearInterpolation(mod.asset_grid[pltages[4]+1,:], Q0[pltages[4],zi,:])(ms.Ap[pltages[4],1, ai,3,zi]))/
                (ms.Ap[pltages[4],2, ai,3,zi]- ms.Ap[pltages[4],1, ai,3,zi]);
        end
        epschng = round((mod.transfer_grid[2] - mod.transfer_grid[1])*100.0)/100.0
        plot(asset_grid_ages,delQages)
        plot!(title = "dQdA change with rebate of $epschng", ylabel="q change / rebate ", xlabel = "a", label=["Age 1" "Age 2" "Age 3" "Age 4"],
            legend=:bottomright, lw=3, size=(1500,600))
        savefig(string("QDif_phi",ϕ_2,"_zi",zi,".png"))



        # experiment giving transitory income shock:
        Cages = hcat((ms.C[pltages[1],2, :,3,zi]- ms.C[pltages[1],1, :,3,zi])./ (mod.transfer_grid[2] - mod.transfer_grid[1]),
                      (ms.C[pltages[2],2,:,3,zi]- ms.C[pltages[2],1,:,3,zi])./ (mod.transfer_grid[2] - mod.transfer_grid[1]),
                      (ms.C[pltages[3],2,:,3,zi]- ms.C[pltages[3],1,:,3,zi])./ (mod.transfer_grid[2] - mod.transfer_grid[1]),
                      (ms.C[pltages[4],2,:,3,zi]- ms.C[pltages[4],1,:,3,zi])./ (mod.transfer_grid[2] - mod.transfer_grid[1]) )
        epschng = round((mod.transfer_grid[2] - mod.transfer_grid[1])*100.0)/100.0

        plot(asset_grid_ages,Cages)
        plot!(title = "C change with rebate of $epschng", ylabel="c change / rebate ", xlabel = "a",
            label=["Age 1" "Age 2" "Age 3" "Age 4"],legend=:bottomright, lw=3, size=(1500,600))
        savefig(string("CDif_ages_phi",ϕ_2,"_zi",zi,".png"))


        Bages = hcat((ms.B[pltages[1],2, :,3,zi]- ms.B[pltages[1],1, :,3,zi])./ (mod.transfer_grid[2] - mod.transfer_grid[1]),
                      (ms.B[pltages[2],2,:,3,zi]- ms.B[pltages[2],1,:,3,zi])./ (mod.transfer_grid[2] - mod.transfer_grid[1]),
                      (ms.B[pltages[3],2,:,3,zi]- ms.B[pltages[3],1,:,3,zi])./ (mod.transfer_grid[2] - mod.transfer_grid[1]),
                      (ms.B[pltages[4],2,:,3,zi]- ms.B[pltages[4],1,:,3,zi])./ (mod.transfer_grid[2] - mod.transfer_grid[1]) )
        epschng = round((mod.transfer_grid[2] - mod.transfer_grid[1])*100.0)/100.0

        plot(asset_grid_ages,Bages)
        plot!(title = "B change with rebate of $epschng", ylabel="b change / rebate ", xlabel = "a",
            label=["Age 1" "Age 2" "Age 3" "Age 4"],legend=:bottomright, lw=3, size=(1500,600))
        savefig(string("BDif_ages_phi",ϕ_2,"_zi",zi,".png"))


        Sages = hcat((ms.S[pltages[1],2, :,3,zi]- ms.S[pltages[1],1, :,3,zi])./ (mod.transfer_grid[2] - mod.transfer_grid[1]),
                      (ms.S[pltages[2],2,:,3,zi]- ms.S[pltages[2],1,:,3,zi])./ (mod.transfer_grid[2] - mod.transfer_grid[1]),
                      (ms.S[pltages[3],2,:,3,zi]- ms.S[pltages[3],1,:,3,zi])./ (mod.transfer_grid[2] - mod.transfer_grid[1]),
                      (ms.S[pltages[4],2,:,3,zi]- ms.S[pltages[4],1,:,3,zi])./ (mod.transfer_grid[2] - mod.transfer_grid[1]) )
        epschng = round((mod.transfer_grid[2] - mod.transfer_grid[1])*100.0)/100.0
        plot(asset_grid_ages,Sages)
        plot!(title = "S change with rebate of $epschng", ylabel="s change / rebate ", xlabel = "a",
            label=["Age 1" "Age 2" "Age 3" "Age 4"],legend=:bottomright, lw=3, size=(1500,600))
        savefig(string("SDif_ages_phi",ϕ_2,"_zi",zi,".png"))
        
    end
end
 #
