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
    println(saveroot)
    println("Beginning structs");
    mutable struct sol
        V::Array{Float64,5}    # defined for (N_ages+1,N_t,N_a,N_ε, N_z)
        C::Array{Float64,5}    # defined for (N_ages,N_t,N_a,N_ε, N_z)
        S::Array{Float64,5}    # defined for (N_ages,N_t,N_a,N_ε, N_z)
        B::Array{Float64,5}    # defined for (N_ages,N_t,N_a,N_ε, N_z)
        Ap::Array{Float64,5}   # defined for (N_ages,N_t,N_a,N_ε, N_z)
        Q0::Array{Float64,3}   # defined for (N_ages,N_z,N_a) known states

        Aεz_dist::Array{Float64,4} #Equilibrium distribution across state (N_ages,a,z,ε)
    end

    function sol()
        V  = zeros(N_ages+1,N_t,N_a,N_ε, N_z);
        C  = zeros(N_ages,N_t,N_a,N_ε, N_z);
        S  = zeros(N_ages,N_t,N_a,N_ε, N_z);
        B  = zeros(N_ages,N_t,N_a,N_ε, N_z);
        Ap = zeros(N_ages,N_t,N_a,N_ε, N_z);
        Q0 = zeros(N_ages, N_z,N_a);
        Aεz_dist = zeros(N_ages,N_a,N_ε,N_z);
        return sol(V,C,S,B,Ap,Q0,Aεz_dist)
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
        ϕ_1::Float64 #average cost level
        ϕ_2::Float64 # curvature
        ϕ_3::Float64 # fixed costs
        reclaimrate::Float64
        r::Float64
        β::Float64
        λ::Float64
        transfer_grid::Array{Float64,1}

        model() = (modd = new();modd.asset_grid=zeros(N_ages,N_a); modd.cbar = 0.0; modd.γ =1; return modd)
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
        mpc_cor_neg_asset = zeros(6,6)
        mpc_cor_pos_asset = zeros(3,3)
        mpc_ac0 = 0.0

        return moments(fr_neg_asset,avg_indiv_debt_income,avg_debt,avg_income,avg_income,extent_delinq,avg_s,fr_s1,mpc_cor_neg_asset,mpc_cor_pos_asset,mpc_ac0)
    end

    println("Done defining structs.")
    include("functions/incomeProcess.jl")
    include("functions/agent.jl")
    include("functions/equilibrium.jl")
    println("Done including other functions.")

    const qtol         = 1e-3; # tolerance for convergence in Q
    const maxiterq     = 30;   # how many iterations?
    const maxiterV     = 30;
    const Vtol         = 1e-3; #can be reasonably loose b/c iterating on Q too
    ## Parameters
    const freq         = 4; # 1=annual, 4=quarterly
    const γ            = 1;


    ## State space

    const N_a          = 60;    # Number of asset states, 40-50
    const N_z          = 5;      # Number of permanent income states, 7-10
    const N_ε          = 4;         # Number of transitory income states, 5-21
    const N_t          = 3;      # Number of transfer levels

    const grid_curvature = 2.0;  #unevenly spaced grid points. >0 means more close to 0

    ## Income Process  - from Kaplan & Violante (2010)
    ρ          =  1.0;
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

    ## Agents begin at  age 25, die at 95, retire at 65

    # Survival Probabilities
    const survivalprobs = ones(N_ages); # N_ages set to 0 in agent.jl


	const fixedassetgrid = false;
	const unevengrid = true;
	const natborrowlimit = true;

    fullcommit = false;


    print("Done defining constants.")
    """
    Takes in parameters, outputs moments / parameters to CSV
    """
    function params_moments(iofilename, parval, j)
        ht = hists();

        mod      = model();
        mod.r    = parval[:r];
        mod.β    = parval[:β];
        mod.ϕ_1  = parval[:ϕ_1];
        mod.ϕ_2  = parval[:ϕ_2];
        mod.ϕ_3  = parval[:ϕ_3];
        mod.λ  = parval[:λ];

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

        println("Solving on PID $pidnum, $workernum at j $j !");
        # optionally, test the timing of this to see how things output
        #sleep(ceil(Int,rand()*5))
        ms.Q0 = solveQ!(mod,ms,ht,mmt) 
        
        saveloc = string(saveroot,"modEnvr_p",j,".jld")
        @save saveloc mod
        saveloc = string(saveroot,"solMats_p",j,".jld")
        @save saveloc ms
        println("Simulating!  $pidnum, $workernum at j $j");        
        sim_hists!(mod, ht, ms, mmt, N_ages, N_a,N_z, N_ε);
        
        #saveloc = string(saveroot,"simHists_p",j,".jld")
        #@save saveloc ht
        #saveloc = string(saveroot,"eqmQ_p",j,".jld")
        #@save saveloc ms.Q0
        saveloc = string(saveroot,"mmts_p",j,".jld")
        @save saveloc mmt
        
        outvec = reshape([j, mod.r, mod.β, mod.ϕ_1, mod.ϕ_2, mod.ϕ_3, mod.λ , mmt.fr_neg_asset, mmt.avg_debt,mmt.avg_s, mmt.extent_delinq,mmt.fr_s1, mmt.avg_income,mmt.avg_income_neg_asset,mmt.mpc_ac0],1,Nparams+Nmoments+1 );
        
        iohr = open(iofilename,"a")
        writedlm(iohr, round.(outvec,digits=5), ',')
        close(iohr);

        #sharedmom[j,1] = mmt.avg_indiv_debt_income;
        #sharedmom[j,2] = mmt.avg_s;
        #sharedmom[j,3] = mmt.fr_neg_asset;
        #sharedmom[j,4] = mmt.fr_s1;

    #    iohr2 = open(string("corMat_negasset_p", j,".csv"),"a");
    #    writedlm(iohr2, round.(mmt.mpc_cor_neg_asset,digits=5), ',');
    #    close(iohr2);
        
        
    #    iohr3 = open(string("corMat_posasset_p", j,".csv"),"a");
    #    writedlm(iohr3, round.(mmt.mpc_cor_pos_asset,digits=5), ',');
    #    close(iohr3);
    

   #     outvec2 = [j reshape(mmt.mpc_cor_neg_asset,1,25 )];
   #     outvec3 = [j reshape(mmt.mpc_cor_pos_asset ,1,9 )];
    
   #     file2 = string(saveroot,"corMat_negasset.csv")
   #     file3 = string(saveroot,"corMat_posasset.csv")

   #     iohr2 = open(file2,"a")
   #     writedlm(iohr2, round.(outvec2,digits=5), ',')
   #     close(iohr2);
        
   #     iohr3 = open(file3,"a")
   #     writedlm(iohr3, round.(outvec3,digits=5), ',')
   #     close(iohr3);
    
        
    end
    
    
end # end parallel everywhere block


function estimate_elasticities(param_grids::OrderedDict{Symbol,Array{Float64,1}}, file::String = string(saveroot,"out.csv"))

    file =string(saveroot,"explore_params.csv");

    params_header = reshape(["idx_num","r","beta","phi_1", "phi_2", "phi_3", "borrow_wdg","fr_neg_asset", "avg_debt", "avg_s","extent_delinq", "fr_s1", "avg_income","avg_income_neg_asset","MPC_a0"],1,Nparams+Nmoments+1 );
    
    # note: this will 'delete/overwrite' the file if it already exists
    writedlm(file,params_header,',') # note: last argument is the delimiter which should be the same as below

    #file2 = string(saveroot,"corMat_negasset.csv")
    #file3 = string(saveroot,"corMat_posasset.csv")
    #writedlm(file2,reshape(["idx_num",
    #    "cor(-a,-a)" ,"cor(-a, mpc)","cor(-a,mpb)" ,"cor(-a,mps)" ,"cor(-a,mpq)" ,
    #    "cor(mpc,-a)","cor(mpc,mpc)","cor(mpc,mpb)","cor(mpc,mps)","cor(mpc,mpq)",
    #    "cor(mpb,-a)","cor(mpb,mpc)","cor(mpb,mpb)","cor(mpb,mps)","cor(mpb,mpq)", 
    #    "cor(mps,-a)","cor(mps,mpc)","cor(mps,mpb)","cor(mps,mps)","cor(mps,mpq)",
    #    "cor(mpq,-a)","cor(mpq,mpc)","cor(mpq,mpb)","cor(mpq,mps)","cor(mpq,mpq)"],1,26), ',') #
    #writedlm(file3,reshape(["idx_num", 
    #    "cor(a,a)"  ,"cor(a, mpc)"  ,"cor(a,mpb)"   ,
    #    "cor(mpc,a)","cor(mpc, mpc)","cor(mpc, mpb)",
    #    "cor(mpb,a)","cor(mpb, mpc)","cor(mpb, mpb)"],1,10),',') # note: last argument is the delimiter which should be the same as below


    parvals = OrderedDict{Int64, OrderedDict{Symbol,Float64}}()
    i =1;
    for ω in values(param_grids[:ω]), β in values(param_grids[:β]), ϕ_1 in values(param_grids[:ϕ_1]), ϕ_2 in values(param_grids[:ϕ_2]), ϕ_3 in values(param_grids[:ϕ_3]), λ in values(param_grids[:λ])
        
        β_f = β^(1/freq); # correct for frequency of model
        temp = OrderedDict{Symbol, Float64}(
                    :β => β_f,
                    :ϕ_1 => ϕ_1/freq,
                    :ϕ_2 => ϕ_2,
                    :ϕ_3 => ϕ_3,
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
        (:ω, collect(LinRange(0.00, 0.50, 3)) ),
        (:β, [0.925, 0.95, 0.975]),
        (:ϕ_1,[0.1, 0.3, 0.5]),
        (:ϕ_2, [1.25, 2.0, 4, 5.5]),
        (:ϕ_3, [0.0, 0.1, 0.2]),
        (:λ , collect(LinRange(0.0, 0.06, 3))) ]);

estimate_elasticities(param_grids)

rmprocs(workers())
