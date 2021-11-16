# run line: JULIA_NUM_THREADS=32 julia main.jl


#request memory: probably just "add procs"
#nworkers = 5;  # num workers -- REQUEST WORKERS INSIDE CODE, NOT IN TERMINAL
#mem      = 1;  # request memory for each worker
#ENV["frbnyjuliamemory"]=string(mem)*"G"
#addprocs_frbny(nworkers)




## Packages

#using Pkg
using Distributed
#when running on the RAN cluster
#using ClusterManagers

const nworkers = 4;
addprocs(nworkers)
#addprocs_slurm(nworkers,partition="long-40core",nodes=4, time="16:50:00",mem="150GB");
#when running on the RAN cluster
#addprocs_frbny(nworkers);


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


    ## State space

    const N_a          = 40;    # Number of asset states, 40-50
    const N_z          = 5;      # Number of permanent income states, 7-10
    const N_ε          = 4;         # Number of transitory income states, 5-21
    const N_t          = 3;      # Number of transfer levels

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

    """
    Takes in parameters, outputs moments / parameters to CSV 
    """
    function params_moments(j)
        
        
        ht = hists();

        mod      = model();
        
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
        #ms.Q0 = solveQ!(mod,ms,ht,mmt) 
        
        saveloc = string(saveroot,"modEnvr_p",j,".jld")
        @load saveloc mod
        saveloc = string(saveroot,"solMats_p",j,".jld")
        @load saveloc ms
        
        sim_hists!(mod, ht, ms, mmt, N_ages, N_a,N_z, N_ε);
        
        saveloc = string(saveroot,"simHists_p",j,".jld")
        @save saveloc ht
        #saveloc = string(saveroot,"eqmQ_p",j,".jld")
        #@save saveloc ms.Q0
        saveloc = string(saveroot,"mmts_p",j,".jld")
        @save saveloc mmt
            
    end

    
end # end parallel everywhere block




function estimate_elasticities()


    # better for heavy lifting, allow julia to distribute
    pmap(j->params_moments(j), 409:411)

    # Save a final output (note: jld2 seems to work better for shared arrays)
    # with additional objects
    #saveloc = string(saveroot,"MPPDD_test.jld2")
    #@save saveloc momentsout
end

estimate_elasticities()

rmprocs(workers())



#++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++
## diagnostic plots:
#=++++++++++++++++++++++++++++

for j = 1:length(Nparams+1)
    
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
