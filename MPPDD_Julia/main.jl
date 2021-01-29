# run line: JULIA_NUM_THREADS=32 julia main.jl


## Packages

#using Pkg
#include("import.jl") # uncomment to add required packages

using LinearAlgebra
using SparseArrays, ParSpMatVec
using Arpack
using SchumakerSpline, Interpolations
using Distributions, Random
using Optim
using JLD
using Plots

import Base.Threads.@spawn
#cd("/Users/meghanagaur/DropBox/Consumption_Model_Codes/julia/")     # personal
#cd("C:\\Users\\rcemxg23\\DropBox\\Consumption_Model_Codes\\julia") # ran virtual desktop
cd("/home/david/workspace/rebate_timing/MPPDD_Julia/")
#cd("/Users/rcedxm19/Dropbox/MPPDD/codes/MPPDD_Julia")

saveroot = "/home/david/workspace/rebate_timing/MPPDD_Julia/"

mutable struct sol
    V::Array{Float64,4}    # defined for (T+1,N_k,N_ε, N_z)
    C::Array{Float64,4}    # defined for (T,N_k,N_ε, N_z)
    S::Array{Float64,4}    # defined for (T,N_k,N_ε, N_z)
    B::Array{Float64,4}    # defined for (T,N_k,N_ε, N_z)
    Ap::Array{Float64,4}   # defined for (T,N_k,N_ε, N_z)
    Q0::Array{Float64,3}   # 
    Aεz_dist::Array{Float64,4} #Equilibrium distribution across state (T,a,z,ε)
    sol() = new()
end


mutable struct model
    #fill this in with model objects....
    zprobs::Array{Float64,3}
    zgrid ::Array{Float64,2}
    εprobs::Array{Float64,2}
    εgrid ::Array{Float64,2}
    Ygrid::Array{Float64,3}
    asset_grid::Array{Float64,2}
    cbar::Float64
    β::Float64
    γ::Int64
    ϕ_1::Float64
    ϕ_2::Float64
        
    rebate:: Float64 # added to BC at t_reb (argument in backwardSolve)
    t_rebate::Int64
    model() = new()
end

mutable struct hists
    #histories for:
    chist::Array{Float64,2}
    bhist::Array{Float64,2}
    shist::Array{Float64,2}
    ahist::Array{Float64,2}
    #should be drawn once:
    εidxhist::Array{Int64,2}
    zidxhist::Array{Int64,2}
    agehist::Array{Int64,2}
    
    a0qtlhist::Array{Float64,1} #quantile of a in the asset distribution for first period.

    hists() = new()
end

include("functions/incomeProcess.jl")
include("functions/agent.jl")

const qtol         = 1e-5; # tolerance for convergence in Q
const maxiterq     = 40;   # how many iterations?

## Parameters
const freq         = 1; # 1=annual, 4=quarterly
const r            = (1 + 0.05)^(1/freq) - 1;
const β            = 0.925^(1/freq);
const γ            = 1;


## State space

const N_k          = 200;      # Number of asset states, 40-50
const N_z          = 7;      # Number of permanent income states, 7-10
const N_ε          = 5;	    # Number of transitory income states, 5-21
const N_s          = 15;      # Number of delinquency grid points, 10-5

## Income Process  - from Kaplan & Violante (2010)
ρ	      =  1.0;
σ_η       = (0.01/freq)^(1/2);  # Variance of permanent income shocks = 0.01
σ_z0      = (0.15/freq)^(1/2);  # Variance of permanent income shocks = 0.15
σ_ε       = (0.05/freq)^(1/2);  # Variance of transitory income shocks = 0.05

## Simulation:
const Nsim      = 500000; #number of Agents
const Tsim      = 10;    #number of periods in the simulation
const Tburnin   = 20;    #number of burnin periods

## Ages:
const T = 71;        #corresponds to 95
const T_retire = 41; #corresponds to 65
const t0 = 25; 
## Agents begin at  age 25, die at 95, retire at 65

# Survival Probabilities
const survivalprobs = ones(T); # T set to 0 in agent.jl


const ϕ_2grid = [1.5];


## Solve model
function solveQ_ϕ2!(mod::model)
    
    updq = 0.25;  #speed to update the price function
    
    ## Initialize asset and repayment grid
    asset_max  = zeros(T+1);
    asset_min  = zeros(T+1);
    mod.asset_grid = zeros(T+1, N_k);
     # Minimum assets (Maximum debt)
    for t = (T-1):-1:2
        income_min   = minimum(mod.Ygrid[t,:,:]);
        asset_min[t] = asset_min[t+1]/(1.0+r) - income_min+ 1e-4;
    end

    for t = 1:T
        income_max     = maximum(mod.Ygrid[t,:,:]);
        asset_max[t+1] = min(asset_max[t] + 0.95*income_max, 100)
    end
    for t = 1:T
        mod.asset_grid[t, :] = LinRange(asset_min[t],asset_max[t],N_k);
    end


    ## Initialize Q0
    Q0 = ones((T,N_z,N_k)).*(1/(1+r));
    for t=1:T
        #a0 = sum( asset_grid[t,:].<0.0 )+1
        for i = 1:N_z
            for api=1:N_k
                Q0[t, i, api] = mod.asset_grid[t,api]<0.0 ? exp(.1*mod.asset_grid[t,api]) / (1.0 + r) : 1.0 / (1.0 + r);
            end
        end
    end


    ϕ_2 = mod.ϕ_2
    println("Starting to compute eqm for ϕ_2 = $ϕ_2")
    println("")

    #++++++++++++++++++++++++++++++++++++++++++
    #allocating stuff:
    ms = sol();
    ms.V = zeros(T+1,N_k,N_ε, N_z);
    ms.C = zeros(T,N_k,N_ε, N_z);
    ms.S = zeros(T,N_k,N_ε, N_z);
    ms.B = zeros(T,N_k,N_ε, N_z);
    ms.Ap = zeros(T,N_k,N_ε, N_z);


    ht = hists();
    ht.a0qtlhist = zeros(Nsim);
    ht.zidxhist  = Array{Int64,2}(undef, Nsim, Tsim + Tburnin);
    ht.εidxhist  = Array{Int64,2}(undef, Nsim, Tsim + Tburnin);
    ht.ahist     = zeros(Nsim,Tsim + Tburnin);
    ht.shist     = zeros(Nsim,Tsim + Tburnin);
    ht.chist     = zeros(Nsim,Tsim + Tburnin);
    ht.bhist     = zeros(Nsim,Tsim + Tburnin);
    ht.agehist   = Array{Int64,2}(undef, Nsim, Tsim + Tburnin);

    

    #= use last Q as a guess or go back to our initial one?
    for qiter = 1:maxiterq
        
        ## Asset and repayment grid : note, this depends on Q
        mod.asset_grid = assetGrid(T, mod.Ygrid, N_k, Q0, mod.asset_grid);
        mina = minimum(mod.asset_grid);
        println("Most possible borrowing: $mina")

        @time backwardSolve!(ms,mod, T, N_k, N_z, N_ε,Q0);

        Q = equilibriumQ(Q0, ms.S, ms.Ap, mod.asset_grid, r, mod.εprobs, mod.zprobs);
        Qdist = maximum( abs.(Q.-Q0) );
        println("sup norm of distance is $Qdist")

        #little diagnostic:
        sinterior = 0
        sbad0 = 0
        sbad0T = 0
        for t=1:T
            for k=1:N_k
                for εi =1:N_ε
                    for zi =1:N_z
                        sinterior = ms.S[t, k, εi, zi]< 0.99 && ms.S[t, k, εi, zi]> 0.01 ? sinterior + 1 : sinterior
                        sbad0     = ms.S[t, k, εi, zi]< 0.99 && mod.asset_grid[t, k]>0.0 ? sbad0 + 1     : sbad0
                        sbad0T    = ms.S[t, k, εi, zi]< 0.99 && mod.asset_grid[t, k]>0.0 && t< T-1 ? sbad0T + 1     : sbad0T
                    end
                end
            end
        end
        #println("The number of bad0 payment shares: $sbad0, $sbad0T were young . Number of interior payment shares: $sinterior")


        #this seems to be an allocation/creation operation rather than element-wise assignment
        # Q0 = deepcopy(Q);
        for t=1:T
            for zi = 1:N_z
                for api = 1:N_k
                    Q0[t,zi,api] = updq * Q[t,zi,api] + (1.0 - updq) * Q0[t,zi,api];
                end
            end
        end
        if Qdist < qtol
            break;
        end
        #if Qdist > .5
        #    updq = .1;
        #end 

    end


    saveloc = string(saveroot,"MPPDD_eqmQ_",ϕ_2,".jld")
    @save saveloc Q0

    ms.Aεz_dist = zeros(T, N_k, N_ε, N_z);
    #steadystateDist!( ms, mod,survivalprobs, T, N_k, N_z, N_ε)

    saveloc = string(saveroot,"MPPDD_solMats_",ϕ_2,".jld")
    @save saveloc ms
    =#
    #ϕ_2 = 1.5;
    saveloc = string(saveroot,"MPPDD_eqmQ_",ϕ_2,".jld")
    @load saveloc Q0
    saveloc = string(saveroot,"MPPDD_solMats_",ϕ_2,".jld")
    @load saveloc ms
    

    println("Drawing shocks!");
    draw_shocks!(ht,mod, 12281951);
    println("Simulating!");
    sim_hists!(ht, ms, mod, T, N_k,N_z, N_ε);
    saveloc = string(saveroot,"MPPDD_simHists_",ϕ_2,".jld")
    @save saveloc ht

end

function setparam_sol()

    # some comparative statics, solving it several times over phi parameters:

    for phi2idx in 1:length(ϕ_2grid)
        
        
        T_peak 	  = 21*freq; # of working years at which the experience profile peaks

        mod = model();
        
        mod.t_rebate = 32 - 25+1; # agent receives rebate at 32 (between 20 and 82)
        mod.rebate = 0.0
        mod.Ygrid, mod.zgrid, mod.εgrid = get_Ygrid(T_retire, T_peak, T, N_z, N_ε, σ_ε, σ_z0, σ_η, ρ); # level grids of income

        mod.zprobs    = get_zprobs(T, N_z, σ_z0, σ_η, ρ); # transition matrices for permanent income
        mod.εprobs    = get_εprobs(T, N_ε, σ_ε);          # transition matrices for transitory income

        mod.ϕ_1 = 0.5; #random, update this
        mod.ϕ_2 = ϕ_2grid[phi2idx];

        solveQ_ϕ2!(mod)
        saveloc = string(saveroot,"MPPDD_modEnvr_",mod.ϕ_2,".jld")
        @save saveloc mod

    end
end


setparam_sol();
#+++++++++++++++++++++++++++++++++
## diagnostic plots:
#=++++++++++++++++++++++++++++++

for phi2idx = 1:length(ϕ_2grid)
    ϕ_2 = ϕ_2grid[phi2idx]

    saveloc = string(saveroot,"MPPDD_eqmQ_",ϕ_2,".jld")
    @load saveloc Q0
    saveloc = string(saveroot,"MPPDD_solMats_",ϕ_2,".jld")
    @load saveloc ms
    saveloc = string(saveroot,"MPPDD_modEnvr_",ϕ_2,".jld")
    @load saveloc mod
    saveloc = string(saveroot,"MPPDD_simHists_",ϕ_2,".jld")
    @load saveloc ht

    for zi=1:N_z
                
        Qages = hcat( Q0[1,zi,:],Q0[15,zi,:],Q0[25,zi,:], Q0[35,zi,:] )
        plot!(hcat(mod.asset_grid[2,:],mod.asset_grid[16,:],mod.asset_grid[26,:],mod.asset_grid[36,:]),Qages,title = "Debt price schedule", ylabel="Equilibrium debt price", xlabel = "a'",
            label=["Age 1" "Age 15" "Age 25" "Age 35"],legend=:bottomright, lw=3, size=(1500,600))
        savefig(string("Q_ages_phi",ϕ_2,"_zi",zi,".png"))

        N_k_neg = sum(mod.asset_grid[2,:].<0)
        Qassets = zeros(T_retire-1, N_k_neg)
        for ai=1:N_k_neg
            for t=2:T_retire
                natborrowing = minimum(mod.asset_grid[t,:]);
                if mod.asset_grid[2,ai] > natborrowing
                    Qassets[t-1,ai] = LinearInterpolation(mod.asset_grid[t,:], Q0[t-1,zi,:])(mod.asset_grid[2,ai]);
                else
                    Qassets[t-1,ai] = NaN
                end
            end
        end
        labs = round.(mod.asset_grid[2,1:N_k_neg]*10.0)/10.0
        plot!(1:(T_retire-1),Qassets,title = "Price as a function of age", ylabel="Equilibrium debt price", xlabel = "Age",
            label= labs', legend=:bottomright, lw=3, size=(1500,600)) #label=["Age 1" "Age 15" "Age 25" "Age 35"],
        savefig(string("Q_assets_phi",ϕ_2,"_zi",zi,".png"))
    
        asset_grid_ages = hcat(mod.asset_grid[2,:],mod.asset_grid[15,:],mod.asset_grid[25,:],mod.asset_grid[35,:]);
        for ei=1:N_ε
            
            Vages = hcat( ms.V[2,:,ei,zi],ms.V[15,:,ei,zi],ms.V[25,:,ei,zi], ms.V[35,:,ei,zi] );
            plot(asset_grid_ages,Vages)
            plot!(title = "Value Function", ylabel="", xlabel = "a",
                label=["Age 2" "Age 15" "Age 25" "Age 35"],legend=:bottomright, lw=3, size=(1500,600))
            savefig(string("V_ages_phi",ϕ_2,"_zi",zi,"_ei",ei,".png"))

            Cages = hcat( ms.C[2,:,ei,zi],ms.C[15,:,ei,zi],ms.C[25,:,ei,zi], ms.C[35,:,ei,zi] );
            plot(asset_grid_ages,Cages)
            plot!(title = "Consumption", ylabel="", xlabel = "a",
                label=["Age 2" "Age 15" "Age 25" "Age 35"],legend=:bottomright, lw=3, size=(1500,600))
            savefig(string("C_ages_phi",ϕ_2,"_zi",zi,"_ei",ei,".png"))

            Sages = hcat( ms.S[2,:,ei,zi],ms.S[15,:,ei,zi],ms.S[25,:,ei,zi], ms.S[35,:,ei,zi] );
            plot(asset_grid_ages,Sages)
            plot!(title = "Delinquency policy", ylabel="", xlabel = "a",
                label=["Age 2" "Age 15" "Age 25" "Age 35"],legend=:bottomright, lw=3, size=(1500,600))
            savefig(string("S_ages_phi",ϕ_2,"_zi",zi,"_ei",ei,".png"))

            Apages = hcat( ms.Ap[2,:,ei,zi]-mod.asset_grid[2,:],ms.Ap[15,:,ei,zi]-mod.asset_grid[15,:],ms.Ap[25,:,ei,zi]-mod.asset_grid[25,:], ms.Ap[35,:,ei,zi]-mod.asset_grid[35,:] );
            plot(asset_grid_ages,Apages)
            plot!(title = "A' policy", ylabel="a'-a", xlabel = "a",
                label=["Age 2" "Age 15" "Age 25" "Age 35"],legend=:bottomright, lw=3, size=(1500,600))
            savefig(string("Ap_ages_phi",ϕ_2,"_zi",zi,"_ei",ei,".png"))
            
            Adistages = hcat( ms.Aεz_dist[2,:,ei,zi],ms.Aεz_dist[15,:,ei,zi],ms.Aεz_dist[25,:,ei,zi], ms.Aεz_dist[35,:,ei,zi] );
            plot(asset_grid_ages,Adistages)
            plot!(title = "A distribution", ylabel="Density", xlabel = "a",
                label=["Age 2" "Age 15" "Age 25" "Age 35"],legend=:bottomright, lw=3, size=(1500,600))
            savefig(string("Adist_ages_phi",ϕ_2,"_zi",zi,"_ei",ei,".png"))

        end

        

        # experiment giving transitory income shock: 
        Apages = hcat((ms.Ap[2, :,3,zi]- ms.Ap[2, :,2,zi])./ (mod.Ygrid[2, 3, zi]  - mod.Ygrid[2 , 2, zi]),
                      (ms.Ap[15,:,3,zi]- ms.Ap[15,:,2,zi])./ (mod.Ygrid[15, 3, zi] - mod.Ygrid[15, 2, zi]),
                      (ms.Ap[25,:,3,zi]- ms.Ap[25,:,2,zi])./ (mod.Ygrid[25, 3, zi] - mod.Ygrid[25, 2, zi]),
                      (ms.Ap[35,:,3,zi]- ms.Ap[35,:,2,zi])./ (mod.Ygrid[35, 3, zi] - mod.Ygrid[35, 2, zi]) )
        epschng = round(mod.Ygrid[2, 3, zi]*100.0 - mod.Ygrid[2, 2, zi]*100.0)/100.0
        plot(asset_grid_ages,Apages)
        plot!(title = "A' change with ε shock of $epschng", ylabel="a' change / ε shock ", xlabel = "a",
            label=["Age 2" "Age 15" "Age 25" "Age 35"],legend=:bottomright, lw=3, size=(1500,600))
        savefig(string("ApDif_ages_phi",ϕ_2,"_zi",zi,".png"))
        
        delQages = zeros(N_k,4 )
        for ai=1:N_k
            delQages[ai,1] = (LinearInterpolation(mod.asset_grid[3,:], Q0[2,zi,:])(ms.Ap[2, ai,3,zi]) - LinearInterpolation(mod.asset_grid[3,:], Q0[2,zi,:])(ms.Ap[2, ai,2,zi]))/
                (ms.Ap[2, ai,3,zi]- ms.Ap[2, ai,2,zi])
            delQages[ai,2] = (LinearInterpolation(mod.asset_grid[16,:], Q0[15,zi,:])(ms.Ap[15, ai,3,zi]) - LinearInterpolation(mod.asset_grid[16,:], Q0[15,zi,:])(ms.Ap[15, ai,2,zi]))/
                (ms.Ap[15, ai,3,zi]- ms.Ap[15, ai,2,zi])
            delQages[ai,3] = (LinearInterpolation(mod.asset_grid[26,:], Q0[25,zi,:])(ms.Ap[25, ai,3,zi]) - LinearInterpolation(mod.asset_grid[26,:], Q0[25,zi,:])(ms.Ap[25, ai,2,zi]))/
                (ms.Ap[2, ai,3,zi]- ms.Ap[2, ai,2,zi])
            delQages[ai,4] = (LinearInterpolation(mod.asset_grid[36,:], Q0[35,zi,:])(ms.Ap[35, ai,3,zi]) - LinearInterpolation(mod.asset_grid[36,:], Q0[35,zi,:])(ms.Ap[35, ai,2,zi]))/
                (ms.Ap[15, ai,3,zi]- ms.Ap[15, ai,2,zi])
        end
        epschng = round(mod.Ygrid[2, 3, zi]*100.0 - mod.Ygrid[2, 2, zi]*100.0)/100.0
        plot(asset_grid_ages,delQages)
        plot!(title = "dQdA change with ε shock of $epschng", ylabel="q change / ε shock ", xlabel = "a", label=["Age 2" "Age 15" "Age 25" "Age 35"],
            legend=:bottomright, lw=3, size=(1500,600))
        savefig(string("QDif_phi",ϕ_2,"_zi",zi,".png"))
    
        

        # experiment giving transitory income shock: 
        Cages = hcat((ms.C[2, :,3,zi]- ms.C[2, :,2,zi])./ (mod.Ygrid[2, 3, zi]  - mod.Ygrid[2 , 2, zi]),
        (ms.C[15,:,3,zi]- ms.C[15,:,2,zi])./ (mod.Ygrid[15, 3, zi] - mod.Ygrid[15, 2, zi]),
        (ms.C[25,:,3,zi]- ms.C[25,:,2,zi])./ (mod.Ygrid[25, 3, zi] - mod.Ygrid[25, 2, zi]),
        (ms.C[35,:,3,zi]- ms.C[35,:,2,zi])./ (mod.Ygrid[35, 3, zi] - mod.Ygrid[35, 2, zi]) )
        epschng = round(mod.Ygrid[2, 3, zi]*100.0 - mod.Ygrid[2, 2, zi]*100.0)/100.0
        plot(asset_grid_ages,Cages)
        plot!(title = "C change with ε shock of $epschng", ylabel="c change / ε shock ", xlabel = "a",
            label=["Age 2" "Age 15" "Age 25" "Age 35"],legend=:bottomright, lw=3, size=(1500,600))
        savefig(string("CDif_ages_phi",ϕ_2,"_zi",zi,".png"))
        
        Bages = hcat((ms.B[2, :,3,zi]- ms.B[2, :,2,zi])./ (mod.Ygrid[2, 3, zi]  - mod.Ygrid[2 , 2, zi]),
        (ms.B[15,:,3,zi]- ms.B[15,:,2,zi])./ (mod.Ygrid[15, 3, zi] - mod.Ygrid[15, 2, zi]),
        (ms.B[25,:,3,zi]- ms.B[25,:,2,zi])./ (mod.Ygrid[25, 3, zi] - mod.Ygrid[25, 2, zi]),
        (ms.B[35,:,3,zi]- ms.B[35,:,2,zi])./ (mod.Ygrid[35, 3, zi] - mod.Ygrid[35, 2, zi]) )
        epschng = round(mod.Ygrid[2, 3, zi]*100.0 - mod.Ygrid[2, 2, zi]*100.0)/100.0
        plot(asset_grid_ages,Bages)
        plot!(title = "B change with ε shock of $epschng", ylabel="b change / ε shock ", xlabel = "a",
            label=["Age 2" "Age 15" "Age 25" "Age 35"],legend=:bottomright, lw=3, size=(1500,600))
        savefig(string("BDif_ages_phi",ϕ_2,"_zi",zi,".png"))

        
        Sages = hcat((ms.S[2, :,3,zi]- ms.S[2, :,2,zi]),
        (ms.S[15,:,3,zi]- ms.S[15,:,2,zi]),
        (ms.S[25,:,3,zi]- ms.S[25,:,2,zi]),
        (ms.S[35,:,3,zi]- ms.S[35,:,2,zi]) )
        epschng = round(mod.Ygrid[2, 3, zi]*100.0 - mod.Ygrid[2, 2, zi]*100.0)/100.0
        plot(asset_grid_ages,Sages)
        plot!(title = "S change with ε shock of $epschng", ylabel="s change / ε shock ", xlabel = "a",
            label=["Age 2" "Age 15" "Age 25" "Age 35"],legend=:bottomright, lw=3, size=(1500,600))
        savefig(string("SDif_ages_phi",ϕ_2,"_zi",zi,".png"))
             
    end
end
=#