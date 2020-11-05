# run line: JULIA_NUM_THREADS=32 julia main.jl


## Packages

#using Pkg
#include("import.jl") # uncomment to add required packages

using LinearAlgebra
using SparseArrays
using Arpack
using Interpolations
using Distributions
using Optim
using JLD
using Plots

import Base.Threads.@spawn
#cd("/Users/meghanagaur/DropBox/Consumption_Model_Codes/julia/")     # personal
#cd("C:\\Users\\rcemxg23\\DropBox\\Consumption_Model_Codes\\julia") # ran virtual desktop
cd("/home/david/workspace/rebate_timing/MPPDD_Julia/")
#cd("/Users/rcedxm19/Dropbox/MPPDD/codes/MPPDD_Julia")

saveroot = "/home/david/workspace/rebate_timing/MPPDD_Julia/"

include("functions/incomeProcess.jl")
include("functions/agent.jl")

qtol         = 1e-5; # tolerance for convergence in Q
maxiterq     = 80;   # how many iterations?
updq         = 0.5;  #speed to update the price function

## Parameters
freq         = 1; # 1=annual, 4=quarterly
r            = (1 + 0.03)^(1/freq) - 1;
β            = 0.95^(1/freq);
γ            = 1;
ϕ_1          = 0.5; #random, update this
ϕ_2          = 1.5; #random, update this

## Agents begin at  age 25, die at 95, retire at 65
t0           = 25; #25;
T            = 95; #95;
T_retire     = 65; #65;

"""
Convert from years to index
"""
function convertt(tin::Int64, freq::Int64, t0::Int64)
        return tout = (tin - t0)*freq + 1
end

T 	         = convertt(T, freq,t0);
T_retire     = convertt(T_retire, freq,t0);

rebate       = 0.0;  # added to BC at t_reb (argument in backwardSolve)
cbar         = 0.0;  # consumption floor

# Survival Probabilities
survivalprobs = ones(T); # T set to 0 in agent.jl

## State space

N_k          = 50;      # Number of asset states, 40-50
N_z          = 7;      # Number of permanent income states, 7-10
N_ε          = 5;	    # Number of transitory income states, 5-21
N_s          = 25;      # Number of delinquency grid points, 10-5

## Income Process  - from Kaplan & Violante (2010)
ρ	      =  1.0;
σ_η       = (0.01/freq)^(1/2);  # Variance of permanent income shocks = 0.01
σ_z0      = (0.15/freq)^(1/2);  # Variance of permanent income shocks = 0.15
σ_ε       = (0.05/freq)^(1/2);  # Variance of transitory income shocks = 0.05


T_peak 	  = 21*freq; # of working years at which the experience profile peaks

Ygrid, Zgrid, εgrid = get_Ygrid(T_retire, T_peak, T, N_z, N_ε, σ_ε, σ_z0, σ_η, ρ); # level grids of income

zprobs    = get_zprobs(T, N_z, σ_z0, σ_η, ρ); # transition matrices for permanent income
εprobs    = get_εprobs(T, N_ε, σ_ε);          # transition matrices for transitory income


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
    N_k::Int64
    N_s::Int64
    N_z::Int64
    N_ε::Int64

    model() = new()
end



## Asset and repayment grid

asset_grid = assetGrid(T, Ygrid, N_k, r);
s_grid     = collect(LinRange(0, 1, N_s));

## Solve model

t_reb    = 32; # actual age that agent receives rebate (between 20 and 82)
t_rebate = convertt(t_reb, freq,t0); # index


Q0 = ones((T,N_z,N_k)).*(1/(1+r));
for t=1:T
    #a0 = sum( asset_grid[t,:].<0.0 )+1
    for i = 1:N_z
        for api=1:N_k
            Q0[t, i, api] = asset_grid[t,api]<0.0 ? exp(asset_grid[t,api]) / (1.0 + r) : 1.0 / (1.0 + r);
        end
    end
end


# some comparative statics, solving it several times over phi parameters:
ϕ_2grid = [1.1 1.5 5.0]
for phi2idx in 1:length(ϕ_2grid)

    ϕ_2 = ϕ_2grid[phi2idx];

    

    println("Starting to compute eqm for ϕ_2 = $ϕ_2")

    #use last Q as a guess or go back to our initial one?
    for qiter = 1:maxiterq
        @time V, B, C, S, Ap = backwardSolve(rebate, t_rebate, β, γ, ϕ_1, ϕ_2, r,
            εprobs, zprobs, Ygrid, s_grid, asset_grid, T, N_k, N_s, N_z, N_ε, cbar,Q0);

        Q = equilibriumQ(Q0, S, Ap, asset_grid, r, εprobs, zprobs);
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
                        sinterior = S[t, k, εi, zi]< 0.99 && S[t, k, εi, zi]> 0.01 ? sinterior + 1 : sinterior
                        sbad0     = S[t, k, εi, zi]< 0.99 && asset_grid[t, k]>0.0 ? sbad0 + 1     : sbad0
                        sbad0T    = S[t, k, εi, zi]< 0.99 && asset_grid[t, k]>0.0 && t< T-1 ? sbad0T + 1     : sbad0T
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
    end


    saveloc = string("/home/david/workspace/rebate_timing/MPPDD_Julia/","MPPDD_eqmQ_",ϕ_2,".jld")
    @save saveloc Q0
    #saveloc = string(saveroot,"MPPDD_solMats_",ϕ_2,".jld")
    #@save saveloc msol

end

#+++++++++++++++++++++++++++++++++
## diagnostic plots:
#+++++++++++++++++++++++++++++++

ϕ_2grid = [1.1 1.5 5.0]
for phi2idx = 1:3
    ϕ_2 = ϕ_2grid[phi2idx]

    saveloc = string(saveroot,"MPPDD_eqmQ_",ϕ_2,".jld")
    @load saveloc Q0
#    saveloc = string(saveroot,"MPPDD_solMats_",ϕ_2,".jld")
#    @load saveloc V B C S Ap

    @time V, B, C, S, Ap = backwardSolve(rebate, t_rebate, β, γ, ϕ_1, ϕ_2, r,
        εprobs, zprobs, Ygrid, s_grid, asset_grid, T, N_k, N_s, N_z, N_ε, cbar,Q0);

    Aεz_dist = zeros(T, N_k, N_ε, N_z);
    steadystateDist!( Aεz_dist, Ap, survivalprobs, T, N_k, N_s, N_z, N_ε, asset_grid, εprobs, zprobs)
    saveloc = string(saveroot,"MPPDD_eqmDist_",ϕ_2,".jld")
    @save saveloc Aεz_dist
    
    for zi=1:N_z
zi=4
        Qages = hcat( Q0[1,zi,:],Q0[15,zi,:],Q0[25,zi,:], Q0[35,zi,:] )
        plot(hcat(asset_grid[2,:],asset_grid[16,:],asset_grid[26,:],asset_grid[36,:]),Qages,title = "Debt price schedule", ylabel="Equilibrium debt price", xlabel = "a'",
            label=["Age 1" "Age 15" "Age 25" "Age 35"],legend=:bottomright, lw=3, size=(1500,600))
        savefig(string("Q_ages_phi",ϕ_2,"_zi",zi,".png"))

        N_k_neg = sum(asset_grid[2,:].<0)
        Qassets = zeros(T_retire-1, N_k_neg)
        for ai=1:N_k_neg
            for t=2:T_retire
                natborrowing = minimum(asset_grid[t,:]);
                if asset_grid[2,ai] > natborrowing
                    Qassets[t-1,ai] = LinearInterpolation(asset_grid[t,:], Q0[t-1,zi,:])(asset_grid[2,ai]);
                else
                    Qassets[t-1,ai] = NaN
                end
            end
        end
        labs = round.(asset_grid[2,1:N_k_neg]*10.0)/10.0
        plot(1:(T_retire-1),Qassets,title = "Price as a function of age", ylabel="Equilibrium debt price", xlabel = "Age",
            label= labs', legend=:bottomright, lw=3, size=(1500,600)) #label=["Age 1" "Age 15" "Age 25" "Age 35"],
        savefig(string("Q_assets_phi",ϕ_2,"_zi",zi,".png"))


        for ei=1:N_ε

            Sages = hcat( S[2,:,ei,zi],S[15,:,ei,zi],S[25,:,ei,zi], S[35,:,ei,zi] );
            plot(hcat(asset_grid[2,:],asset_grid[15,:],asset_grid[25,:],asset_grid[35,:]),Sages,title = "Delinquency policy", ylabel="Optimal s", xlabel = "a",
                label=["Age 2" "Age 15" "Age 25" "Age 35"],legend=:bottomright, lw=3, size=(1500,600))
            savefig(string("S_ages_phi",ϕ_2,"_zi",zi,"_ei",ei,".png"))

            Apages = hcat( Ap[2,:,ei,zi]-asset_grid[2,:],Ap[15,:,ei,zi]-asset_grid[15,:],Ap[25,:,ei,zi]-asset_grid[25,:], Ap[35,:,ei,zi]-asset_grid[35,:] )
            plot(hcat(asset_grid[2,:],asset_grid[15,:],asset_grid[25,:],asset_grid[35,:]),Apages,title = "A' policy", ylabel="a'-a", xlabel = "a",
                label=["Age 2" "Age 15" "Age 25" "Age 35"],legend=:bottomright, lw=3, size=(1500,600))
            savefig(string("Ap_ages_phi",ϕ_2,"_zi",zi,"_ei",ei,".png"))

            Distages = hcat( Aεz_dist[2,:,ei,zi],Aεz_dist[15,:,ei,zi],Ap[25,:,ei,zi], Aεz_dist[35,:,ei,zi])
            plot(hcat(asset_grid[2,:],asset_grid[15,:],asset_grid[25,:],asset_grid[35,:]),
                Distages,title = "Distribution of Assets", ylabel="Density(a)", xlabel = "a",
                label=["Age 2" "Age 15" "Age 25" "Age 35"],legend=:bottomright, lw=3, size=(1500,600))
            savefig(string("Adist_ages_phi",ϕ_2,"_zi",zi,"_ei",ei,".png"))
        end


        # experiment giving transitory income shock: 
        Apages = hcat( (Ap[2, :,3,zi]- Ap[2, :,2,zi])./ (Ygrid[2, 3, zi]  - Ygrid[2 , 2, zi]),
                       (Ap[15,:,3,zi]- Ap[15,:,2,zi])./ (Ygrid[15, 3, zi] - Ygrid[15, 2, zi]),
                       (Ap[25,:,3,zi]- Ap[25,:,2,zi])./ (Ygrid[25, 3, zi] - Ygrid[25, 2, zi]),
                       (Ap[35,:,3,zi]- Ap[35,:,2,zi])./ (Ygrid[35, 3, zi] - Ygrid[35, 2, zi]) )
        epschng = round(εgrid[2,N_ε]/Ygrid[2, 1, zi]*100.0 - εgrid[2,1]/Ygrid[2, 1, zi]*100.0)/100.0
        plot(hcat(asset_grid[2,:],asset_grid[15,:],asset_grid[25,:],asset_grid[35,:]),Apages,title = "A' change with ε shock of $epschng", ylabel="a' change / ε shock ", xlabel = "a",
            label=["Age 2" "Age 15" "Age 25" "Age 35"],legend=:bottomright, lw=3, size=(1500,600))
        savefig(string("ApDif_ages_phi",ϕ_2,"_zi",zi,".png"))

    end
end
