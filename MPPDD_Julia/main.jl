## Packages

using Pkg
#include("import.jl") # uncomment to add required packages

using LinearAlgebra
using Interpolations
using Distributions
using BenchmarkTools
using Optim
using JLD2
import Base.Threads.@spawn
#cd("/Users/meghanagaur/DropBox/Consumption_Model_Codes/julia/")     # personal
#cd("C:\\Users\\rcemxg23\\DropBox\\Consumption_Model_Codes\\julia") # ran virtual desktop
cd("/home/david/workspace/rebate_timing/MPPDD_Julia/")

include("functions/incomeProcess.jl")
include("functions/agent.jl")

## Parameters
freq         = 1; # 1=annual, 4=quarterly
r            = (1 + 0.03)^(1/freq) - 1;
β            = 0.95^(1/freq);
γ            = 1;
ϕ_1          = 0.5; #random, update this
ϕ_2          = 1.5; #random, update this

## Agents begin at  age 25, die at 95, retire at 65
t0           = 25;
T            = 95;
T_retire     = 65;

"""
Convert from years to index
"""
function convertt(tin::Int64, freq::Int64)
        return tout = (tin - 20)*freq + 1
end

T 	     = convertt(T, freq);
T_retire     = convertt(T_retire, freq);

rebate       = 0.0;  # added to BC at t_reb (argument in backwardSolve)
cbar         = 0.0;  # conusmption floor

# Survival Probabilities
survivalprobs = ones(T); # T set to 0 in agent.jl

## State space

N_k          = 40;      # Number of asset states, 50
N_z          = 7;      # Number of permanent income states, 10
N_ε          = 7;	    # Number of transitory income states, 21
N_s          = 6;       # Number of delinquency grid points, 5

## Income Process  - from Kaplan & Violante (2010)
ρ	      =  1.0;
σ_η       = (0.01/freq)^(1/2);  # Variance of permanent income shocks = 0.01
σ_z0      = (0.15/freq)^(1/2);  # Variance of permanent income shocks = 0.15
σ_ε       = (0.05/freq)^(1/2);  # Variance of transitory income shocks = 0.05

survivalprobs = ones(T)         # demographics

T_peak 	  = 21*freq; # of working years at which the experience profile peaks

Ygrid, Zgrid, εgrid = get_Ygrid(T_retire, T_peak, T, N_z, N_ε, σ_ε, σ_z0, σ_η, ρ); # level grids of income

zprobs    = get_zprobs(T, N_z, σ_z0, σ_η, ρ); # transition matrices for permanent income
εprobs    = get_εprobs(T, N_ε, σ_ε);          # transition matrices for transitory income

#= NOTE: the below lines define income grids for debugging;
# I will add incomeProcess.jl later to produce the actual grids + transition matrices

# Permanent income grid and transition matrix
M       = [0.8 0.15 0.05; 0.10 0.8 0.10; 0.05 0.15 0.80];
zprobs  = zeros(T,3,3);

for t = 1:T_ret-1
        zprobs[t,:,:]  = M;
end

for t = T_ret:T
        zprobs[t,:,:] = Matrix(0.1I, 3,3);
end

zgrid  = zeros(T, 3)

for t = 1:T_ret-1
        zgrid[t,:] = [0.10 0.20 0.30];
end

for t = T_ret:T
        zgrid[t,:]  = [0.10 0.20 0.30]*0.5;
end

# Transitory income grid
εgrid  = zeros(T, 1);
εprobs = ones(T, 1, 1)

# Total income = permanent income
Ygrid         = zeros(T, 1, 3);
Ygrid[:,1,:] .= z_grid

N_z          = size(Ygrid,3);     # Number of permanent income states
N_ε          = size(Ygrid,2);	# Number of transitory income states
=#

## Asset and repayment grid

asset_grid = assetGrid(T, Ygrid, N_k, r);
s_grid     = collect(LinRange(0, 1, N_s));

## Solve model

t_reb    = 32; # actual age that agent receives rebate (between 20 and 82)
t_rebate = convertt(t_reb, freq); # index


# two possible guess for Q
Q0 = ones((T,N_z,N_k)).*(1/(1+r));
Q1 = ones((T,N_z,N_k)).*(1/(1+r));
for t=1:T
    #a0 = sum( asset_grid[t,:].<0.0 )+1
    for i = 1:N_z
        for api=1:N_k
            Q0[t, i, api] = asset_grid[t,api]<0.0 ? exp(asset_grid[t,api]) / (1.0 + r) : 1.0 / (1.0 + r);
        end
    end
end


for qiter = 1:40
    V, B, C, S, Ap = backwardSolve(rebate, t_rebate, β, γ, ϕ_1, ϕ_2, r,
        εprobs, zprobs, Ygrid, s_grid, asset_grid, T, N_k, N_s, N_z, N_ε, cbar,Q0);

    Q = equilibriumQ(Q0, S, Ap, asset_grid, r, εprobs, zprobs);
    Qdist = maximum( abs.(Q.-Q0) );
    println("sup norm of distance is $Qdist")

    #this seems to be an allocation/creation operation rather than element-wise assignment
    # Q0 = deepcopy(Q);
    for t=1:T
        for i = 1:N_z
            for api = 1:N_k
                Q0[t,i,api] = Q[t,i,api];
            end
        end
    end
end
@save "MPPDD_eqmq_newphi.jld2" Q0

#quick plot for check:
using Plots
Q1 = hcat( Q0[1,1,:],Q0[1,3,:],Q0[1,6,:], )
default(titlefont = (14, "times"), legendfont = (10,"times"), guidefont = (12, :darkgreen,"times"), tickfont = (12, :orange), framestyle = :zerolines, yminorgrid = true)
plot(asset_grid[2,:],Q1,title = "Debt price schedule", ylabel="Equilibrium debt price", xlabel = "a'",
    label=["Low z" "Medium z" "High z"],legend=:bottomright, lw=3)
savefig("Q_t1.png")

D1 = hcat( D[2,:,1,3],D[2,:,3,3],D[2,:,6,3], )
default(titlefont = (14, "times"), legendfont = (10,"times"), guidefont = (12, :darkgreen,"times"), tickfont = (12, :orange), framestyle = :zerolines, yminorgrid = true)
plot(asset_grid[2,:],D1,title = "Debt price schedule", ylabel="D", xlabel = "a",
    label=["Low z" "Medium z" "High z"],legend=:bottomright, lw=3)


Q10 = hcat( Q0[10,1,:],Q0[10,3,:],Q0[10,6,:], )
plot(asset_grid[10,:],Q10,title = "Debt price schedule", ylabel="Equilibrium debt price", xlabel = "a'",
    label=["Low z" "Medium z" "High z"],legend=:bottomright, lw=3)

Q20 = hcat( Q0[20,1,:],Q0[20,3,:],Q0[20,6,:], )
plot(asset_grid[20,:],Q20,title = "Debt price schedule", ylabel="Equilibrium debt price", xlabel = "a'",
    label=["Low z" "Medium z" "High z"],legend=:bottomright, lw=3)

Qages = hcat( Q0[1,3,:],Q0[15,3,:],Q0[25,3,:], Q0[35,3,:] )
plot(hcat(asset_grid[2,:],asset_grid[16,:],asset_grid[26,:],asset_grid[36,:]),Qages,title = "Debt price schedule", ylabel="Equilibrium debt price", xlabel = "a'",
    label=["Age 1" "Age 15" "Age 25" "Age 35"],legend=:bottomright, lw=3, size=(1500,600))
savefig("Q_ages.png")

Sages = hcat( S[2,:,3,3],S[15,:,3,3],S[25,:,3,3], S[35,:,3,3] )
plot(hcat(asset_grid[2,:],asset_grid[15,:],asset_grid[25,:],asset_grid[35,:]),Sages,title = "Delinquency policy", ylabel="Optimal d", xlabel = "a",
    label=["Age 1" "Age 15" "Age 25" "Age 35"],legend=:bottomright, lw=3, size=(1500,600))
savefig("S_ages_newphi.png")

Apages = hcat( Ap[2,:,3,3]-asset_grid[2,:],Ap[15,:,3,3]-asset_grid[15,:],Ap[25,:,3,3]-asset_grid[25,:], Ap[35,:,3,3]-asset_grid[35,:] )
plot(hcat(asset_grid[2,:],asset_grid[15,:],asset_grid[25,:],asset_grid[35,:]),Apages,title = "A' policy", ylabel="a'-a", xlabel = "a",
    label=["Age 1" "Age 15" "Age 25" "Age 35"],legend=:bottomright, lw=3, size=(1500,600))
