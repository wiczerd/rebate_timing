#Pkg.add("MAT")
#Pkg.add("JLD")

#nworkers = 3;  # num workers -- REQUEST WORKERS INSIDE CODE, NOT IN TERMINAL
#mem      = 5; # request memory for each worker

#ENV["frbnyjuliamemory"]=string(mem)*"G"
#addprocs_frbny(nworkers)
cd("/home/david/workspace/rebate_timing/MPPDD_Julia/")
#@everywhere begin
#@everywhere begin
    #using JLD2
    using JLD
    using MAT
#    using Distributed
#end

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


#function separ()

#ss = string("MPPDD_solMats_1.5.jld")
#@load ss ms
#@load "MPPDD_solMats_1.5.jld"

#=

@load "MPPDD_solMats_1.5.jld"
matf = matopen("MPPDD_solMatsC_1.5.mat","w")
write(matf,"C",ms.C)
close(matf)

matf = matopen("MPPDD_solMatsB_1.5.mat","w")
write(matf,"B",ms.B)
close(matf)

matf = matopen("MPPDD_solMatsAp_1.5.mat","w")
write(matf,"Ap",ms.Ap)
close(matf)

matf = matopen("MPPDD_solMatsS_1.5.mat","w")
write(matf,"S",ms.S)
close(matf)


@load "MPPDD_eqmQ_1.5.jld"  Q0

matf = matopen("MPPDD_eqmQ_1.5.mat","w")
write(matf,"Q",Q0)
close(matf)

@load "MPPDD_modEnvr_1.5.jld"  mod


matf = matopen("MPPDD_modEnvrAssetGrid_1.5.mat","w")
write(matf,"asset_grid",mod.asset_grid)
close(matf)

matf = matopen("MPPDD_modEnvrYgrid_1.5.mat","w")
write(matf,"Ygrid",mod.Ygrid)
close(matf)



=#

@load "MPPDD_simHists_1.5.jld"  ht


matf = matopen("MPPDD_simHistsAhist_1.5.mat","w")
write(matf,"ahist",ht.ahist)
close(matf)

matf = matopen("MPPDD_simHistsChist_1.5.mat","w")
write(matf,"chist",ht.chist)
close(matf)

matf = matopen("MPPDD_simHistsShist_1.5.mat","w")
write(matf,"shist",ht.shist)
close(matf)

matf = matopen("MPPDD_simHistsBhist_1.5.mat","w")
write(matf,"bhist",ht.bhist)
close(matf)

matf = matopen("MPPDD_simHistsAgehist_1.5.mat","w")
write(matf,"agehist",ht.agehist)
close(matf)

matf = matopen("MPPDD_simHistsEpshist_1.5.mat","w")
write(matf,"epsidxhist",ht.εidxhist)
close(matf)

matf = matopen("MPPDD_simHistsZhist_1.5.mat","w")
write(matf,"zidxhist",ht.zidxhist)
close(matf)

