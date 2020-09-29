"""
Grid of permanent income
"""
function get_zgrid(T::Int64, N_z::Int64, σ_z0::Float64, σ_η::Float64, ρ::Float64)
     zgrid = zeros(T + 2, N_z);
     σ     = ones(T+1)*σ_z0;

     for t = 1:T
         zgrid[t,:] = collect(LinRange(-3.5*σ[t], 3.5*σ[t], N_z))
         σ[t+1]     = (ρ^2)*(σ[t]^2) + σ_η^2
     end

     return zgrid
end

"""
Grid of log transitory income shock
"""
function get_εgrid(N_ε::Int64, σ_ε::Float64)
  return collect(LinRange(-3.5*σ_ε, 3.5*σ_ε, N_ε))
end

"""
Define a pension function
"""
function pension(x::Float64)
    return 0.5*x
    # Retirement income (pension) in years after retirement = 1/2
    # in their final working period.
    # simple for now, can define this more formally, as in K & V 2010
end

"""
Return deterministic experience profile of log income, common to all households
"""
function pincome_profile(T_retire::Int64, T_peak::Int64)
    years_worked = cumsum(ones(T_retire)) .-1;
    return (2 .- ((years_worked .- T_peak).^2 ./ T_peak^2))
end

"""
Return level grid of total income
"""
function get_Ygrid(T_retire::Int64, T_peak::Int64, T::Int64, N_z::Int64,
    N_ε::Int64, σ_ε::Float64, σ_z0::Float64, σ_z::Float64, ρ::Float64)

    # Log components
    zgrid        = get_zgrid(T, N_z, σ_z0, σ_η, ρ);
    κ            = pincome_profile(T_retire, T_peak);
    ε            = get_εgrid(N_ε, σ_ε);

    # Level grids
    Ygrid        = zeros(T, N_ε, N_z);
    Zgrid        = zeros(T, N_z);
    εgrid        = zeros(T, N_ε);

    for t = 1:T
        for i = 1:N_z
            for j = 1:N_ε
                if t >= T_retire
                    Ygrid[t, j, i]    = pension(exp(zgrid[T_retire - 1, i]));
                    Zgrid[t, i]       = pension(exp(κ[T_retire - 1] + zgrid[T_retire - 1, i]));
                    εgrid[t, j]       = 1;
                else
                    Ygrid[t, j, i]    = exp(κ[t] + zgrid[t, i] + ε[j]);
                    Zgrid[t, i]       = exp(zgrid[t, i]);
                    εgrid[t, j]       = exp(ε[j]);
                end
            end
        end
    end
    return Ygrid, Zgrid, εgrid
end

function get_zprobs(T::Int64, N_z::Int64, σ_z0::Float64, σ_η::Float64, ρ::Float64)

    zgrid        = get_zgrid(T, N_z, σ_z0, σ_η, ρ);
    tprobs       = zeros(T, N_z, N_z);

    function Pr(x)
        d = Normal(0, σ_η)
        return cdf(d, x)
    end

    # Get transition matrices over time
    for t = 1:T
        if t < T_retire - 1
            zp       = zgrid[t + 1,:];
            z0       = zgrid[t,:];
            for j = 1:N_z
                for i = 1:N_z
                    # Interior points
                    if j > 1 && j < N_z
                        # Probability of moving within range of zp[j] in the next period, given being in z0[i] in the current period
                        tprobs[t, i, j] = Pr(0.5*(zp[j] + zp[j + 1]) - ρ*z0[i]) - Pr(0.5*(zp[j - 1] + zp[j]) - ρ*z0[i])
                    elseif j == 1
                        # Transition to minimum
                        tprobs[t, i, j] = Pr(0.5*(zp[j] + zp[j + 1]) - ρ*z0[i])
                    else
                        # Transitions to maximum
                        tprobs[t, i, j] = 1 - Pr(0.5*(zp[j - 1] + zp[j]) - ρ*z0[i])
                    end
                end
            end
        else
            # Retirement transition matrix is just identity
            tprobs[t, :,:] .= Matrix(1.0I, N_z, N_z);
        end
    end
    return tprobs
end


function get_εprobs(T::Int64, N_ε::Int64, σ_ε::Float64)

    # agent does not actually face transitory shocks after retirement
    # but ignore this fact for now.

    ε      = get_εgrid(N_ε, σ_ε);
    tprobs = zeros(T, N_ε);
    idx0   = argmin(ε.^2);

    function Pr(x)
        d = Normal(0, σ_ε)
        return cdf(d, x)
    end


    for j = 1:N_ε
        # Interior points
        if j > 1 && j < N_ε
            # Probability of moving within range of tshock[j] in the next period
            tprobs[:, j] .= Pr(0.5*(ε[j] + ε[j + 1]) ) - Pr(0.5*(ε[j - 1] + ε[j]));
        elseif j == 1
            # Transition to minimum
            tprobs[:, j] .= Pr(0.5*(ε[j] + ε[j + 1]) );
        else
            # Transitions to maximum
            tprobs[:, j] .= 1 - Pr(0.5*(ε[j - 1] + ε[j]) );
        end
    end

    #= not necessary
    tprobs[T_retire - 1, :, idx0] .= 1

    for t = T_retire:T
        tprobs[t, :, :]    .=  Matrix(1.0I, N_ε, N_ε);
    end
    =#

    return tprobs
end


function get_εprobs_trx(T::Int64, N_ε::Int64, σ_ε::Float64)

    # agent does not actually face transitory shocks after retirement
    # but ignore this fact for now.

    ε      = get_εgrid(N_ε, σ_ε);
    tprobs = zeros(T, N_ε, N_ε);
    idx0   = argmin(ε.^2);

    function Pr(x)
        d = Normal(0, σ_ε)
        return cdf(d, x)
    end

    for i = 1:N_ε
        for j = 1:N_ε
            # Interior points
            if j > 1 && j < N_ε
                # Probability of moving within range of tshock[j] in the next period
                tprobs[:, i, j] .= Pr(0.5*(ε[j] + ε[j + 1]) - ε[i]) - Pr(0.5*(ε[j - 1] + ε[j]) - ε[i]);
            elseif j == 1
                # Transition to minimum
                tprobs[:, i, j] .= Pr(0.5*(ε[j] + ε[j + 1]) - ε[i]);
            else
                # Transitions to maximum
                tprobs[:, i, j] .= 1 - Pr(0.5*(ε[j - 1] + ε[j]) - ε[i]);
            end
        end
    end

    #= not necessary
    tprobs[T_retire - 1, :, idx0] .= 1

    for t = T_retire:T
        tprobs[t, :, :]    .=  Matrix(1.0I, N_ε, N_ε);
    end
    =#

    return tprobs
end
