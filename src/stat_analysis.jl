using Random
using Distributions
import AffineInvariantMCMC
using DelimitedFiles
include("super_rad.jl")
include("tde_input.jl")
using Dates
# using PyCall


# Random.seed!(1234)


LowMass = false
input_data = "Andy" # Doddy, Masha, or Andy
alpha_max_cut = 1.0
alpha_min_cut = 0.01

Ftag = "_"
tau_max_override = nothing

numwalkers=10
thinning=1
numsamples_perwalker=2000
burnin=500


if LowMass
    ###### Low Mass Region
    lg_m_low = -13
    lg_m_high = -10
    lg_f_high = 19
    lg_f_low = 9
    
    if isnothing.(tau_max_override)
        tau_max = 1e8
    else
        tau_max = tau_max_override
    end
    
    if input_data == "Doddy"
        data = open(readdlm, "BH_data/Doddy_full.dat")
        use_input_table = true
    elseif input_data == "Masha"
        data = open(readdlm, "BH_data/Masha_Vals.dat")
        use_input_table = true
    else
        print("No Data! \n\n")
    end
    
    Fname = "LowMassRegion_TauMax_"*string(round(tau_max, sigdigits=2))
    Fname *= "_alpha_maxmin_"*string(round(alpha_max_cut, sigdigits=2))*"_"*string(round(alpha_min_cut, sigdigits=2))
    Fname *= "_"*input_data*Ftag
else
    ###### High Mass Region
    lg_m_low = -19.3
    lg_f_high = 19
    lg_f_low = 13
    lg_m_high = -16.0
    
    if isnothing.(tau_max_override)
        tau_max = 1e10
    else
        tau_max = tau_max_override
    end
    
    if input_data == "Doddy"
        lg_m_high = -15.5
        data = open(readdlm, "BH_data/Doddy_SMBH.dat")
        use_input_table = true
    elseif input_data == "Masha"
        lg_m_high = -15.5
        data = open(readdlm, "BH_data/Masha_SMBH.dat")
        use_input_table = true
    elseif input_data == "Andy"
        lg_m_high = -18.0
        # data = open(readdlm, "BH_data/TDE_test.dat")
        data = open(readdlm, "BH_data/TDE_twoBest.dat")
        use_input_table = false
        
        N_psi = 200  # Number of angles I throw stars in at (linearly spaced between 0 and pi/2).
        N_spin = 300  # Number of absolute values of spins linearly spaced between 0 and 0.9999
        nameMassMatrix = "input_info/max_mass_matrix_$(N_psi)_$(N_spin).txt"
        if isfile(nameMassMatrix)
            max_mass_matrix = readdlm(nameMassMatrix)
        else
            _, max_mass_matrix = get_hills_masses(N_psi, N_spin)  # Need this matrix for tidal force likelihood
            writedlm(nameMassMatrix, max_mass_matrix)
        end
    else
        print("No Data! \n\n")
    end
    
    Fname = "HighMassRegion_TauMax_"*string(round(tau_max, sigdigits=2))
    Fname *= "_alpha_maxmin_"*string(round(alpha_max_cut, sigdigits=2))*"_"*string(round(alpha_min_cut, sigdigits=2))
    Fname *= "_"*input_data*Ftag
end
solve_322 = true




function prior(theta, lg_m_low, lg_m_high, lg_f_low, lg_f_high)
    log_m, log_f = theta
    
    if (lg_m_low < log_m < lg_m_high)&&(lg_f_low < log_f < lg_f_high)
        return 0.0
    end
    return -Inf
end

function log_probability(theta, data, lg_m_low, lg_m_high, lg_f_low, lg_f_high; tau_max=1e4, alpha_max_cut=0.2, use_input_table=true, solve_322=true, impose_low_cut=1e-100, max_mass_matrix=nothing)

    lp = prior(theta, lg_m_low, lg_m_high, lg_f_low, lg_f_high)
    if !isfinite.(lp)
        return -Inf
    end

    return lp + log_likelihood(theta, data, tau_max=tau_max, alpha_max_cut=alpha_max_cut, use_input_table=use_input_table, solve_322=solve_322, impose_low_cut=impose_low_cut, max_mass_matrix=max_mass_matrix)
end

function log_likelihood(theta, data; tau_max=1e4, alpha_max_cut=0.2, use_input_table=true, solve_322=true, impose_low_cut=1e-100, max_mass_matrix=nothing)
    
    log_m, log_f = theta
    sum_loglike = 0.0
    
    if use_input_table
        MassBH_c = data[:, 1]
        MassBH_errU = data[:, 2]
        MassBH_errD = data[:, 3]
        SpinBH_c = data[:, 4]
        SpinBH_errU = data[:, 5]
        SpinBH_errD = data[:, 6]
        
        num_data = length(MassBH_c)
        
        for i in 1:num_data
        
            val_found = false
            p_or_neg = Int(round(rand()))
            MassBH = nothing
            d_mass = nothing
            if p_or_neg == 0
                d_mass = Normal(MassBH_c[i], MassBH_errU[i])
            else
                d_mass = Normal(MassBH_c[i], MassBH_errD[i])
            end
            
            # print("looking for mass \n")
            while !val_found
                MassBH = rand(d_mass,1)[1]
                if (p_or_neg == 0) && (MassBH >= MassBH_c[i])
                    val_found = true
                elseif (p_or_neg == 1) && (MassBH <= MassBH_c[i]) && (MassBH > 0)
                    val_found = true
                end
            end
            
            # print("looking for spin \n")
            SpinBH = nothing
            p_or_neg = Int(round(rand()))
            SpinBH = nothing
            val_found = false
            if p_or_neg == 0
                SpinBH = rand(Uniform(SpinBH_c[i], 0.998))
            else
                d_spin = Normal(SpinBH_c[i], SpinBH_errD[i])
                while !val_found
                    SpinBH = rand(d_spin,1)[1]
                    if (SpinBH <= SpinBH_c[i]) && (SpinBH >= 0.0)
                        val_found = true
                    end
                end
            end
                

            if SpinBH > 0.998
                SpinBH = 0.998
            end
            if SpinBH < 0.0
                SpinBH = 0.0
            end
            print(MassBH, "\t", SpinBH, "\t", 10 .^log_m, "\t", 10 .^log_f, "\n")
            final_spin = super_rad_check(MassBH, SpinBH, 10 .^log_m, 10 .^log_f, tau_max=tau_max, alpha_max_cut=alpha_max_cut, debug=false, solve_322=solve_322, impose_low_cut=alpha_min_cut)
            # print(final_spin, "\n\n")
            if final_spin > SpinBH_c[i]
                sum_loglike += -0.5 * (SpinBH_c[i] - final_spin).^2 / SpinBH_errU[i].^2
            else
                sum_loglike += -0.5 * (SpinBH_c[i] - final_spin).^2 / SpinBH_errD[i].^2
            end
        end
    else
        # assume data format: [Indx DataEntry Label], Label=0 [velocity disp] Label=1 Galaxy Mass
        # Indx = data[:, 1]
        # DataEntry = data[:, 2]
        # DataLabel = data[:, 3]
        ##### LOG 10 Values ########
        MassBH_c = data[:, 1]
        MassBH_errU = data[:, 2]
        MassBH_errD = data[:, 3]
        
    
        num_data = length(MassBH_c)
        
        for i in 1:num_data
            val_found = false
            p_or_neg = Int(round(rand()))
            MassBH = nothing
            d_mass = nothing
            if p_or_neg == 0
                d_mass = Normal(MassBH_c[i], MassBH_errU[i])
            else
                d_mass = Normal(MassBH_c[i], MassBH_errD[i])
            end
            
            while !val_found
                MassBH = 10 .^rand(d_mass,1)[1]
                if (p_or_neg == 0) && (MassBH >= 10 .^MassBH_c[i])
                    val_found = true
                elseif (p_or_neg == 1) && (MassBH <= 10 .^MassBH_c[i]) && (MassBH > 0)
                    val_found = true
                end
            end
            
            prior_spins(a) = ones(size(a))  # agnostic spin prior.
            min_spin = one_d_spin_fixed_mass(MassBH .* Ms, prior_spins, max_mass_matrix; return_all=false)
                                
            spinBH_sample = rand() .* (0.998 .- min_spin) .+ min_spin
            
            final_spin = super_rad_check(MassBH, spinBH_sample, 10 .^log_m, 10 .^log_f, tau_max=tau_max, alpha_max_cut=alpha_max_cut, debug=false, solve_322=solve_322, impose_low_cut=alpha_min_cut)
            loglike = tde_like(MassBH, final_spin, max_mass_matrix; plot=false)
            print(MassBH, "\t", 10 .^log_m, "\t", 10 .^log_f, "\t", spinBH_sample, "\t", final_spin, "\t", loglike, "\n")
            sum_loglike += loglike
            
        end
        
    
    end
    return sum_loglike
end

function spin_prior(x)
    return 1.0
end

function mass_prior(x, log10_M, sig)
    return exp.(-(x .- log10_M).^2 ./ (2 .* sig.^2)) ./ (sqrt.(2 .* pi) .* sig)
end

function BH_mass_from_obs(val, label)
    if label == 1
        MBH = 10 .^ (7.43 .+ 1.61 .+ log10.(val ./ 3e10))
        return MBH
    else
        MBH = 10 .^ (7.87 .+ 4.384 .* log10.(val ./ 160.0))
        return MBH
    end
end

function initialize_walkers(numwalkers, lg_m_low, lg_m_high, lg_f_low, lg_f_high)
    len_mass = abs.(lg_m_high - lg_m_low )
    len_fa = abs.(lg_f_high .- lg_f_low)
    
    x0 = rand(2, numwalkers)
    x0[1, :] .*= len_mass
    x0[1, :] .+= lg_m_low
    x0[2, :] .*= len_fa
    x0[2, :] .+= lg_f_low
    return x0
end

function mcmc_func_minimize(data, Fname; lg_m_low=-20, lg_m_high=-18, lg_f_high=19, lg_f_low=18, tau_max=1e4, alpha_max_cut=0.2, alpha_min_cut=1e-100, use_input_table=true, solve_322=true, numwalkers=10, thinning=1, numsamples_perwalker=2000, burnin=500, max_mass_matrix=nothing)
    numdims = 2


    function llhood(x)
        return log_probability(x, data, lg_m_low, lg_m_high, lg_f_low, lg_f_high, tau_max=tau_max, alpha_max_cut=alpha_max_cut, use_input_table=use_input_table, solve_322=solve_322, impose_low_cut=alpha_min_cut, max_mass_matrix=max_mass_matrix)
    end

    x0 = initialize_walkers(numwalkers, lg_m_low, lg_m_high, lg_f_low, lg_f_high)
    print("Starting burn-in...\n")
    chain, llhoodvals = AffineInvariantMCMC.sample(llhood, numwalkers, x0, burnin, 1)
    print("Starting main run...\n")
    chain, llhoodvals = AffineInvariantMCMC.sample(llhood, numwalkers, chain[:, :, end], numsamples_perwalker, thinning)
    print("Finished main run...\n")
    flatchain, flatllhoodvals = AffineInvariantMCMC.flattenmcmcarray(chain, llhoodvals)

    writedlm("output_mcmc/"*Fname*"_mcmc.dat", flatchain')
    writedlm("output_mcmc/"*Fname*"_likevals.dat", flatllhoodvals')
    
end
        
mcmc_func_minimize(data, Fname, lg_m_low=lg_m_low, lg_m_high=lg_m_high, lg_f_high=lg_f_high, lg_f_low=lg_f_low, tau_max=tau_max, alpha_max_cut=alpha_max_cut, alpha_min_cut=alpha_min_cut, use_input_table=use_input_table, solve_322=solve_322, numwalkers=numwalkers, thinning=thinning, numsamples_perwalker=numsamples_perwalker, burnin=burnin, max_mass_matrix=max_mass_matrix)
