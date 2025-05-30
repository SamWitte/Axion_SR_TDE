using ArgParse
include("stat_analysis.jl")


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
    
        "--alpha_max_cut"
            arg_type = Float64
            default = 1e3

        "--alpha_min_cut"
            arg_type = Float64
            default = 0.001

        "--tau_max_override"
            arg_type = Float64
            default = 5.0e6
            
        "--stop_on_a"
            arg_type = Float64
            default = 0.0
            
        "--eq_threshold"
            arg_type = Float64
            default = 1e-100
        
        "--abstol"
            arg_type = Float64
            default = 1e-30
            
        "--Ftag"
            arg_type = String
            default = "_"

        "--one_BH"
            arg_type = String
            default = "" # Cyg_X1, M33X7, nothing
            
        "--LowMass"
            arg_type = Bool
            default = true

        "--non_rel"
            arg_type = Bool
            default = false
                        
        "--numwalkers"
            arg_type = Int
            default = 10

        "--thinning"
            arg_type = Int
            default = 1
            
        "--numsamples_perwalker"
            arg_type = Int
            default = 10000
            
        "--burnin"
            arg_type = Int
            default = 1000
            
        "--Nmax"
            arg_type = Int
            default = 3
            
        "--cheby"
            arg_type = Bool
            default = false

    end

    return parse_args(s)
end

parsed_args = parse_commandline()

alpha_max_cut = parsed_args["alpha_max_cut"];
alpha_min_cut = parsed_args["alpha_min_cut"];
tau_max_override = parsed_args["tau_max_override"];
Ftag = parsed_args["Ftag"];
one_BH = parsed_args["one_BH"];
if one_BH == ""
    one_BH = nothing
end
LowMass = parsed_args["LowMass"];
numwalkers = parsed_args["numwalkers"];
thinning = parsed_args["thinning"];
numsamples_perwalker = parsed_args["numsamples_perwalker"];
burnin = parsed_args["burnin"];
non_rel = parsed_args["non_rel"]
stop_on_a = parsed_args["stop_on_a"]
eq_threshold = parsed_args["eq_threshold"]
abstol = parsed_args["abstol"]

Nmax = parsed_args["Nmax"]
cheby = parsed_args["cheby"]

print("Deets...\n\n")
println("Nmax: ", Nmax)
println("Use Cheby (if not, Leaver): ", cheby)
print("Alpha cuts \t", alpha_min_cut, "\t", alpha_max_cut, "\n")
print("Tau max \t ", tau_max_override, "\n")
print("Ftag \t", Ftag, "\n")
print("One BH? \t", one_BH, "\n")
print("Low mass BHs? \t", LowMass, "\n")
print("Non Rel? \t", non_rel, "\n")
print("Nwalkers, Nsamples, Nburn \t", [numwalkers numsamples_perwalker burnin], "\n\n\n")

# Random.seed!(2084339465932781399)

max_mass_matrix=nothing

if LowMass
    ###### Low Mass Region
    lg_m_low = -13
    lg_m_high = log10(3e-11)
    lg_f_high = 19
    lg_f_low = 9.5

#    lg_m_low = log10(3e-13)
#    lg_m_high = log10(2e-12)
#    lg_f_high = log10(1e19)
#    lg_f_low = log10(1e15)
    
    if isnothing.(tau_max_override)
        tau_max = 1e8
    else
        tau_max = tau_max_override
    end

    if !isnothing(one_BH)
        data = open(readdlm, "BH_data/"*one_BH*".dat")[2:end, :]
    else
        data = open(readdlm, "BH_data/Masha_Vals.dat")[2:end, :]
    end
    use_input_table = true
    
    Fname = "LowMassRegion_TauMax_"*string(round(tau_max, sigdigits=2))
    Fname *= "_alpha_maxmin_"*string(round(alpha_max_cut, sigdigits=2))*"_"*string(round(alpha_min_cut, sigdigits=2))

    if !isnothing(one_BH)
        Fname *= "_"*one_BH
    end
    Fname *= "_Nsamps_$(numsamples_perwalker)_"
    
    Fname *= "_Nmax_$(Nmax)_"
    if cheby
        Fname *= "_cheby_"
    end
 
    if non_rel
        Fname *= "_NonRel_"
    else
        Fname *= "_FullRel_"
    end
    
else
    ###### High Mass Region
    #### DONT RUN
    
    lg_m_low = -21.0
    lg_f_high = 19
    lg_f_low = 13
    lg_m_high = -16.0
    
    if isnothing.(tau_max_override)
        tau_max = 1e10
    else
        tau_max = tau_max_override
    end
    
    
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
    
    
    Fname = "HighMassRegion_TauMax_"*string(round(tau_max, sigdigits=2))
    Fname *= "_alpha_maxmin_"*string(round(alpha_max_cut, sigdigits=2))*"_"*string(round(alpha_min_cut, sigdigits=2))
    
    Fname *= "_Nmax_$(Nmax)_"
    if cheby
        Fname *= "_cheby_"
    end
    
    Fname *= "_Nsamps_$(numsamples_perwalker)_"
    
end


time0=Dates.now()
@inbounds @fastmath mcmc_func_minimize(data, Fname, lg_m_low=lg_m_low, lg_m_high=lg_m_high, lg_f_high=lg_f_high, lg_f_low=lg_f_low, tau_max=tau_max, alpha_max_cut=alpha_max_cut, alpha_min_cut=alpha_min_cut, use_input_table=use_input_table, numwalkers=numwalkers, thinning=thinning, numsamples_perwalker=numsamples_perwalker, burnin=burnin, max_mass_matrix=max_mass_matrix, stop_on_a=stop_on_a, eq_threshold=eq_threshold, abstol=abstol, non_rel=non_rel, Nmax=Nmax, cheby=cheby)

time1=Dates.now()
print("\n\n Run time: ", time1-time0, "\n")

