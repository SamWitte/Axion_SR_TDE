using SpecialFunctions
using Random
using OrdinaryDiffEq
using Statistics
using Distributions
using ForwardDiff: gradient, derivative, Dual, Partials, hessian
using DelimitedFiles
using NLsolve
using DifferentialEquations
using Interpolations
using HypergeometricFunctions
using SpinWeightedSpheroidalHarmonics
# using MCIntegration
# using Dates
# using LinearAlgebra
using Optim
include("Constants.jl")


function ergL(n, l, m, massB, MBH, a; full=true)
    # Key to next level
    alph = GNew * MBH * massB
    if full
        if l > 0
            return massB .* (1.0 .- alph.^2 ./ (2 .* n.^2) - alph.^4 ./ (8 * n.^4) + alph.^4 ./ n^4 .* (2 * l - 3 * n + 1) ./ (l + 0.5) + 2 * a * m * alph.^5 ./ n^3 ./ (l * (l + 0.5) * (l+1)))
        else
            return massB .* (1.0 .- alph.^2 ./ (2 .* n.^2) - alph.^4 ./ (8 * n.^4) + alph.^4 ./ n^4 .* (2 * l - 3 * n + 1) ./ (l + 0.5))
        end
    else
        return massB .* (1.0 .- alph.^2 ./ (2 .* n.^2) - alph.^4 ./ (8 * n.^4))
    end
end


function sr_rates(n, l, m, massB, MBH, aBH; impose_low_cut=0.001, solve_322=true)
    if (n==3)&&(l==2)&&(m==1)&&(solve_322==false)
        return 0.0
    end
    

    alph = GNew .* MBH .* massB
    
#    if (alph ./ l < impose_low_cut)&&(MBH < 1e2)
#        # We expect binaries to disrupt.
#        return 0.0
#    end
    
    
    rP = nothing
    if aBH .> maxSpin
        rP = 1.0 .+ sqrt.(1 - maxSpin .^2)
        aBH = maxSpin
    elseif aBH .< 0.0
        rP = 2.0
        aBH = 0.0
    else
        rP = 1.0 .+ sqrt.(1 - aBH.^2)
    end
    rP *= (GNew .* MBH)
    OmegaH = aBH ./ (2 .* (GNew .* MBH) .* (1 .+ sqrt.(1 .- aBH.^2)))
    Anl = 2 .^(4 .* l .+ 1) .* factorial(big(Int(l .+ n))) ./ (n.^(2 .* l .+ 4) .* factorial(big(n .- l .- 1)))
    Anl *= (factorial(Int(l)) ./ (factorial(Int(2 .* l)) .* factorial(Int(2 .* l .+ 1)))).^2
    Chilm = 1.0
    # erg = ergL(n, l, m, massB, MBH, aBH)
    for k in 1:Int(l)
        Chilm *= (k.^2 .* (1.0 .- aBH.^2) .+ (aBH * m .- 2 .* (rP ./ (GNew .* MBH)) .* alph).^2)
        # Chilm *= (k.^2 .* (1.0 .- aBH.^2) .+ (aBH * m .- 2 .* (rP ./ (GNew .* MBH)) .* alph).^2)
    end
    Gamma_nlm = 2 * massB .* rP .* (m .* OmegaH .- ergL(n, l, m, massB, MBH, aBH)) .* alph.^(4 .* l + 4) .* Anl .* Chilm
    

    if Gamma_nlm > 0.0
        return Gamma_nlm
    elseif (Gamma_nlm < 0.0)&&(m .* OmegaH .> ergL(n, l, m, massB, MBH, aBH))
        return 0.0
    else
        return Gamma_nlm
    end
end

function radial_bound_NR(n, l, m, mu, M, r; physU=false)
    # r in units r -> r / (GM)
    alph = GNew * M * mu
    
    if !physU
        a0 = 1 / alph.^2
    else
        a0 = 1 / (mu * alph)
    end
    rF = sqrt.((2 / (n * a0)).^3 * factorial(big(n - l - 1)) / (2 * n * factorial(big(n + l)))) .* exp.( - r ./ (n .* a0)) .* (2 .* r ./ (n .* a0)).^l .* generalized_laguerre(n - l - 1, 2 * l + 1, 2 .* r ./ (n .* a0)) # I've thrown through a factor of (GNew * M)^{3/2}
    return rF
end

function radial_inf_NR(k, l, mu, M, r)
    # r assumed to be in [G M]
    # k unitless
    alph = GNew * M * mu
    a0 = 1 / alph.^2
    rF = zeros(Complex, length(r))
    for i in 1:length(r)
        rF[i] = 2 .* k .* exp.(pi ./ (2 .* k .* a0)) .* abs.(gamma(l + 1 - im ./ (k .* a0))) ./ factorial(2 * l + 1) .* (2 .* k .* r[i]).^l .* exp.(-im * k * r[i]) .* pFq((im ./ (k .* a0) + l + 1,), (2 * l + 2, ), 2 * im .* k .* r[i])
    end
    return rF # unitless
end


function freq_shifts(mu, M, a, n1, l1, m1, n2, l2, m2;  rpts=500, rmaxT=100, Nang=100000, epsil_2=1.0, Npts_Bnd=500)
    rp = 1 + sqrt.(1 - a.^2)
    alph = mu * GNew * M
    

    erg_1G = ergL(n1, l1, m1, mu, M, a)
    q = abs.(real(- sqrt.(alph.^2 .- (erg_1G .* GNew .* M).^2)))
    rmax = rmaxT ./ q
    rlist = 10 .^(range(log10.(rp), log10.(rmax), rpts))
    
    rl, r1, erg_1 = solve_radial(mu, M, a, n1, l1, m1; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true)
    itp = LinearInterpolation(log10.(rl), log10.(r1), extrapolation_bc=Line())
    rf_1 = 10 .^itp(log10.(rlist))
    Z1 = spheroidals(l1, m1, a, erg_1 ./ (GNew .* M))
    
    # rf_1 = radial_bound_NR(n1, l1, m1, mu, M, rlist)

    
    if (n2 == n1)&&(l2==l1)&&(m2==m1)
        rf_2 = rf_1
        erg_2 = erg_1
        Z2 = Z1
        mult_fac = 1.0 ./ 2.0
    else
        rl, r2, erg_2 = solve_radial(mu, M, a, n2, l2, m2; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true)
        itp = LinearInterpolation(log10.(rl), log10.(r2), extrapolation_bc=Line())
        rf_2 = 10 .^itp(log10.(rlist))
        Z2 = spheroidals(l2, m2, a, erg_2 ./ (GNew .* M))
        mult_fac = 1.0
    end
    
    function func_ang(x)
        return real(Z1.(x[1], x[2]) .* Z2.(x[1], x[2]) .* conj(Z1.(x[1], x[2])) .* conj(Z2.(x[1], x[2])))
    end
    trapz(y,x) = @views sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
    
    thetaV = acos.(1.0 .- 2.0 .* rand(Nang))
    phiV =  2.0 * pi .* rand(Nang)
    CG = 0.0
    for i in 1:Nang
        CG += func_ang([thetaV[i], phiV[i]])
    end
    CG *= 4 * pi / Nang
    
    r_integrd = rf_1 .* rf_2 .* conj(rf_1) .* conj(rf_2) .* rlist.^2
    radial_int = trapz(r_integrd, rlist) ./ (GNew * M).^3
    
    delt_om = -1.0 ./ (4 * mu^2) .* radial_int .* CG .* mult_fac
    
    return Float64(real(alph.^2 .* delt_om .* epsil_2)) # \delta Omega [Needs to be multiplied by (M_pl / f_a)^2], do dF / (alpha^5 * mu) to get coefficient hitting epsilon_x

end

function s_rate_bnd(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; kpts=10, rpts=2000, rmaxT=100, inf_nr=true, Nang=100000, Npts_Bnd=1000, debug=false, include_cont=true, Ntot_safe=5000, sve_for_test=false, bnd_thresh=1e-3, use_analytic=false, eps_r=1e-10)
    
    rp = 1 + sqrt.(1 - a.^2)
    alph = mu * GNew * M
    
    # approx energy estimates
    erg_1G = ergL(n1, l1, m1, mu, M, a)
    erg_2G = ergL(n2, l2, m2, mu, M, a)
    erg_3G = ergL(n3, l3, m3, mu, M, a)
    
    q1 = abs.(real(- sqrt.(alph.^2 .- (erg_1G .* GNew .* M).^2)))
    q2 = abs.(real(- sqrt.(alph.^2 .- (erg_2G .* GNew .* M).^2)))
    q3 = abs.(real(- sqrt.(alph.^2 .- (erg_3G .* GNew .* M).^2)))
    
    qmax = maximum([q1, q2, q3])
    qmin = minimum([q1, q2, q3])
    rmax = rmaxT ./ qmax
    rlist = 10 .^(range(log10.(rp .* (1.0 .+ eps_r)), log10.(rmax), rpts))
    
    if !use_analytic
        rl, r1, erg_1 = solve_radial(mu, M, a, n1, l1, m1; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe)
        itp = LinearInterpolation(log10.(rl), log10.(r1), extrapolation_bc=Line())
        rf_1 = 10 .^itp(log10.(rlist))
    else
        rf_1 = radial_bound_NR(n1, l1, m1, mu, M, rlist)
        erg_1 = erg_1G .* GNew * M
    end
    

    
    if (n2 == n1)&&(l2==l1)&&(m2==m1)
        rf_2 = rf_1
        erg_2 = erg_1
    else
        if !use_analytic
            rl, r2, erg_2 = solve_radial(mu, M, a, n2, l2, m2; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe)
            itp = LinearInterpolation(log10.(rl), log10.(r2), extrapolation_bc=Line())
            rf_2 = 10 .^itp(log10.(rlist))
        else
            rf_2 = radial_bound_NR(n2, l2, m2, mu, M, rlist)
            erg_2 = erg_2G .* GNew * M
        end
    end

    if !use_analytic
        rl, r3, erg_3 = solve_radial(mu, M, a, n3, l3, m3; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe)
        itp = LinearInterpolation(log10.(rl), log10.(r3), extrapolation_bc=Line())
        rf_3 = 10 .^itp(log10.(rlist))
    else
        rf_3 = radial_bound_NR(n3, l3, m3, mu, M, rlist)
        erg_3 = erg_3G .* GNew * M
    end
    

    erg_ind = erg_1 .+ erg_2 - erg_3
    k_ind_2 = mu.^2 .- (erg_ind ./ (GNew .* M)).^2
    
    # erg_ind = erg_1G .+ erg_2G - erg_3G
    # k_ind_2 = mu.^2 .- erg_ind.^2
    

    Z1 = spheroidals(l1, m1, a, erg_1 ./ (GNew .* M))
    Z2 = spheroidals(l2, m2, a, erg_2 ./ (GNew .* M))
    Z3 = spheroidals(l3, m3, a, erg_3 ./ (GNew .* M))
    
    
    trapz(y,x) = @views sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
        

    # compute bound contribution
    kmin = 0.02 .* alph^2
    kmax = 2.0 .* alph^2
    kk_list = 10 .^LinRange(log10.(kmin), log10.(kmax),  kpts)
    # print(kmin, "\t", kmax, "\n")
    ck_list = zeros(Complex, kpts)
 
    UnF = GNew .* M
    bound_c = 0.0
    contin_c = 0.0
    done_nmax = false
    n = 1
    
    function func_ang(x, Zf1, Zf2, Zf3, Zf4)
        return real(Zf1.(x[1], x[2]) .* Zf2.(x[1], x[2]) .* conj(Zf3.(x[1], x[2])) .* conj(Zf4.(x[1], x[2])))
    end
        
        
    while !done_nmax
        
        rl, r4, erg_4 = solve_radial(mu, M, a, n, 0, 0; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe)
        itp = LinearInterpolation(log10.(rl), log10.(r4), extrapolation_bc=Line())
        rf_4 = 10 .^itp(log10.(rlist))
        # rf_4 = radial_bound_NR(n, 0, 0, mu, M, rlist)
        
        # erg_4G = ergL(n, 0, 0, mu, M, a)
        Z4 = spheroidals(0, 0, a, erg_4 ./ (GNew .* M))
        
        thetaV = acos.(1.0 .- 2.0 .* rand(Nang))
        phiV = rand(Nang) .* 2*pi
        CG = 0.0
        for i in 1:Nang
            CG += func_ang([thetaV[i], phiV[i]], Z1, Z2, Z3, Z4)
        end
        CG *= 4*pi / Nang
       
       
        if sve_for_test&&(n==2)
            writedlm("test_store/check_real_1.dat", hcat(rlist, float(real(rf_1))))
            writedlm("test_store/check_imag_1.dat", hcat(rlist, float(imag(rf_1))))
            writedlm("test_store/check_real_2.dat", hcat(rlist, float(real(rf_2))))
            writedlm("test_store/check_imag_2.dat", hcat(rlist, float(imag(rf_2))))
            writedlm("test_store/check_real_3.dat", hcat(rlist, float(real(rf_3))))
            writedlm("test_store/check_imag_3.dat", hcat(rlist, float(imag(rf_3))))
            writedlm("test_store/check_real_4.dat", hcat(rlist, float(real(rf_4))))
            writedlm("test_store/check_imag_4.dat", hcat(rlist, float(imag(rf_4))))
        end
        
        r_integrd = (rf_1 ./ UnF^(3/2)) .* (rf_2 / UnF^(3/2)) .* conj(rf_3 / UnF^(3/2)) .* conj(rf_4 / UnF^(3/2)) .* (rlist * UnF).^2
        radial_int = trapz(r_integrd, rlist * UnF)
        
        
        kdiff_sq = (erg_4.^2 - erg_ind^2) ./ (GNew .* M).^2
        
        ff = CG .* radial_int ./ (2 .* mu).^(3 / 2) ./ (kdiff_sq)
        
        ### double factor 1/2
        if (n1==n2)&&(l1==l2)&&(m1==m2)
            ff /= 2
        end
        res_n = (rf_4[1] ./ UnF^(3/2)) .* ff
        if debug
            print(n, "\t", Float64(abs.(res_n)), "\n")
        end
        
        bound_c += res_n
        if abs.(res_n / bound_c) < bnd_thresh
            done_nmax = true
        end
        n += 1
    end
    


    if include_cont
        # compute continuous contribution
        for i in 1:kpts
            k = kk_list[i] ./ (GNew .* M) # physical units
            erg_New = sqrt.(k.^2 .+ mu.^2)
            Z4 = spheroidals(0, 0, a, erg_New)
            


            thetaV = acos.(1.0 .- 2.0 .* rand(Nang))
            phiV = rand(Nang) .* 2*pi

            CG = 0.0
            
            for i in 1:Nang
                CG += func_ang([thetaV[i], phiV[i]], Z1, Z2, Z3, Z4)
            end
            CG *= 4*pi / Nang
        
            if inf_nr
                out_goingR = radial_inf_NR(k .* GNew .* M, 0, mu, M, rlist)
            else
                rl, r4 = radial_inf(erg_New .* GNew .* M, mu, M, a, 0, 0; rmax_val=rmax, rpts=rpts)
                # itp = LinearInterpolation(log10.(rl[3:end]), log10.(r4[3:end]), extrapolation_bc=Line())
                itp = LinearInterpolation(log10.(rl[3:end]), r4[3:end], extrapolation_bc=Line())
                # out_goingR = 10 .^itp(log10.(rlist))
                out_goingR = itp(log10.(rlist))
            end
            
            r_integrd = (rf_1 / UnF^(3/2)) .* (rf_2 / UnF^(3/2)) .* conj(rf_3 / UnF^(3/2)) .* conj(out_goingR ./ UnF) .* (rlist * UnF).^2
            radial_int = trapz(r_integrd, rlist * UnF)
            
            
            # print("radial int \t ", k, "\t", radial_int, "\n\n")
            
            
            ff = CG .* radial_int ./ (2 .* mu).^(3 / 2) ./ (k_ind_2 .+ k.^2) ./ (2 * pi)
            
            ### double factor 1/2
            if (n1==n2)&&(l1==l2)&&(m1==m2)
                ff /= 2
            end
            
            if !isnan.(out_goingR[1] .* ff ./ UnF)
                ck_list[i] = out_goingR[1] .* ff ./ UnF
            else
                print("getting nan... \t", i, "\t", out_goingR[1:3], "\t", ff, "\n")
            end
            if debug
                print("cont \t", i, "\t", kk_list[i], "\t", Float64(abs.(ck_list[i])), "\t", ff, "\n")
            end
        end
        contin_c += trapz(ck_list, kk_list / UnF)
    end
    
    psi_1 = bound_c .+ contin_c
    # print("bound contribution:  ", Float64(abs.(bound_c)), "\t contin: ", Float64(abs.(contin_c)), "\n")
        
    
    lam = (mu ./ (M_pl .* 1e9))^2
    rate_out = 4 .* alph.^2 .* (1 .+ sqrt.(1 - a.^2)) .* Float64(real(psi_1 .* conj(psi_1))) .* lam^2
    # rate_out = 2 * alph.^2 .* rp * Float64(real(psi_1 .* conj(psi_1))) .* lam^2
    
    return rate_out ./ mu^2 .* (GNew * M^2 * M_to_eV)^2 # unitless [gamma / mu]
end

function s_rate_inf(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3, lF_min; rpts=4000, rmaxT=90,  sve_for_test=false, inf_nr=false, Npts_Bnd=1000, Nang=300000, debug=false, Ntot_safe=2000, xtol=1e-2, ftol=1e-2, iter=20)
    # lF_min is l relative to minimal value needed for non-zero m
    # returns scattering re-normalized \gamma (basically just radial integral ratio)
    

    trapz(y,x) = @views sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
    
    rp = (1.0 .+ sqrt.(1 - a.^2))
    alph = mu .* GNew .* M
    
    # rough estimate energies
    erg_1 = ergL(n1, l1, m1, mu, M, a)
    erg_2 = ergL(n2, l2, m2, mu, M, a)
    erg_3 = ergL(n3, l3, m3, mu, M, a)
    erg_New = erg_1 .+ erg_2 - erg_3
    # k = sqrt.(erg_New.^2 .- mu.^2)
    if debug
        print("erg new \t", erg_New .* GNew * M, "\n")
    end
    k = sqrt.((erg_New .* GNew * M).^2 .- alph.^2)
    mF = (m1 .+ m2 - m3)
    lF = abs.(mF) + lF_min
    
    
    
    q1 = abs.(real(- sqrt.(alph.^2 .- (GNew .* M .* erg_1).^2)))
    q2 = abs.(real(- sqrt.(alph.^2 .- (GNew .* M .* erg_2).^2)))
    q3 = abs.(real(- sqrt.(alph.^2 .- (GNew .* M .* erg_3).^2)))
    q = maximum([q1, q2, q3])
    qcheck = 6 ./ (k .* GNew .* M)
    rmax = rmaxT ./ q
    
    rlist = 10 .^(range(log10.(rp), log10.(rmax), rpts))
    # print("New erg \t", erg_New .* GNew .* M, "\n")


    Z1 = spheroidals(l1, m1, a, erg_1)
    Z2 = spheroidals(l2, m2, a, erg_2)
    Z3 = spheroidals(l3, m3, a, erg_3)
    Z4 = spheroidals(lF, mF, a, erg_New)
   
    thetaV = acos.(1.0 .- 2.0 .* rand(Nang))
    phiV = rand(Nang) .* 2*pi

    function func_ang(x)
        return real(Z1.(x[1], x[2]) .* Z2.(x[1], x[2]) .* conj(Z3.(x[1], x[2])) .* conj(Z4.(x[1], x[2]))) #
    end
    CG = 0.0
    ang_fac2 = 0.0
    
    for i in 1:Nang
        CG += func_ang([thetaV[i], phiV[i]])
    end
    CG *= 4*pi / Nang
   
    
    rl, r1 = solve_radial(mu, M, a, n1, l1, m1; rpts=Npts_Bnd, rmaxT=rmaxT, Ntot_safe=Ntot_safe)
    itp = LinearInterpolation(log10.(rl), log10.(r1), extrapolation_bc=Line())
    rf_1 = 10 .^itp(log10.(rlist))
    
    if sve_for_test
        writedlm("test_store/check_real_0.dat", hcat(rl, float(real(r1))))
        writedlm("test_store/check_imag_0.dat", hcat(rl, float(imag(r1))))
    end
    
    if (n2 == n1)&&(l2==l1)&&(m2==m1)
        rf_2 = rf_1
    else
        rl, r2 = solve_radial(mu, M, a, n2, l2, m2; rpts=Npts_Bnd, rmaxT=rmaxT, Ntot_safe=Ntot_safe)
        itp = LinearInterpolation(log10.(rl), log10.(r2), extrapolation_bc=Line())
        rf_2 = 10 .^itp(log10.(rlist))
    end

    
    rl, r3 = solve_radial(mu, M, a, n3, l3, m3; rpts=Npts_Bnd, rmaxT=rmaxT, Ntot_safe=Ntot_safe)
    itp = LinearInterpolation(log10.(rl), log10.(r3), extrapolation_bc=Line())
    rf_3 = 10 .^itp(log10.(rlist))
    
    
    if inf_nr
        out_goingR = radial_inf_NR(k, lF, mu, M, rlist)
    else
        rl, r4 = radial_inf(erg_New .* GNew .* M, mu, M, a, lF, mF; rmax_val=rmax, rpts=rpts, debug=debug, xtol=xtol, ftol=ftol, iter=iter)
        # itp = LinearInterpolation(log10.(rl[3:end]), log10.(r4[3:end]), extrapolation_bc=Line())
        # out_goingR = 10 .^itp(log10.(rlist))
        itp = LinearInterpolation(log10.(rl[3:end]), r4[3:end], extrapolation_bc=Line())
        out_goingR = itp(log10.(rlist))
    end
    
    if sve_for_test
        writedlm("test_store/check_real_1.dat", hcat(rlist, float(real(rf_1))))
        writedlm("test_store/check_imag_1.dat", hcat(rlist, float(imag(rf_1))))
        writedlm("test_store/check_real_2.dat", hcat(rlist, float(real(rf_2))))
        writedlm("test_store/check_imag_2.dat", hcat(rlist, float(imag(rf_2))))
        writedlm("test_store/check_real_3.dat", hcat(rlist, float(real(rf_3))))
        writedlm("test_store/check_imag_3.dat", hcat(rlist, float(imag(rf_3))))
        writedlm("test_store/check_real_4.dat", hcat(rlist, float(real(out_goingR))))
        writedlm("test_store/check_imag_4.dat", hcat(rlist, float(imag(out_goingR))))
    end
    
    unitsF = GNew * M
    r_integrd = (rf_1 ./ unitsF^(3/2))  .* (rf_2 ./ unitsF^(3/2)) .* conj(rf_3 / unitsF^(3/2)) .* conj(out_goingR / unitsF) .* (rlist * unitsF).^2
    radial_int = trapz(r_integrd, rlist * unitsF)
    
    # Rbound [1/GM]^{3/2}, R_c [eV], r [GM]
    ff = CG .* 4 * pi * (- im).^lF ./ (2 .* k ./ unitsF) .* radial_int ./ (2 .* mu).^(3 / 2)
    
    ### double factor 1/2
    if (n1==n2)&&(l1==l2)&&(m1==m2)
        ff /= 2
    end
    
    lam = (mu ./ (M_pl .* 1e9))^2
    rate_out = Float64(real(2 .* erg_New * (k / unitsF) ./ (4 * pi)^2 .* ff .* conj(ff) .* lam.^2))
    
    return rate_out ./ mu^2 .* (GNew * M^2 * M_to_eV)^2 # unitless [gamma / mu]
    # Return: rate  (evaluated at f_a = Mpl)

    
end

function solve_radial(mu, M, a, n, l, m; rpts=1000, rmaxT=50, debug=false, iter=500, xtol=1e-20, ftol=1e-90, sve=false, fnm="test_store/WF_", return_erg=false, Ntot_safe=5000, eps_r=1e-10)
    ### dolan 2007
    # everything normalized
    alph = GNew .* M .* mu
    
    if m == 1
        Ntot = 1500
    elseif (m==2)
        Ntot = 4250
    elseif (m==3)
        Ntot = 5000
    else
        Ntot = Ntot_safe
    end
    
    if a > 0.95
        Ntot *= 2.5
        Ntot = Int(Ntot)
    end
    if Ntot_safe > Ntot
        Ntot = Ntot_safe
    end
    
    wR, wI = find_im_part(mu, M, a, n, l, m; debug=debug, iter=iter, xtol=xtol, ftol=ftol, return_both=true, for_s_rates=true)
    erg = wR .+ im .* wI
    if debug
        print("erg \t", erg, "\n")
    end
    
    q = - sqrt.(alph.^2 .- erg.^2)
    
    # rmax = rmaxT * 1.0 ./ abs.(real(q))
    rmax = Float64.(100 * n ./ alph.^2)
    # print("Check \t", rmax, "\t", rmaxT * 1.0 ./ abs.(real(q)), "\n")

    
    
    rm = (1.0 .- sqrt.(1 - a.^2))
    rp = (1.0 .+ sqrt.(1 - a.^2))
    OmH = a  ./ (2 .* rp)
    

    
    
    sigm = 2 .* rp .* (erg .- OmH .* m) ./ (rp .- rm)
    chi = (alph.^2 .- 2 .* erg.^2) ./ q
    b = sqrt.(1 - a.^2)
    
    gam = im * a * sqrt.(erg.^2 .- alph.^2)
    cc2 = a.^2 .* (erg.^2 .- alph.^2)
    LLM = l * (l + 1)
    # LLM += 2 .* cc2 .* (m.^2 .- l .* (l + 1) + 0.5) ./ ((2 * l - 1) .* (2 * l + 3)) # alternative describ 2nd order
    LLM += (-1 + 2 * l * (l + 1) - 2 * m.^2) * gam.^2 ./ (-3 + 4 * l * (l + 1))
    LLM += ((l - m - 1 * (l - m) * (l + m) * (l + m - 1)) ./ ((-3 + 2 * l) * (2 * l - 1).^2) - (l + 1 - m) * (2 * l - m) * (l + m + 1) * (2 + l + m) ./ ((3 + 2 * l).^2 * (5 + 2 * l))) * gam.^4 ./ (2 * (1 + 2 * l))
    LLM += (4 * ((-1 + 4 * m^2) * (l * (1 + l) * (121 + l * (1 + l) * (213 + 8 * l * (1 + l) * (-37 + 10 * l * (1 + l)))) - 2 * l * (1 + l) * (-137 + 56 * l * (1 + l) * (3 + 2 * l * (1 + l))) * m^2 + (705 + 8 * l * (1 + l) * (125 + 18 * l * (1 + l))) * m^4 - 15 * (1 + 46 * m^2))) * gam^6) / ((-5 + 2 * l) * (-3 + 2 * l) * (5 + 2 * l) * (7 + 2 * l) * (-3 + 4 * l * (1 + l))^5)
    
    c0 = 1.0 .- 2.0 * im * erg - 2 * im ./ b .* (erg .- a .* m ./ 2.0)
    c1 = -4.0 .+ 4 * im * (erg - im * q * (1.0 + b)) + 4 * im / b * (erg - a * m / 2.0) .- 2.0 * (erg.^2 .+ q.^2) ./ q
    c2 = 3.0 - 2 * im * erg - 2.0 * (q.^2 - erg.^2) ./ q - 2.0 * im / b * (erg - a * m ./ 2)
    c3 = 2.0 * im * (erg - im * q).^3 ./ q .+ 2 * (erg .- im * q).^2 .* b + q.^2 * a.^2 .+ 2 * im * q * a * m - LLM - 1 - (erg - im * q).^2 ./ q .+ 2 * q * b + 2 * im / b * ( (erg - im * q).^2 ./ q + 1.0) * (erg .- a * m / 2)
    c4 = (erg .- im * q).^4 ./ q.^2 .+ 2 * im * erg * (erg .- im * q).^2 ./ q .- 2 * im ./ b .* (erg .- im * q).^2 ./ q .* (erg .- a * m ./ 2)
    
    function alphaN(nn)
        return nn.^2 .+ (c0 .+ 1) * nn + c0
    end
    function betaN(nn)
        return -2 * nn.^2 .+ (c1 .+ 2) * nn + c3
    end
    function gammaN(nn)
        return nn.^2 .+ (c2 .- 3) * nn + c4
    end
    
    
    rlist = 10 .^(range(log10.(rp  .* (1.0 .+ eps_r)), log10.(rmax), rpts))
    if (wR == 0)&&(wI == 0)
        return rlist, zeros(length(rlist))
    end
    
    preF = (rlist .- rp).^(- im .* sigm) .* (rlist .- rm).^(im .* sigm .+ chi .- 1.0) .* exp.(q .* rlist)
    
    nfactor = ((rlist .- rp) ./ (rlist .- rm))
    Rout = zeros(Complex, length(rlist))
    bk = zeros(Complex, Ntot)
    bk[1] = 1.0
    bk[2] = -betaN(0) ./ alphaN(0) .* bk[1]
    
    for i in 2:(Ntot-1)
        bk[i+1] = (-betaN(i-1) .* bk[i] .- gammaN(i-1) .* bk[i-1]) ./ alphaN(i-1)
    end
    
    for i in 1:Ntot
        Rout .+= preF .* bk[i] .* nfactor.^(i-1)
    end
    
    # Normalize integration
    trapz(y,x) = @views sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
    
    nm = trapz(Rout .* conj(Rout) .* rlist.^2, rlist)
    
    # print(Float64.(real(Rout ./ sqrt.(nm))), "\n")
    if sve
        normOut = Rout .* conj(Rout) ./ nm
        writedlm(fnm*"_$(alph)_n_$(n)_l_$(l)_m_$(m)_.dat", hcat(float(real(rlist)), float(real(normOut))))
    else
        if return_erg
            return rlist, Rout ./ sqrt.(nm), erg # erg normalized by GM
        else
            return rlist, Rout ./ sqrt.(nm)
        end
    end
end

function radial_inf(erg, mu, M, a, l, m; rpts=1000, rmax_val=1e4, debug=false, iter=50, xtol=1e-120, ftol=1e-120, sve_for_test=false, fnm="test_store/test_radial", is_infin=true)
    # r, erg, unitless
    
    
    rp = 1.0 .+ sqrt.(1.0 .- a.^2)
    r = LinRange(rp .* 1.01, rmax_val .* rp, rpts)
    h = r[2] .- r[1]
    
    alph = mu .* GNew .* M
    k = sqrt.(Complex.(erg.^2 .- alph.^2))
    b = sqrt.(1 - a.^2)
    gam2 = - a * sqrt.(erg.^2 .- alph.^2)
    LLM = l * (l + 1)
    LLM += (-1 + 2 * l * (l + 1) - 2 * m.^2) * gam2 ./ (-3 + 4 * l * (l + 1))
    LLM += ((l - m - 1 * (l - m) * (l + m) * (l + m - 1)) ./ ((-3 + 2 * l) * (2 * l - 1).^2) - (l + 1 - m) * (2 * l - m) * (l + m + 1) * (2 + l + m) ./ ((3 + 2 * l).^2 * (5 + 2 * l))) * gam2.^2 ./ (2 * (1 + 2 * l))
    LLM += (4 * ((-1 + 4 * m^2) * (l * (1 + l) * (121 + l * (1 + l) * (213 + 8 * l * (1 + l) * (-37 + 10 * l * (1 + l)))) - 2 * l * (1 + l) * (-137 + 56 * l * (1 + l) * (3 + 2 * l * (1 + l))) * m^2 + (705 + 8 * l * (1 + l) * (125 + 18 * l * (1 + l))) * m^4 - 15 * (1 + 46 * m^2))) * gam2^3) / ((-5 + 2 * l) * (-3 + 2 * l) * (5 + 2 * l) * (7 + 2 * l) * (-3 + 4 * l * (1 + l))^5)
    
    
    delt = r.^2 .- 2 .* r .+ a.^2
    dr_delt = 2 .* r .- 2
    # rhs = Float64.((erg.^2 .* (r.^2 .+ a.^2).^2 .- 4 .* a .* m .* erg .* r .+ m.^2 .* a.^2) ./ delt .- (erg.^2 .* a.^2 .+ alph.^2 .* r.^2 .+ LLM))
    # rhs = real((erg.^2 .* (r.^2 .+ a.^2).^2 .- 4 .* a .* m .* erg .* r .+ m.^2 .* a.^2) ./ delt .- (erg.^2 .* a.^2 .+ alph.^2 .* r.^2 .+ LLM))
    rhs = ((erg.^2 .* (r.^2 .+ a.^2).^2 .- 4 .* a .* m .* erg .* r .+ m.^2 .* a.^2) ./ delt .- (erg.^2 .* a.^2 .+ alph.^2 .* r.^2 .+ LLM))
    
    # Create the second derivative matrix (Laplacian matrix), and first derivative
    D2 = zeros(rpts, rpts)

    # Fill the diagonal
    for i in 2:rpts-1
        D2[i, i-1] += 1.0 .* delt[i] ./ h^2
        D2[i, i] += -2.0 .* delt[i] ./ h^2
        D2[i, i+1] += 1.0 .* delt[i] ./ h^2
        
        D2[i, i+1] += 1.0 .* dr_delt[i] ./ (2 * h)
        D2[i, i-1] += -1.0 .* dr_delt[i] ./ (2 * h)
    end

    # BNDRY TERMS
    D2[1, 1] += -3.0 .* dr_delt[1] ./ (2 * h)
    D2[1, 2] += 4.0 .* dr_delt[1] ./ (2 * h)
    D2[1, 3] += -1.0 .* dr_delt[1] ./ (2 * h)
    D2[rpts, rpts] += 3.0 .* dr_delt[rpts] ./ (2 * h)
    D2[rpts, rpts - 1] += -4.0 .* dr_delt[rpts] ./ (2 * h)
    D2[rpts, rpts - 2] += 1.0 .* dr_delt[rpts] ./ (2 * h)
    
    D2[1, 1] += 2.0 .* delt[1] ./ h.^3
    D2[1, 2] += -5.0 .* delt[1] ./ h.^3
    D2[1, 3] += 4.0 .* delt[1] ./ h.^3
    D2[1, 4] += -1.0 .* delt[1] ./ h.^3
    D2[rpts, rpts] += 2.0 .* delt[rpts] ./ h.^3
    D2[rpts, rpts - 1] += -5.0 .* delt[rpts] ./ h.^3
    D2[rpts, rpts - 2] += 4.0 .* delt[rpts] ./ h.^3
    D2[rpts, rpts - 3] += -1.0 .* delt[rpts] ./ h.^3
    

    
    function wrapper!(F, x)
        # temp = (D2 * x) .+ rhs .* x
        temp = x.^(-1) .* (D2 * x) .+ rhs
        # F .= real(temp) .+ imag(temp)
        F .= abs.(temp)
        
        
        # real_r = x[1:Int(length(x)/2)]
        # imag_r = x[Int(length(x)/2 + 1):end]
        # r_tot = real_r .+ im .* imag_r
        # temp = (D2 * r_tot) .+ rhs .* r_tot
        # F[1:Int(length(x)/2)] .= real(temp) .* 1e60
        # F[Int(length(x)/2 + 1):end] .= imag(temp) .* 1e60
    end
    
    k = sqrt.(erg.^2 .- alph.^2)
    
    rGuess = radial_inf_NR(k, l, mu, M, r)
    
    # normalize solution!
    trapz(y,x) = @views sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
  
    r_in = real(rGuess)
    
    sol = nlsolve(wrapper!, r_in, show_trace=debug, autodiff = :forward, xtol=xtol, ftol=ftol, iterations=iter)
    full_out = sol.zero
    
   
    if is_infin
        ## need to normalize
        done_search=false
        found1=false
        found2=false
        cnt1 = 0
        cnt2 = 0
        
        indx = rpts
        while !done_search
            if sign(real(full_out[indx])) != sign(real(full_out[indx-1]))
                cnt1 += 1
            end
            if sign(real(r_in[indx])) != sign(real(r_in[indx-1]))
                cnt2 += 1
            end
            
            if cnt1 > 4
                found1 = true
            end
            if cnt2 > 4
                found2 = true
            end
                
            indx -= 1
            if found1&&found2
                done_search = true
            end
        end
        nm1 = maximum(abs.(full_out[indx:end]))
        nm2 = maximum(abs.(r_in[indx:end]))
        
        full_out .*= nm2 ./ nm1
        
        if sve_for_test
            writedlm(fnm*"_real.dat", hcat(r, real(full_out)))
            writedlm(fnm*"_imag.dat", hcat(r, imag(full_out)))
        end
        
        return r, full_out
    else
        
    end
end

function spheroidals(l, m, a, erg)
    # pass erg in normalized units
    Zlm = spin_weighted_spheroidal_harmonic(0, l, m, a .* erg)
    return Zlm
end

function find_im_part(mu, M, a, n, l, m; debug=false, Ntot_force=200, iter=10000, xtol=1e-20, ftol=1e-90, return_both=false, for_s_rates=false, QNM=false, QNM_E=1.0, erg_Guess=nothing)
    
    OmegaH = a ./ (2 .* (GNew .* M) .* (1 .+ sqrt.(1 .- a.^2)))
    alph = mu * GNew * M

    
    if (ergL(n, l, m, mu, M, a) < m .* OmegaH)||(for_s_rates==true)||QNM
        

        if (alph < 0.03)
            alph_ev = 0.03
        else
            alph_ev = alph
        end
        
        if m == 1
            Ntot = 1500
        elseif (m==2)
            Ntot = 4250
        elseif (m==3)
            Ntot = 10000
        else
            Ntot = 4000
        end
        
        if a > 0.95
            Ntot *= 2.5
            Ntot = Int(Ntot)
        end
        
        if Ntot_force > Ntot
            Ntot = Ntot_force
        end

        rescale = (alph ./ alph_ev).^(4 .* l + 5)
        
        b = sqrt.(1 - a.^2)
        
        
        if !QNM
            if (n <= 5)&&(n >= 1)&&(l==0)&&(m==0)&&(alph > 0.03)&&(alph < 2.0)
                erg0_r = open(readdlm, "input_info/ErgEvol_n_$(n)_l_0_m_0_.dat")
                erg0_i = open(readdlm, "input_info/ErgEvol_I_n_$(n)_l_0_m_0_.dat");
                itp = LinearInterpolation(erg0_r[:, 1], erg0_r[:, 2], extrapolation_bc=Line())
                itpI = LinearInterpolation(erg0_i[:, 1], erg0_i[:, 2], extrapolation_bc=Line())
                w0 =  itp(alph) .+ im .* itpI(alph)
            else
                if isnothing(erg_Guess)
                    SR211_g = sr_rates(n, l, m, alph_ev ./ (GNew * M), M, a)
                    # print("test \t ", SR211_g, "\n")
                    if SR211_g == 0
                        SR211_g = -1e-10 .* mu
                    end
                    w0 = (ergL(n, l, m, alph_ev ./ (GNew * M), M, a; full=false) .+ im * SR211_g) .* GNew * M
                else
                    w0 = erg_Guess
                end
            end
        else
            w0 = QNM_E .+ 0.0 * im
        end
        
        
        
        function wrapper!(F, x)
            # wR = 10 .^ x[1]  # real
            wR = x[1]
            wI = x[2]  # imag
            erg = wR + im * wI
    
            if !QNM
                q = - sqrt.(alph_ev.^2 .- erg.^2)
            else
                q = sqrt.(alph_ev.^2 .- erg.^2)
            end
           
            
            gam = im * a * sqrt.(erg.^2 .- alph_ev.^2)
            LLM = l * (l + 1)
            LLM += (-1 + 2 * l * (l + 1) - 2 * m.^2) * gam.^2 ./ (-3 + 4 * l * (l + 1))
            LLM += ((l - m - 1 * (l - m) * (l + m) * (l + m - 1)) ./ ((-3 + 2 * l) * (2 * l - 1).^2) - (l + 1 - m) * (2 * l - m) * (l + m + 1) * (2 + l + m) ./ ((3 + 2 * l).^2 * (5 + 2 * l))) * gam.^4 ./ (2 * (1 + 2 * l))
            LLM += (4 * ((-1 + 4 * m^2) * (l * (1 + l) * (121 + l * (1 + l) * (213 + 8 * l * (1 + l) * (-37 + 10 * l * (1 + l)))) - 2 * l * (1 + l) * (-137 + 56 * l * (1 + l) * (3 + 2 * l * (1 + l))) * m^2 + (705 + 8 * l * (1 + l) * (125 + 18 * l * (1 + l))) * m^4 - 15 * (1 + 46 * m^2))) * gam^6) / ((-5 + 2 * l) * (-3 + 2 * l) * (5 + 2 * l) * (7 + 2 * l) * (-3 + 4 * l * (1 + l))^5)

            
            c0 = 1.0 .- 2.0 * im * erg - 2 * im ./ b .* (erg .- a .* m ./ 2.0)
            c1 = -4.0 .+ 4 * im * (erg - im * q * (1.0 + b)) + 4 * im / b * (erg - a * m / 2.0) .- 2.0 * (erg.^2 .+ q.^2) ./ q
            c2 = 3.0 - 2 * im * erg - 2.0 * (q.^2 - erg.^2) ./ q - 2.0 * im / b * (erg - a * m ./ 2)
            c3 = 2.0 * im * (erg - im * q).^3 ./ q .+ 2 * (erg .- im * q).^2 .* b + q.^2 * a.^2 .+ 2 * im * q * a * m - LLM - 1 - (erg - im * q).^2 ./ q .+ 2 * q * b + 2 * im / b * ( (erg - im * q).^2 ./ q + 1.0) * (erg .- a * m / 2)
            c4 = (erg .- im * q).^4 ./ q.^2 .+ 2 * im * erg * (erg .- im * q).^2 ./ q .- 2 * im ./ b .* (erg .- im * q).^2 ./ q .* (erg .- a * m ./ 2)
           
            function alphaN(nn)
                return nn.^2 .+ (c0 .+ 1) * nn + c0
            end
            function betaN(nn)
                return -2 * nn.^2 .+ (c1 .+ 2) * nn + c3
            end
            function gammaN(nn)
                return nn.^2 .+ (c2 .- 3) * nn + c4
            end
            
            temp = 1
            for i in Ntot:-1:1
                temp = gammaN(i) ./ (betaN(i) .- alphaN(i) .* temp)
            end
            recur = betaN(0) ./ alphaN(0) .- temp
            
            F[1] = real(recur)
            F[2] = imag(recur)
            
            # print(x, "\t", F, "\n")
        end
        
        if real(w0) < 0
            w0 = alph_ev + imag(w0)
        end
        
        # sol = nlsolve(wrapper!, [BigFloat(log10.(real(w0))), BigFloat(imag(w0))], autodiff = :forward, xtol=xtol, ftol=ftol, iterations=iter)
        sol = nlsolve(wrapper!, [BigFloat((real(w0))), BigFloat(imag(w0))], autodiff = :forward, xtol=xtol, ftol=ftol, iterations=iter)
        
        if debug
            print(sol, "\n\n")

          
            Ff = zeros(2)
            wrapper!(Ff, sol.zero)
            print("sanity check \t", Ff, "\n")
        end
        
        if ((!sol.x_converged) && (!sol.f_converged))||(sol.zero[1] .< 0.0)
            if return_both
                return 0.0, 0.0
            else
                return 0.0
            end
        end

        if !return_both
            return sol.zero[2] .* rescale  # im(G * M * omega)
        else
            return sol.zero[1], sol.zero[2]
        end
    else
        if !return_both
            return 0.0
        else
            return 0.0, 0.0
        end
    end
end

function compute_gridded(mu, M, a, n, l, m; Ntot=200, iter=50, xtol=1e-7, npts=30, amin=0.0)
    a_max = a * 1.1 # safety, just in case upward fluctuation
    if a_max > maxSpin
        a_max = maxSpin
    end
    
    output = zeros(npts)
    if amin < a_max
        alist = LinRange(amin, a_max, npts);
    
        alph = GNew * M * mu

        for i in 1:length(alist)
            output[i] = find_im_part(mu, M, alist[i], n, l, m, Ntot_force=Ntot, iter=iter, xtol=xtol) ./ (GNew * M)
        end
        condit = output .<= 0.0
        output[condit] .= 1e-100
    else
        alist = LinRange(a_max, amin, npts);
    end
    return alist, output
end

function generalized_laguerre(n, α, x)
    return (-1).^ n ./ factorial(big(n)) .* HypergeometricFunctions.U.(-n, α+1 ,x)
    # return binomial(n + α, n) .* HypergeometricFunctions.M.(-n, α+1 , x)
end

function test_projection_scatter(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=20000, Npts_Bnd=4000, rmaxT=100, debug=false, Ntot_safe=5000, sve_for_test=false)
    
    rp = 1 + sqrt.(1 - a.^2)
    alph = mu * GNew * M
    
    # approx energy estimates
    # erg_1G = ergL(n1, l1, m1, mu, M, a)
    # erg_2G = ergL(n2, l2, m2, mu, M, a)
    # erg_3G = ergL(n3, l3, m3, mu, M, a)
    # erg_pxy = (erg_1G + erg_2G - erg_3G) * GNew * M
    
    wR1, wI1 = find_im_part(mu, M, a, n1, l1, m1; return_both=true, for_s_rates=true)
    wR2, wI2 = find_im_part(mu, M, a, n2, l2, m2; return_both=true, for_s_rates=true)
    wR3, wI3 = find_im_part(mu, M, a, n3, l3, m3; return_both=true, for_s_rates=true)
    erg_pxy = (wR1 + wR2 - wR3)
    
    qq = abs.(real(- sqrt.(alph.^2 .- erg_pxy.^2)))
    rmax = rmaxT ./ qq
    
    # rlist = 10 .^(range(log10.(rp), log10.(rmax), rpts))
    rlist = range(rp, rmax, rpts)
    
    theta_rand = 2 * acos.(sqrt.(rand(rpts)))
    # r_rand = rand(rpts).^(1/3) * (rmax^3 ./ 3 - 1.0 / 3)^(1/3)
    sampls = rand(rpts)
    r_rand = rp .+ sampls .* (rmax .- rp)
    exF = (rmax .- rp)
    
    alph = mu .* GNew .* M
    
    
    
    erg_out = ergL(1, 0, 0, mu, M, a)
    erg_diff = abs.(-(erg_pxy.^2 .- (erg_out .* GNew .* M).^2))
    
    idx_run = 1
    for i in 2:30
        erg_out = ergL(i, 0, 0, mu, M, a)
        erg_diff_hold = abs.(-(erg_pxy.^2 .- (erg_out .* GNew .* M).^2))
        if erg_diff_hold < erg_diff
            erg_diff = erg_diff_hold
            idx_run = i
        end
    end
    
    
    for i in idx_run:idx_run
        rl, r1, erg = solve_radial(mu, M, a, i, 0, 0; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe)
        
        itp_r = LinearInterpolation(log10.(Float64.(rl)), real(r1), extrapolation_bc=Line())
        itp_im = LinearInterpolation(log10.(Float64.(rl)), imag(r1), extrapolation_bc=Line())
        function rad_comp(r)
            return itp_r(log10.(r)) .+ im * itp_im(log10.(r))
        end
        

        #### TESTING
        trapz(y,x) = @views sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
        # erg = ergL(i, 0, 0, mu, M, a) .* GNew .* M
        erg = erg_pxy
        
        r = LinRange(rp .* 1.01, rmax, rpts)
        h = r[2] .- r[1]
        
        alph = mu .* GNew .* M
        k = sqrt.(Complex.(erg.^2 .- alph.^2))
                
        delt = r.^2 .- 2 .* r .+ a.^2
        dr_delt = 2 .* r .- 2
        m = 0
        # rhs = real((erg.^2 .* (r.^2 .+ a.^2).^2 .- 4 .* a .* m .* erg .* r .+ m.^2 .* a.^2) ./ delt .- (erg.^2 .* a.^2 .+ alph.^2 .* r.^2))
        rhs = ((erg.^2 .* (r.^2 .+ a.^2).^2 .- 4 .* a .* m .* erg .* r .+ m.^2 .* a.^2) ./ delt .- (erg.^2 .* a.^2 .+ alph.^2 .* r.^2))
        
        # Create the second derivative matrix (Laplacian matrix), and first derivative
        D2 = zeros(rpts, rpts)

        # Fill the diagonal
        for i in 2:rpts-1
            D2[i, i-1] += 1.0 .* delt[i] ./ h^2
            D2[i, i] += -2.0 .* delt[i] ./ h^2
            D2[i, i+1] += 1.0 .* delt[i] ./ h^2
            
            D2[i, i+1] += 1.0 .* dr_delt[i] ./ (2 * h)
            D2[i, i-1] += -1.0 .* dr_delt[i] ./ (2 * h)
        end

        # BNDRY TERMS
        D2[1, 1] += -1.0 .* dr_delt[1] ./ (2 * h)
        D2[1, 2] += 1.0 .* dr_delt[1] ./ (2 * h)
        # D2[1, 3] += -1.0 .* dr_delt[1] ./ (2 * h)
        D2[rpts, rpts] += 1.0 .* dr_delt[rpts] ./ (2 * h)
        D2[rpts, rpts - 1] += -1.0 .* dr_delt[rpts] ./ (2 * h)
        # D2[rpts, rpts - 2] += 1.0 .* dr_delt[rpts] ./ (2 * h)
        
        D2[1, 1] += 1.0 .* delt[1] ./ h.^2
        D2[1, 2] += -2.0 .* delt[1] ./ h.^2
        D2[1, 3] += 1.0 .* delt[1] ./ h.^2
        # D2[1, 4] += -1.0 .* delt[1] ./ h.^2
        D2[rpts, rpts] += 1.0 .* delt[rpts] ./ h.^2
        D2[rpts, rpts - 1] += -2.0 .* delt[rpts] ./ h.^2
        D2[rpts, rpts - 2] += 1.0 .* delt[rpts] ./ h.^2
        # D2[rpts, rpts - 3] += -1.0 .* delt[rpts] ./ h.^2
        
        r_TEST = rad_comp(r)
        temp = (D2 * r_TEST) .+ rhs .* r_TEST
      
        out = Float64(real(-trapz(temp .* conj(r_TEST) , r)))
 
        ################
        
        
        # writedlm("test_store/check_real_1.dat", hcat(real(rlist), float(real(itp_dr(log10.(rlist))))))
        # writedlm("test_store/check_im_1.dat", hcat(real(rlist), float(imag(itp_dr(log10.(rlist))))))
        # writedlm("test_store/check_real_2.dat", hcat(real(rlist), float(real(itp_dr2(log10.(rlist))))))
        # writedlm("test_store/check_im_2.dat", hcat(real(rlist), float(imag(itp_dr2(log10.(rlist))))))
        

        # print(integrd, "\n")
        if sve_for_test
            writedlm("test_store/check_real_1.dat", hcat(real(r_rand), float(real(d_t2phi .+ dr_term .+ theta_term .+ mu_term .* Z_outC .* conj(rad_comp(r_rand))))))
        end

        
        erg_out = ergL(i, 0, 0, mu, M, a)
        erg_diff_pxy = -(erg_pxy.^2 .- (erg_out .* GNew .* M).^2)
        
        
        
#        if debug
#            if i == 1
#                print(i, "\t", -11 * alph^4 ./ 18, "\t", erg_diff_pxy, "\t", out, "\n")
#            elseif i == 2
#                print(i, "\t", 5 * alph^4 ./ 36, "\t", erg_diff_pxy, "\t", out, "\n")
#            end
#        end


        return erg_diff_pxy, out
        
    end
end

function direct_solve_radialL1(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=1000, Npts_Bnd=1000, rmaxT=100, debug=false, Ntot_safe=5000, sve_for_test=false,  iter=10, xtol=1e-10, ftol=1e-10, tag="_", Nang=100000, eps_r=1e-10)
    
    rp = 1 + sqrt.(1 - a.^2)
    alph = mu * GNew * M
    
    erg_1G = ergL(n1, l1, m1, mu, M, a)
    erg_2G = ergL(n2, l2, m2, mu, M, a)
    erg_3G = ergL(n3, l3, m3, mu, M, a)
    erg_pxy = erg_1G + erg_2G - erg_3G
   
    q = abs.(real(- sqrt.(abs.(alph.^2 .- (erg_pxy .* GNew .* M).^2))))
    rmax = rmaxT ./ q
    rlist = range(rp .* (1.0 .+ eps_r), rmax, rpts)
    h = rlist[2] - rlist[1]
    
    rl, r1, erg_1 = solve_radial(mu, M, a, n1, l1, m1; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true)
    itp = LinearInterpolation(log10.(rl), r1, extrapolation_bc=Line())
    rf_1 = itp(log10.(rlist))
    
    
    rl, r2, erg_2 = solve_radial(mu, M, a, n1, l1, m1; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true)
    itp = LinearInterpolation(log10.(rl), r2, extrapolation_bc=Line())
    rf_2 = itp(log10.(rlist))
    
    rl, r3, erg_3 = solve_radial(mu, M, a, n1, l1, m1; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true)
    itp = LinearInterpolation(log10.(rl), r3, extrapolation_bc=Line())
    rf_3 = itp(log10.(rlist))
    
    erg = (erg_1 + erg_2 - erg_3)
    
    
    Z1 = spheroidals(l1, m1, a, erg_1)
    Z2 = spheroidals(l2, m2, a, erg_2)
    Z3 = spheroidals(l3, m3, a, erg_3)
    Z4 = spheroidals(0, 0, a, erg)
   
    thetaV = acos.(1.0 .- 2.0 .* rand(Nang))
    phiV = rand(Nang) .* 2*pi

    function func_ang(x)
        return real(Z1.(x[1], x[2]) .* Z2.(x[1], x[2]) .* conj(Z3.(x[1], x[2])) .* conj(Z4.(x[1], x[2]))) #
    end
    CG = 0.0
    for i in 1:Nang
        CG += func_ang([thetaV[i], phiV[i]])
    end
    CG *= 4*pi / Nang
    
    
    delt = rlist.^2 .- 2 .* rlist .+ a.^2
    dr_delt = 2 .* rlist .- 2
    m = 0
    # assumes l = m = 0
    rhs = ((erg.^2 .* (rlist.^2 .+ a.^2).^2 .- 4 .* a .* m .* erg .* rlist .+ m.^2 .* a.^2) ./ delt .- (erg.^2 .* a.^2 .+ alph.^2 .* rlist.^2))
    
    # Create the second derivative matrix (Laplacian matrix), and first derivative
    D2 = zeros(rpts, rpts)

    # Fill the diagonal
    for i in 2:rpts-1
        D2[i, i-1] += 1.0 .* delt[i] ./ h^2
        D2[i, i] += -2.0 .* delt[i] ./ h^2
        D2[i, i+1] += 1.0 .* delt[i] ./ h^2
        
        D2[i, i+1] += 1.0 .* dr_delt[i] ./ (2 * h)
        D2[i, i-1] += -1.0 .* dr_delt[i] ./ (2 * h)
    end

    # BNDRY TERMS
    D2[1, 1] += -1.0 .* dr_delt[1] ./ (2 * h)
    D2[1, 2] += 1.0 .* dr_delt[1] ./ (2 * h)
    # D2[1, 3] += -1.0 .* dr_delt[1] ./ (2 * h)
    D2[rpts, rpts] += 1.0 .* dr_delt[rpts] ./ (2 * h)
    D2[rpts, rpts - 1] += -1.0 .* dr_delt[rpts] ./ (2 * h)
    # D2[rpts, rpts - 2] += 1.0 .* dr_delt[rpts] ./ (2 * h)
    
    D2[1, 1] += 1.0 .* delt[1] ./ h.^2
    D2[1, 2] += -2.0 .* delt[1] ./ h.^2
    D2[1, 3] += 1.0 .* delt[1] ./ h.^2
    # D2[1, 4] += -1.0 .* delt[1] ./ h.^2
    D2[rpts, rpts] += 1.0 .* delt[rpts] ./ h.^2
    D2[rpts, rpts - 1] += -2.0 .* delt[rpts] ./ h.^2
    D2[rpts, rpts - 2] += 1.0 .* delt[rpts] ./ h.^2
    # D2[rpts, rpts - 3] += -1.0 .* delt[rpts] ./ h.^2
    
    
    max_states = 2
    hold_r_list = zeros(Complex, rpts, max_states)
    for i in 1:max_states
        rl, rout, erg_out = solve_radial(mu, M, a, i, 0, 0; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true)
        itp = LinearInterpolation(log10.(rl), rout, extrapolation_bc=Line())
        r_on_targ = itp(log10.(rlist))
        hold_r_list[:, i] .= r_on_targ
    end
    
    unitMatch = 1.0 ./ (GNew .* M .* alph.^(3/2))
    
  
    function wrapper(x)
        WF  = zeros(Complex, rpts)
        for i in 1:max_states
            WF .+= hold_r_list[:, i] .* x[i]
        end
        temp = ((D2 * WF) .+ rhs .* WF) ./ rlist.^2 .+ CG .* (rf_1 .* rf_2 .* rf_3) .* unitMatch ./ 6.0
        # F .= abs.(temp)
        return sum(abs.(temp))
    end
    
    # rGuess = rf_1
    # r_in = real(rGuess)
    c_in = ones(max_states) .* 1e-5
    
    # test
#    c_in = ones(max_states)
#    c_in[1] = -0.0002 .* alph.^(9/2)
#    c_in[2] = -0.00004 .* alph.^(9/2)
#    c_in[3] = -0.00001 .* alph.^(9/2)
#    c_in[4] = -0.000008 .* alph.^(9/2)
#    c_in[5] = -0.000005 .* alph.^(9/2)
#    c_in[6] = -0.000004 .* alph.^(9/2)
#    c_in[7] = -0.000003 .* alph.^(9/2)
#    c_in[8] = -0.000002 .* alph.^(9/2)
#    c_in[9] = -0.000002 .* alph.^(9/2)
#    c_in[10] = -0.000001 .* alph.^(9/2)
    
    print("initial \t", c_in, "\n")
    # sol = nlsolve(wrapper!, c_in, show_trace=debug, autodiff = :forward, xtol=xtol, ftol=ftol, iterations=iter)
    # full_out = sol.zero
    sol = optimize(wrapper, c_in, LBFGS(), Optim.Options(g_tol = 1e-11, iterations=iter, show_trace=true, show_every=1))
    full_out = Optim.minimizer(sol)
    print(full_out, "\n")
    
    Psi0 = 0.0
    for i in 1:max_states
        Psi0 += hold_r_list[1, i] .* full_out[i]
    end
    return Psi0 .* conj(Psi0)
    
#    trapz(y,x) = @views sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
#    # nm = trapz(full_out .* conj(full_out) .* rlist.^2, rlist)
#    # print(trapz(full_out .* conj(full_out) .* rlist.^2, rlist), "\n")
#    nm = 1.0
    
#    zero_pt = (full_out ./ sqrt.(nm))
    
    writedlm("test_store/test_radial_induced"*tag*".dat", hcat(real(rlist), float((full_out))))
    
    
#    itp = LinearInterpolation(log10.(rlist[2:end]), log10.(abs.(zero_pt[2:end])), extrapolation_bc=Line())
#    psi_1 = 10 .^itp(log10.(rp))
    
#    return psi_1.^2
    
end

function integrate_radialEq(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=1000, Npts_Bnd=1000, rmaxT=100, debug=false, Ntot_safe=5000, sve_for_test=false,  iter=10, xtol=1e-10, ftol=1e-10, tag="_", Nang=100000, n_times=300000, add_source=true, eps_fac = 1e-10, dtmax_slv, m=0, l=0, reltol=1e-8, NON_REL=false)
    
    rp = BigFloat(1.0 .+ sqrt.(1.0 .- a.^2))
    rmm = BigFloat(1.0 .- sqrt.(1.0 .- a.^2))
    alph = mu * GNew * M
    

    
    erg_1G = ergL(n1, l1, m1, mu, M, a; full=false)
    erg_2G = ergL(n2, l2, m2, mu, M, a; full=false)
    erg_3G = ergL(n3, l3, m3, mu, M, a; full=false)
    erg_pxy = erg_1G + erg_2G - erg_3G
    
    maxN = maximum([n1 n2 n3])
    
    rmax = Float64.(50 ./ alph.^2)
    
    rlist = range(rp .* (1.0 .+ eps_fac), rmax, rpts)
    
    # rl, r1, erg_1 = solve_radial(mu, M, a, n1, l1, m1; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe)
    # itp = LinearInterpolation(log10.(rl .* (GNew .* M)), r1 ./ (GNew .* M).^(3/2), extrapolation_bc=Line())
    # rf_1 = itp(log10.(rlist))
    rf_1 = radial_bound_NR(n1, l1, m1, mu, M, rlist; physU=false)
    val1 = sr_rates(n1, l1, m1,  mu, M, a, impose_low_cut=1e-10) * (GNew .* M)
    erg_1 = erg_1G * GNew * M .+ val1 .* im
    
    
    # rl, r2, erg_2 = solve_radial(mu, M, a, n2, l2, m2; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe)
    # itp = LinearInterpolation(log10.(rl .* (GNew .* M)), r2 ./ (GNew .* M).^(3/2), extrapolation_bc=Line())
    # rf_2 = itp(log10.(rlist))
    rf_2 = radial_bound_NR(n2, l2, m2, mu, M, rlist; physU=false)
    val2 = sr_rates(n2, l2, m2,  mu, M, a, impose_low_cut=1e-10) * (GNew .* M)
    erg_2 = erg_2G * GNew * M .+ val2 .* im
    
    
    # rl, r3, erg_3 = solve_radial(mu, M, a, n3, l2, m3; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe)
    # itp = LinearInterpolation(log10.(rl .* (GNew .* M)), r3 ./ (GNew .* M).^(3/2), extrapolation_bc=Line())
    # rf_3 = itp(log10.(rlist))
    rf_3 = radial_bound_NR(n3, l3, m3, mu, M, rlist; physU=false)
    val3 = sr_rates(n3, l3, m3,  mu, M, a, impose_low_cut=1e-10) * (GNew .* M)
    erg_3 = erg_3G * GNew * M .+ val3 .* im


    # wR1, wI1 = find_im_part(mu, M, a, n1, l1, m1; debug=false, return_both=true, for_s_rates=true, Ntot_force=10000)
    # wR2, wI2 = find_im_part(mu, M, a, n2, l2, m2; debug=false, return_both=true, for_s_rates=true, Ntot_force=10000)
    # wR3, wI3 = find_im_part(mu, M, a, n3, l3, m3; debug=false, return_both=true, for_s_rates=true, Ntot_force=10000)
    
    # erg = ((wR1 .+ wR2 .- wR3) .+ im * (wI1 .+ wI2 .- wI3))   #
    erg = (erg_1 + erg_2 - conj(erg_3))  #
    print("ERG \t", Float64.(real(erg)), "\t [Im]", Float64.(imag(erg)), "\n")

    Z1 = spheroidals(l1, m1, a, erg_1)
    Z2 = spheroidals(l2, m2, a, erg_2)
    Z3 = spheroidals(l3, m3, a, erg_3)
    Z4 = spheroidals(l, m, a, erg)
   
    thetaV = acos.(1.0 .- 2.0 .* rand(Nang))
    phiV = rand(Nang) .* 2*pi

    function func_ang(x)
        return real(Z1.(x[1], x[2]) .* Z2.(x[1], x[2]) .* conj(Z3.(x[1], x[2])) .* conj(Z4.(x[1], x[2])))
    end
    function func_ang_2(x)
        return real(Z1.(x[1], x[2]) .* Z2.(x[1], x[2]) .* conj(Z3.(x[1], x[2])) .* conj(Z4.(x[1], x[2]))) .* a.^2 .* cos.(x[1]).^2
    end
    CG = 0.0
    CG_2 = 0.0
    for i in 1:Nang
        CG += func_ang([thetaV[i], phiV[i]])
        CG_2 += func_ang_2([thetaV[i], phiV[i]])
    end
    CG *= 4*pi / Nang
    CG_2 *= 4*pi / Nang
    
    print("CG / CG2 \t", CG, "\t", CG_2, "\n")
    print("rp  ", Float64.(rp), "\t alpha \t", alph, "\n")
    
    lam_eff = 1.0
    if (n1==n2)&&(l1==l2)&&(m1==m2)
        preFac = 3
    else
        preFac = 6
    end
    unitMatch = -1.0 ./ (2 .* alph).^(3/2) .* lam_eff
    gammaT = (preFac .* (rf_1 .* rf_2 .* conj(rf_3)) ./ 6.0 .* unitMatch) ./ (GNew .* M)
    print("Max test \t", Float64.(maximum(real(abs.(gammaT)))), "\t", Float64.(rlist[argmax(real(abs.(gammaT)))]), "\n")
    print("r1/r2/r3 \t", Float64(maximum(real(abs.(rf_1)))), "\t", Float64(maximum(real(abs.(rf_2)))), "\t", Float64(maximum(real(abs.(rf_3)))), "\n")
    
    itpG = LinearInterpolation(log10.(rlist), Float64.(real.(gammaT)), extrapolation_bc=Line())
    itpGI = LinearInterpolation(log10.(rlist), Float64.(imag.(gammaT)), extrapolation_bc=Line())
    # writedlm("test_store/source_test.dat", hcat(rlist, Float64.(real.(gammaT))))
    # writedlm("test_store/source_test_i.dat", hcat(rlist, Float64.(imag.(gammaT))))

    Mvars = [1.0]
    
    if sve_for_test
        nn = l + 1
        wR, wI = find_im_part(mu, M, a, nn, l, m; debug=false, iter=500, xtol=1e-20, ftol=1e-90, return_both=true, for_s_rates=true)
        erg = wR .+ im .* wI
    end
    
    
    gam = im * a * sqrt.(erg.^2 .- alph.^2)
    LLM = l * (l + 1)
    LLM += (-1 + 2 * l * (l + 1) - 2 * m.^2) * gam.^2 ./ (-3 + 4 * l * (l + 1))
    LLM += ((l - m - 1 * (l - m) * (l + m) * (l + m - 1)) ./ ((-3 + 2 * l) * (2 * l - 1).^2) - (l + 1 - m) * (2 * l - m) * (l + m + 1) * (2 + l + m) ./ ((3 + 2 * l).^2 * (5 + 2 * l))) * gam.^4 ./ (2 * (1 + 2 * l))
    LLM += (4 * ((-1 + 4 * m^2) * (l * (1 + l) * (121 + l * (1 + l) * (213 + 8 * l * (1 + l) * (-37 + 10 * l * (1 + l)))) - 2 * l * (1 + l) * (-137 + 56 * l * (1 + l) * (3 + 2 * l * (1 + l))) * m^2 + (705 + 8 * l * (1 + l) * (125 + 18 * l * (1 + l))) * m^4 - 15 * (1 + 46 * m^2))) * gam^6) / ((-5 + 2 * l) * (-3 + 2 * l) * (5 + 2 * l) * (7 + 2 * l) * (-3 + 4 * l * (1 + l))^5)
    
    
    higherO = false
    
    function radial_diff!(du, u, Mvars, rstar)
        Rr = u[1]
        dRr = u[2]
        # dRr = 0.0
        
        if NON_REL
            
            # NR - limit, rstar = r
            r = rstar
            SS = (itpG(log10.(r)) + im * itpGI(log10.(r)))
            
            delt = (r.^2 .- 2 .* r .+ a.^2)
        
            # secDer = 2 .* alph.^2 .* Rr ./ r .- 4 .* erg.^2 .* Rr ./ r .+ (alph.^2 .- erg.^2) .* Rr  .- 2 .* dRr ./ r .+ LLM ./ (2 .* r) .* Rr
            secDer = -2 .* alph.^2 .* Rr ./ r .+ (alph.^2 .- erg.^2) .* Rr  .- 2 .* dRr ./ r .+ LLM ./ (2 .* r) .* Rr
            
            if add_source
                secDer += SS .* CG .* (1 .+ 2.0 ./ r)
                # secDer += SS .* CG
            end
            
            
            
            if higherO
                if add_source
                    secDer += SS .* ((4 .- a.^2) .* CG  .+ a.^2 .* CG_2) ./ r.^2
                end
                secDer += (4 .* alph.^2 .- a.^2 .* alph.^2 .- 12 .* erg.^2 .+ a.^2 .* erg.^2) .* Rr ./ r.^2
                secDer += - 2 * dRr ./ r.^2
                
                # 3rd ord
                if add_source
                    secDer += SS .* ((8 .- 8 .* a.^2) .* CG  .+ 2 .* a.^2 .* CG_2) ./ r.^3
                end
                secDer += (8 .* alph.^2 .- 4 .* a.^2 .* alph.^2 .- 32 .* erg.^2 .+ 6 .* a.^2 .* erg.^2) .* Rr ./ r.^3
                secDer += (-4 .+ 2 .* a.^2) .* dRr ./ r.^3
                
                # 4th ord
                if add_source
                    secDer += SS  .* ((16 .- 12 .* a.^2 .+ a.^4) .* CG  .+ (4 .* a.^2 .- a.^4) .* CG_2) ./ r.^4
                end
                secDer += (16 .* alph.^2 .- 12 .* a.^2 .* alph.^2 .+ a.^4 .* alph.^2 .- 80 .* erg.^2 .+ 28 .* a.^2 .* erg.^2 .- a.^4 .* erg.^2) .* Rr ./ r.^4
                secDer += (-8 .+ 6 .* a.^2) .* dRr ./ r.^4
                
                # 5th ord
                if add_source
                    secDer += SS  .* ((32 .- 32 .* a.^2 .+ 6 .* a.^4) .* CG  .+ (8 .* a.^2 .- 4 .* a.^4) .* CG_2) ./ r.^5
                end
                secDer += (32 .* alph.^2 .- 32 .* a.^2 .* alph.^2 .+ 6 .* a.^4 .* alph.^2 .- 192 .* erg.^2 .+ 104 .* a.^2 .* erg.^2 .- 8 .* a.^4 .* erg.^2) .* Rr ./ r.^5
                secDer += (-16 .+ 16 .* a.^2 .- 2 .* a.^4) .* dRr ./ r.^5
                
                # 6th ord
                if add_source
                    secDer += SS .* ((64 .- 80 .* a.^2 .+ 24 .* a.^4 .- a.^6) .* CG  .+ (16 .* a.^2 .- 12 .* a.^4 .+ a.^6) .* CG_2) ./ r.^6
                end
                secDer += (64 .* alph.^2 .- 80 .* a.^2 .* alph.^2 .+ 24 .* a.^4 .* alph.^2 .- a.^6 .* alph.^2 .- 448 .* erg.^2 .+ 336 .* a.^2 .* erg.^2 .- 48 .* a.^4 .* erg.^2 .+ a.^6 .* erg.^2) .* Rr ./ r.^6
                secDer += (-32 .+ 40 .* a.^2 .- 10 .* a.^4) .* dRr ./ r.^6
                
            end
            
        else
            
            # Full R - limit
            r = itp_rrstar(rstar)
            delt = (r.^2 .- 2 .* r .+ a.^2)
            SS = (itpG(log10.(r)) + im * itpGI(log10.(r)))
            ff = (r.^2 .+ a.^2)
            secDer = (delt .* (LLM .+ r.^2 .* alph.^2) .- r .* (r.^3 .+ a.^2 .* (r .+ 2)) .* erg.^2) ./ delt.^2
            secDer += -2 .* (r .- 1) ./ delt .* dRr

            if add_source
                secDer += SS .* (CG .* r.^2 .+ a.^2 .* CG_2) .* delt ./ ff.^(3/2)
            end
            
        end
        du[1] = u[2]
        du[2] = secDer
    end
  
    # rstart = (itpG(log10.(rmax)) + im * itpGI(log10.(rmax))) .* CG .* (1 .+ 2.0 ./ rmax) ./ rmax.^2
    rstart = 1e-50 .* (1 + im)
    y0 = [-rstart rstart ./ 1e3]
    print("Init conds \t", y0, "\n")
    
    tspan_r = (Float64.(rmax), Float64.(rp * (1.0 .+ eps_fac)))
    
    if !NON_REL
        tspan_r = (BigFloat(rmax), BigFloat(rp * (1.0 .+ eps_fac)))
        r_list_map = 10 .^LinRange(log10.(tspan_r[2]), log10.(tspan_r[1]), 100000)

        rout_star = r_list_map .+ 2.0 .* rp ./ (rp .- rmm) .* log.(r_list_map ./ rp .- 1.0) .- 2.0 .* rmm .* log.(r_list_map ./ rmm .- 1.0) ./ (rp .- rmm)
        itp_rrstar = LinearInterpolation(rout_star, r_list_map, extrapolation_bc=Line())
        tspan_r = (rout_star[end], rout_star[1])
        saveat = (tspan_r[2] .- tspan_r[1]) ./ n_times
    end
    print("Range \t ", tspan_r, "\n")
    saveat = (tspan_r[2] .- tspan_r[1]) ./ n_times
    
    
    function check_dt(u, t, integrator)
        if t < 3.0 * rp
            return true
        else
            return false
        end
    end
    function affect_dt!(integrator)
        if integrator.opts.dtmax > (1e-2 * rp)
           integrator.opts.dtmax *= 0.9
        end
    end
    cbackdt = DiscreteCallback(check_dt, affect_dt!, save_positions=(false, true))
    cbset = CallbackSet(cbackdt)
    
    prob = ODEProblem(radial_diff!, y0, tspan_r, Mvars, reltol=reltol, abstol=1e-80)
    # sol = solve(prob, Rodas4P(autodiff=false), saveat=saveat, maxiters=5e6, dtmax=dtmax_slv, callback=cbset)
    # sol = solve(prob, Rodas4P(autodiff=false), saveat=saveat, maxiters=5e6, dtmax=rp/10.0)
    # prob = ODEProblem(radial_diff!, y0, tspan_r, Mvars, reltol=reltol, abstol=1e-80)
    # sol = solve(prob, Vern6(), saveat=saveat, maxiters=5e6, dtmax=rp/1000.0)
    sol = solve(prob, KenCarp4(autodiff=false), saveat=saveat, maxiters=5e6, dtmax=rp/100.0)
    
    
    
    
  
    radVals = reverse(sol.t)
    if !NON_REL
        # WF = reverse([sol.u[i][1] ./ sqrt.(itp_rrstar(sol.t[i]).^2 .+ a.^2) ./ (GNew .* M) for i in 1:length(sol.u)])
        WF = reverse([sol.u[i][1]   for i in 1:length(sol.u)])
        
        
        writedlm("test_store/test_RAD.dat", hcat(itp_rrstar.(radVals), float(real(WF .* conj.(WF)))))
        # writedlm("test_store/test_RAD_rel_2.dat", hcat(radVals, float(real(WF .* conj.(WF)))))
        
        # writedlm("test_store/test_Rreal.dat", hcat(radVals, float(real(WF))))
        # dWF = reverse([sol.u[i][2]  for i in 1:length(sol.u)])
        # writedlm("test_store/test_dRreal.dat", hcat(radVals, float(real(dWF))))
        
        maxV = maximum(real(WF .* conj.(WF)))
        
        out = Float64.(real(8 * pi .* maxV .* erg.^2 ./ mu.^2 .* (1 .+ sqrt.(1 .- a.^2)) ./ lam_eff.^2))
        
    else
        ## TESTING
        WF = reverse([sol.u[i][1] for i in 1:length(sol.u)])
        # WF = reverse([sol.u[i][1] ./ sol.t[i]  for i in 1:length(sol.u)])
        # WF_new = zeros(Complex, size(WF))
        # WF_new[end] = WF_p[end]
        # WF_new[end-1] = WF_p[end-1]
        # for i in 1:(length(radVals) - 2)
        #     h1 = radVals[end - i - 1] - radVals[end - i]
        #     h2 = radVals[end - i] - radVals[end - i + 1]
        #     SS = (itpG(log10.(radVals[end - i])) + im * itpGI(log10.(radVals[end - i])))
        #
        #     WF_new[end - 1 - i] = WF_new[end - i] .* (1 .+ h1 ./ h2) .- WF_new[end + 1 - i] .* (h1 ./ h2) .+ (h1.^2 .+ h1 .* h2) .* ((mu.^2 .- (erg ./ (GNew .* M)).^2) .* WF_new[end - i] .- 2 .* alph .* mu ./ radVals[end - i] .* WF_new[end - i] .- SS .* CG .* radVals[end - i]) ./ 2
            
        #     test = -WF_p[end - 1 - i] .+  WF_p[end - i] .* (1 .+ h1 ./ h2)  .- WF_p[end + 1 - i] .* (h1 ./ h2) .+ (h1.^2 .+ h1 .* h2) .* ( (mu.^2 .- (erg ./ (GNew .* M)).^2) .* WF_p[end - i] .- 2 .* alph .* mu ./ radVals[end - i] .* WF_p[end - i] .- SS .* CG .* radVals[end - i] ) ./ 2
            
        #     test2 = -WF_new[end - 1 - i] .+  WF_new[end - i] .* (1 .+ h1 ./ h2)  .- WF_new[end + 1 - i] .* (h1 ./ h2) .+ (h1.^2 .+ h1 .* h2) .* ( (mu.^2 .- (erg ./ (GNew .* M)).^2) .* WF_new[end - i] .- 2 .* alph .* mu ./ radVals[end - i] .* WF_new[end - i] .- SS .* CG .* radVals[end - i] ) ./ 2
            
        #     print(i, "\t", radVals[end-i], "\t", Float64(real(WF_p[end - i])), "\t", Float64(real(WF_new[end - 1 - i])), "\t", Float64(real(test)), "\t", Float64(real(test2)), "\n")
        # end
        # WF_new ./= radVals
        
    
        writedlm("test_store/test_RAD_NR.dat", hcat(radVals, float(real(WF .* conj.(WF)))))
        
        # print("Deets \t", mu.^2 .- (erg ./ GNew .* M).^2, "\n")
        writedlm("test_store/test_RAD.dat", hcat(radVals, float(real(WF))))
        # writedlm("test_store/test_RAD_i.dat", hcat(radVals, float(real(WF_new))))
        
        
        maxV = maximum(real(WF .* conj.(WF)))
        # maxV2 = maximum(real(WF_new .* conj.(WF_new)))
        # print("Max ratio \t ", maxV2 ./ maxV, "\n")
        
        out = Float64.(real(8 * pi .* maxV .* erg.^2 ./ mu.^2 .* (1 .+ sqrt.(1 .- a.^2)) ./ lam_eff.^2))
        
    end
    
    if sve_for_test
        nn = l + 1
        rr_test, outtest, erg_3 = solve_radial(mu, M, a, nn, l, m; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe, eps_r=eps_fac)
        rNR = radial_bound_NR(nn, l, m, mu, M, rr_test)
        writedlm("test_store/test_NR.dat", hcat(rr_test, rNR.^2))
        writedlm("test_store/test_R.dat", hcat(rr_test, real(outtest .* conj.(outtest))))
        trapz(y,x) = @views sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
        if !NON_REL
            # nm = trapz(real(WF .* conj.(WF)) .* itp_rrstar.(radVals).^2, itp_rrstar.(radVals))
            # writedlm("test_store/test_Out.dat", hcat(itp_rrstar.(radVals), real(WF .* conj.(WF)) ./ nm))
            nm = trapz(real(WF .* conj.(WF)) .* radVals.^2, radVals)
            writedlm("test_store/test_Out.dat", hcat(radVals , real(WF .* conj.(WF)) ./ nm))
        else
            nm = trapz(real(WF .* conj.(WF)) .* radVals.^2, radVals)
            writedlm("test_store/test_Out_NR.dat", hcat(radVals, real(WF .* conj.(WF)) ./ nm))
        end
    end
    
    
    print("Output: ", out, "\n")
    
    # theory-NR
    test = 4.3e-7 * alph.^(11 - 4) .* (1 .+ sqrt.(1 .- a.^2))
    # test = 1.1e-9 * alph.^(7 - 4) .* (1 .+ sqrt.(1 .- a.^2))
    # test2 = 7.8e-11 * alph.^(7 - 4) .* (1 .+ sqrt.(1 .- a.^2))
    print("NR Rate: ", test, "\n")
    #print("NR Rate 2: ", test2, "\n")
    print("Ratio : ", out ./ test, "\n")
    #print("Ratio : ", out ./ test2, "\n")
    
    return out
    
end


function integrate_radialEq_2(mu, M, a, n1, l1, m1, n2, l2, m2, n3, l3, m3; rpts=1000, Npts_Bnd=1000, rmaxT=100, debug=false, Ntot_safe=5000, sve_for_test=false,  iter=10, xtol=1e-10, ftol=1e-10, tag="_", Nang=100000, n_times=300000, add_source=true, eps_fac = 1e-10, dtmax_slv, m=0, l=0, reltol=1e-8, NON_REL=false, h_mve=100)
    
    rp = BigFloat(1.0 .+ sqrt.(1.0 .- a.^2))
    rmm = BigFloat(1.0 .- sqrt.(1.0 .- a.^2))
    alph = mu * GNew * M
    

    
    erg_1G = ergL(n1, l1, m1, mu, M, a; full=false)
    erg_2G = ergL(n2, l2, m2, mu, M, a; full=false)
    erg_3G = ergL(n3, l3, m3, mu, M, a; full=false)
    erg_pxy = erg_1G + erg_2G - erg_3G
    
    maxN = maximum([n1 n2 n3])
    
    
    rmax = Float64.(70 ./ alph.^2)
    
    rlist = range(rp .* (1.0 .+ eps_fac), rmax, rpts)
    
    # rl, r1, erg_1 = solve_radial(mu, M, a, n1, l1, m1; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe)
    # itp = LinearInterpolation(log10.(rl .* (GNew .* M)), r1 ./ (GNew .* M).^(3/2), extrapolation_bc=Line())
    # rf_1 = itp(log10.(rlist))
    rf_1 = radial_bound_NR(n1, l1, m1, mu, M, rlist; physU=false)
    val1 = sr_rates(n1, l1, m1,  mu, M, a, impose_low_cut=1e-10) * (GNew .* M)
    erg_1 = erg_1G * GNew * M .+ val1 .* im
    
    
    # rl, r2, erg_2 = solve_radial(mu, M, a, n2, l2, m2; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe)
    # itp = LinearInterpolation(log10.(rl .* (GNew .* M)), r2 ./ (GNew .* M).^(3/2), extrapolation_bc=Line())
    # rf_2 = itp(log10.(rlist))
    rf_2 = radial_bound_NR(n2, l2, m2, mu, M, rlist; physU=false)
    val2 = sr_rates(n2, l2, m2,  mu, M, a, impose_low_cut=1e-10) * (GNew .* M)
    erg_2 = erg_2G * GNew * M .+ val2 .* im
    
    
    # rl, r3, erg_3 = solve_radial(mu, M, a, n3, l2, m3; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe)
    # itp = LinearInterpolation(log10.(rl .* (GNew .* M)), r3 ./ (GNew .* M).^(3/2), extrapolation_bc=Line())
    # rf_3 = itp(log10.(rlist))
    rf_3 = radial_bound_NR(n3, l3, m3, mu, M, rlist; physU=false)
    val3 = sr_rates(n3, l3, m3,  mu, M, a, impose_low_cut=1e-10) * (GNew .* M)
    erg_3 = erg_3G * GNew * M .+ val3 .* im


    # wR1, wI1 = find_im_part(mu, M, a, n1, l1, m1; debug=false, return_both=true, for_s_rates=true, Ntot_force=10000)
    # wR2, wI2 = find_im_part(mu, M, a, n2, l2, m2; debug=false, return_both=true, for_s_rates=true, Ntot_force=10000)
    # wR3, wI3 = find_im_part(mu, M, a, n3, l3, m3; debug=false, return_both=true, for_s_rates=true, Ntot_force=10000)
    
    # erg = ((wR1 .+ wR2 .- wR3) .+ im * (wI1 .+ wI2 .- wI3))   #
    # erg = ergL(1, 0, 0, mu, M, a; full=false) * GNew * M .+ 0 * im ## TEST!!!!
    erg = (erg_1 + erg_2 - erg_3)  #
    print("ERG \t", Float64.(real(erg)), "\t [Im]", Float64.(imag(erg)), "\n")

    Z1 = spheroidals(l1, m1, a, erg_1)
    Z2 = spheroidals(l2, m2, a, erg_2)
    Z3 = spheroidals(l3, m3, a, erg_3)
    Z4 = spheroidals(l, m, a, erg)
   
    thetaV = acos.(1.0 .- 2.0 .* rand(Nang))
    phiV = rand(Nang) .* 2*pi

    function func_ang(x)
        return real(Z1.(x[1], x[2]) .* Z2.(x[1], x[2]) .* conj(Z3.(x[1], x[2])) .* conj(Z4.(x[1], x[2])))
    end
    function func_ang_2(x)
        return real(Z1.(x[1], x[2]) .* Z2.(x[1], x[2]) .* conj(Z3.(x[1], x[2])) .* conj(Z4.(x[1], x[2]))) .* a.^2 .* cos.(x[1]).^2
    end
    CG = 0.0
    CG_2 = 0.0
    for i in 1:Nang
        CG += func_ang([thetaV[i], phiV[i]])
        CG_2 += func_ang_2([thetaV[i], phiV[i]])
    end
    CG *= 4*pi / Nang
    CG_2 *= 4*pi / Nang
    
    print("CG / CG2 \t", CG, "\t", CG_2, "\n")
    print("rp  ", Float64.(rp), "\t alpha \t", alph, "\n")
    
    lam_eff = 1.0e0
    if (n1==n2)&&(l1==l2)&&(m1==m2)
        preFac = 3
    else
        preFac = 6
    end
    unitMatch = -1.0 ./ (2 .* alph).^(3/2) .* lam_eff
    gammaT = (preFac .* (rf_1 .* rf_2 .* conj(rf_3)) ./ 6.0 .* unitMatch)
    print("Max test \t", Float64.(maximum(real(abs.(gammaT)))), "\t", Float64.(rlist[argmax(real(abs.(gammaT)))]), "\n")
    print("r1/r2/r3 \t", Float64(maximum(real(abs.(rf_1)))), "\t", Float64(maximum(real(abs.(rf_2)))), "\t", Float64(maximum(real(abs.(rf_3)))), "\n")
    
    itpG = LinearInterpolation(log10.(rlist), Float64.(real.(gammaT)), extrapolation_bc=Line())
    itpGI = LinearInterpolation(log10.(rlist), Float64.(imag.(gammaT)), extrapolation_bc=Line())
    
    # writedlm("test_store/source_test.dat", hcat(rlist, Float64.(real.(gammaT))))
    # writedlm("test_store/source_test_i.dat", hcat(rlist, Float64.(imag.(gammaT))))

    Mvars = [1.0]
    
    if sve_for_test
        nn = l + 1
        wR, wI = find_im_part(mu, M, a, nn, l, m; debug=false, iter=500, xtol=1e-20, ftol=1e-90, return_both=true, for_s_rates=true)
        erg = wR .+ im .* wI
    end
    
    
    
    gam = im * a * sqrt.(erg.^2 .- alph.^2)
    LLM = l * (l + 1)
    LLM += (-1 + 2 * l * (l + 1) - 2 * m.^2) * gam.^2 ./ (-3 + 4 * l * (l + 1))
    LLM += ((l - m - 1 * (l - m) * (l + m) * (l + m - 1)) ./ ((-3 + 2 * l) * (2 * l - 1).^2) - (l + 1 - m) * (2 * l - m) * (l + m + 1) * (2 + l + m) ./ ((3 + 2 * l).^2 * (5 + 2 * l))) * gam.^4 ./ (2 * (1 + 2 * l))
    LLM += (4 * ((-1 + 4 * m^2) * (l * (1 + l) * (121 + l * (1 + l) * (213 + 8 * l * (1 + l) * (-37 + 10 * l * (1 + l)))) - 2 * l * (1 + l) * (-137 + 56 * l * (1 + l) * (3 + 2 * l * (1 + l))) * m^2 + (705 + 8 * l * (1 + l) * (125 + 18 * l * (1 + l))) * m^4 - 15 * (1 + 46 * m^2))) * gam^6) / ((-5 + 2 * l) * (-3 + 2 * l) * (5 + 2 * l) * (7 + 2 * l) * (-3 + 4 * l * (1 + l))^5)
    
    
    
    r_list_map = 10 .^LinRange(log10.(rp * (1.0 .+ eps_fac)), log10.(rmax), 100000)
    rout_star = r_list_map .+ 2.0 .* rp ./ (rp .- rmm) .* log.(r_list_map ./ rp .- 1.0) .- 2.0 .* rmm .* log.(r_list_map ./ rmm .- 1.0) ./ (rp .- rmm)
    itp_rrstar = LinearInterpolation(rout_star, r_list_map, extrapolation_bc=Line())
    
    
    

    h_step = rp ./ h_mve
    if add_source
        rstart = 0.0 .* (1 + im)
    else
        rstart = 1e-70 .* (1 + im)
    end
    outWF = []
    rvals = []
    SS_out = []
    
    rr = itp_rrstar(rmax)


    append!(outWF, rstart)
    append!(rvals, rr)
    append!(SS_out, 0.0)
    rr -= h_step
    append!(outWF, rstart)
    append!(rvals, rr)
    append!(SS_out, 0.0)
    rr -= h_step
    
    HigherO = true
    run_it = true
    idx = 2
    while run_it
        r_input = itp_rrstar(rr)
        h1 = rvals[idx] - rvals[idx - 1]
        h2 = rr - rvals[idx]
        SS = (itpG(log10.(r_input)) + im * itpGI(log10.(r_input)))
        
        if NON_REL
            newV = outWF[idx] .* (1 .+ h1 ./ h2) .- outWF[idx - 1] .* (h1 ./ h2) .+ (h1.^2 .+ h1 .* h2) .* ((alph.^2 .- erg.^2) .* outWF[idx] .- 2 .* alph.^2 ./ r_input .* outWF[idx]) ./ 2.0
            
            append!(SS_out, -SS .* CG)
            if add_source
                newV += - (h1.^2 .+ h1 .* h2) .* SS .* (CG .* (r_input .- 2 .- a.^2 ./ (2 .* r_input)) .+ CG_2 .* a.^2 ./ r_input) ./ 2.0
                
            end
            
            if HigherO
                newV += (h1.^2 .+ h1 .* h2) .* outWF[idx] .* (a.^2 .* (- alph.^2 .+ erg.^2)) ./ r_input.^2 ./ 2.0
                if add_source
                    newV += - (h1.^2 .+ h1 .* h2) .* SS .* (CG .* 3 ./ r_input.^2 .- 2 .* CG_2 ./ r_input.^2 ) ./ 2.0
                end
                
                newV += (h1.^2 .+ h1 .* h2) .* outWF[idx] .* (2 .+ a.^2 .* (4 .* alph.^2 .- 2 .* erg.^2)) ./ r_input.^3 ./ 2.0
                newV += (h1.^2 .+ h1 .* h2) .* outWF[idx] .* (-4 .+ a.^2 .+ a.^4 .* (alph.^2 .- erg.^2)) ./ r_input.^4 ./ 2.0
              
            end
            newV_r = Float64.(real(newV))
            newV_i = Float64.(imag(newV))
            append!(outWF, newV_r + im * newV_i)
        else
            
            delt = (r_input.^2 .- 2 .* r_input .+ a.^2)
            ff = (r_input.^2 .+ a.^2)
            Vv = delt .* alph.^2 ./ ff .+ delt .* (LLM .+ a.^2 .* (erg.^2 .- alph.^2)) ./ ff.^2 .+ delt .* (3 .* r_input.^2 .- 4 .* r_input .+ a.^2) ./ ff.^3 .- 3 .* delt.^2 .* r_input.^2 ./ ff.^4
            newV = outWF[idx] .* (1 .+ h1 ./ h2) .- outWF[idx - 1] .* (h1 ./ h2) .+ (h1.^2 .+ h1 .* h2) .* (Vv .- erg.^2) .* outWF[idx] ./ 2.0
            if add_source
                newV += - (h1.^2 .+ h1 .* h2) .* SS .* (CG .* r_input.^2 .+ CG_2 .* a.^2) .* delt ./ ff.^(3/2) ./ 2.0
                
            end
            
            newV_r = Float64.(real(newV))
            newV_i = Float64.(imag(newV))
            
            
            SS_r = Float64.(real(-SS .* CG))
            SS_i = Float64.(real(-SS .* CG))
            append!(SS_out, SS_r + im * SS_i)
            append!(outWF, newV_r + im * newV_i)
        end
        
        append!(rvals, rr)

        rr -= h_step
        idx += 1
        
        if (itp_rrstar(rr) < rp .*  (1.0 .+ eps_fac))
            run_it = false
        end
        
        
        
    end
    
    rvals = reverse(rvals)
    outWF = reverse(outWF)
    SS_out = reverse(SS_out)
    outWF .*= 1.0 ./ (sqrt.(itp_rrstar(rvals).^2 .+ a.^2)  .* (GNew .* M) )
    
    ### TESTING NR
    testit = false
    sv1 = []
    sv2 = []
    rt = []
    ss1 = []
    ss2 = []
    if testit
        uu = outWF .* (sqrt.(itp_rrstar(rvals).^2 .+ a.^2)  .* (GNew .* M) )
        for i in 2:(length(rvals) - 1)
            r_input = itp_rrstar(rvals[i])
            ## first test diff eq im solving...
            d2u = (uu[i+1] .- 2 .* uu[i] .+ uu[i-1]) ./ h_step.^2
            potV = (alph.^2 .- erg.^2) .* uu[i] .- 2 .* alph.^2 .* uu[i] ./ r_input
            
            
            ## now do variable transform check
            potV2 = (alph.^2 .- erg.^2 .- 2 .* alph.^2 ./ r_input) .* outWF[i] .* (GNew .* M)
            secD = - (outWF[i+1] .- 2 .* outWF[i] .+ outWF[i-1]) .* (GNew .* M) ./ h_step.^2
            firD = - 2 .* (outWF[i+1] .- outWF[i-1]) .* (GNew .* M) ./ (2 .* h_step .* r_input)
            
            rpt = Float64(r_input)
            ot1 = Float64(real((d2u .- potV )))
            ot2 = Float64(real((secD .+ firD .+ potV2) ))
            # print("Check: \t", rpt, "\t", ot1, "\t", ot2, "\n")
            append!(rt, rpt)
            append!(sv1, ot1)
            append!(sv2, ot2)
            
            append!(ss1, SS_out[i] .* (r_input .- 2))
            append!(ss2, SS_out[i])
            
        end
        writedlm("test_store/check1.dat", hcat(rt, sv1))
        writedlm("test_store/check2.dat", hcat(rt, sv2))
        
        writedlm("test_store/ss_check1.dat", hcat(rt, Float64.(real(ss1))))
        writedlm("test_store/ss_check2.dat", hcat(rt, Float64.(real(ss2))))
    end
    ###
    
    

    maxV = maximum(real(outWF .* conj.(outWF)))
    
    # lam = (mu ./ (M_pl .* 1e9))^2
    # rate_out = 4 .* alph.^2 .* (1 .+ sqrt.(1 - a.^2)) .* (maxV ./ lam_eff.^2) .* lam^2 ./ mu^2 .* (GNew * M^2 * M_to_eV)^2
    out = Float64.(real(8 * pi .* maxV .* erg.^2 ./ mu.^2 .* (1 .+ sqrt.(1 .- a.^2)) ./ lam_eff.^2))
        
    
    trapz(y,x) = @views sum(((y[1:end-1].+y[2:end])/2).*(x[2:end].-x[1:end-1]))
    rreal = itp_rrstar(rvals)
    nm = trapz(real(outWF .* conj.(outWF)) .* rreal.^2, rreal)
    if NON_REL && sve_for_test
        writedlm("test_store/test_Out_NR.dat", hcat(rreal, real(outWF .* conj.(outWF)) ./ nm))
    elseif !NON_REL && sve_for_test
        writedlm("test_store/test_Out.dat", hcat(rreal, real(outWF .* conj.(outWF)) ./ nm))
        
        nn = l + 1
        rr_test, outtest, erg_3 = solve_radial(mu, M, a, nn, l, m; rpts=Npts_Bnd, rmaxT=rmaxT, return_erg=true, Ntot_safe=Ntot_safe, eps_r=eps_fac)
        rNR = radial_bound_NR(nn, l, m, mu, M, rr_test)
        writedlm("test_store/test_NR.dat", hcat(rr_test, rNR.^2))
        writedlm("test_store/test_R.dat", hcat(rr_test, real(outtest .* conj.(outtest))))
    end
    
    if NON_REL && add_source
        writedlm("test_store/test_RAD_NR_v2.dat", hcat(itp_rrstar(rvals), Float64.(real(outWF .* conj.(outWF)))))
        if testit
            writedlm("test_store/scs_test.dat", hcat(itp_rrstar(rvals), Float64.(real(SS_out))))
        end
    elseif add_source
        writedlm("test_store/test_RAD_v2.dat", hcat(itp_rrstar(rvals), Float64.(real(outWF .* conj.(outWF)))))
    end
    

    print("Output: ", out, "\n")
    # print("Output [2]: ", rate_out, "\n")
    
    # theory-NR
    test = 4.3e-7 * alph.^(11 - 4) .* (1 .+ sqrt.(1 .- a.^2))
    print("NR Rate: ", test, "\n")
    # test2 = 4.3e-7 * alph.^(11) .* (1 .+ sqrt.(1 .- a.^2))
    # print("NR Rate [2]: ", test2, "\n")
    print("Ratio : ", out ./ test, "\n")
    # print("Ratio [2] : ", rate_out ./ test2, "\n")

    
    return out
    
end
