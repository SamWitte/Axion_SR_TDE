
function load_rate_coeffs(mu, M, a, f_a; non_rel=true, input_data="Me", solve_n4=true, solve_n5=false)
    alph = mu * GNew * M
    rP = 1 + sqrt.(1 - a^2)
    faFac = (M_pl ./ f_a)^4
    
    Drate = Dict()
    
    Drate["211_322^GW"] = 0.0 * alph^16
    
    if input_data == "Doddy"

        Drate["211_211^322^BH"] = 0.0 # format gamma_{A x B}^{C x D}
        Drate["322_322^211^Inf"] = 0.0
        Drate["211_211^GW"] = 0.0
        Drate["322_211^GW"] = 0.0
        Drate["211_211_211^Inf"] = 0.0
        Drate["322_322^GW"] = 0.0
        
        # n = 4
        if solve_n4
            Drate["211_411^322^BH"] = 0.0
            Drate["411_411^322^BH"] = 0.0
            Drate["322_411^211^Inf"] = 0.0
            Drate["322_411^433^BH"] = 0.0
            Drate["211_211^422^BH"] = 0.0
            Drate["411_422^211^Inf"] = 0.0
            Drate["411_411^422^BH"] = 0.0
            Drate["211_422^433^BH"] = 0.0
            Drate["422_322^211^Inf"] = 0.0
            Drate["433_433^211^Inf"] = 0.0
            Drate["322_433^211^Inf"] = 0.0
            Drate["211_322^433^BH"] = 0.0
            
            if solve_n5
                Drate["322_322^544^BH"] = 0.0
                Drate["322_411^533^BH"] = 0.0
                Drate["322_422^544^BH"] = 0.0
                Drate["411_411^522^BH"] = 0.0
                Drate["322_522^544^BH"] = 0.0
                Drate["422_422^544^BH"] = 0.0
                Drate["422_522^544^BH"] = 0.0
                Drate["522_522^544^BH"] = 0.0
                
                Drate["422_544^322^Inf"] = 0.0
                Drate["433_544^322^Inf"] = 0.0
                Drate["422_533^322^Inf"] = 0.0
                Drate["433_533^322^Inf"] = 0.0
                Drate["433_522^322^Inf"] = 0.0
                Drate["422_522^322^Inf"] = 0.0
            end
        end
            
            
        if non_rel
            Drate["211_211^322^BH"] = 4.2e-7 .* alph^11 .* faFac * rP
            Drate["322_322^211^Inf"] = 1.1e-8 * alph^8 .* faFac
            Drate["211_211^GW"] = 1.0e-2 * alph^14
            Drate["322_211^GW"] = 5.0e-6 * alph^10
            Drate["211_211_211^Inf"] = 1.5e-8 * alph^21 .* faFac
            Drate["322_322^GW"] = 3.0e-8 * alph^18
            
            # n = 4
            if solve_n4
                Drate["211_411^322^BH"] = 2.5e-8 * alph^11 * faFac * rP
                Drate["411_411^322^BH"] = 1.7e-11 * alph^11 * faFac * rP ### Disagree
                Drate["322_411^211^Inf"] = 3.7e-9 * alph^8 * faFac
                Drate["211_211^422^BH"] = 1.5e-7 * alph^11 * faFac * rP
                Drate["411_422^211^Inf"] = 2.2e-9 * alph^8 * faFac
                Drate["411_411^422^BH"] = 2.2e-11 * alph^7 * faFac * rP
                Drate["211_422^433^BH"] = 7.83e-11 * alph^7 * faFac * rP ### Disagree
                Drate["422_322^211^Inf"] = 1.6e-8 * alph^8 * faFac
                Drate["433_433^211^Inf"] = 9.2e-11 * alph^8 * faFac
                Drate["322_433^211^Inf"] = 2.6e-9 * alph^8 * faFac
                Drate["211_322^433^BH"] = 9.1e-8 * alph^11 * faFac * rP
                Drate["322_411^433^BH"] = 3.8e-11 * alph^7 * faFac * rP
                
                if solve_n5
                    Drate["322_322^544^BH"] = 1.9e-9 * alph^11 * faFac * rP
                    Drate["322_411^533^BH"] = 2.0e-8 * alph^11 * faFac * rP
                    Drate["322_422^544^BH"] = 1.3e-11 * alph^11 * faFac * rP #
                    Drate["411_411^522^BH"] = 9.0e-11 * alph^11 * faFac * rP #
                    Drate["322_522^544^BH"] = 3.4e-12 * alph^7 * faFac * rP
                    Drate["422_422^544^BH"] = 2.3e-9 * alph^11 * faFac * rP #
                    Drate["422_522^544^BH"] = 3.7e-14 * alph^7 * faFac * rP
                    Drate["522_522^544^BH"] = 2.7e-13 * alph^7 * faFac * rP #
                    
                    Drate["422_544^322^Inf"] = 2.5e-11 * alph^8 * faFac
                    Drate["433_544^322^Inf"] = 7.8e-10 * alph^8 * faFac
                    Drate["422_533^322^Inf"] = 1.2e-9 * alph^8 * faFac
                    Drate["433_533^322^Inf"] = 2.8e-9 * alph^8 * faFac
                    Drate["433_522^322^Inf"] = 6.3e-10 * alph^8 * faFac
                    Drate["422_522^322^Inf"] = 1.6e-9 * alph^8 * faFac
                    
                    Drate["522_544^322^Inf"] = 2.2e-11 * alph^8 * faFac
                    Drate["533_544^322^Inf"] = 1.8e-10 * alph^8 * faFac
                    Drate["544_544^322^Inf"] = 4.3e-11 * alph^8 * faFac
                    Drate["522_533^322^Inf"] = 4.4e-10 * alph^8 * faFac
                    Drate["533_533^322^Inf"] = 3.1e-10 * alph^8 * faFac
                    Drate["522_522^322^Inf"] = 1.6e-10 * alph^8 * faFac
                end
            end
        end
            
    else
        Drate["211_211^GW"] = 1.0e-2 * alph^14
        Drate["322_211^GW"] = 5.0e-6 * alph^10
        Drate["211_211_211^Inf"] = 1.5e-8 * alph^21 .* faFac
        Drate["322_322^GW"] = 3.0e-8 * alph^18
        
        
        # rates computed for fixed rP(a=0.9)
        rP_ratio = rP / (1 + sqrt.(1.0 - 0.9^2))
        
        # load files and get rates
        data = open(readdlm, "rate_sve/211_211_322_BH_NR.dat")
        itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
        Drate["211_211^322^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
        
        data = open(readdlm, "rate_sve/322_322_211_Inf_NR.dat")
        itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
        Drate["322_322^211^Inf"] = 10 .^itp(log10.(alph)) .* faFac
        
        if solve_n4
            data = open(readdlm, "rate_sve/211_411_322_BH_NR.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["211_411^322^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
            
            
            data = open(readdlm, "rate_sve/411_411_322_BH_NR.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["411_411^322^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
        
            
            data = open(readdlm, "rate_sve/322_411_211_Inf_NR.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["322_411^211^Inf"] = 10 .^itp(log10.(alph)) .* faFac
            
            data = open(readdlm, "rate_sve/211_211_422_BH_NR.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["211_211^422^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
            
            data = open(readdlm, "rate_sve/411_422_211_Inf_NR.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["411_422^211^Inf"] = 10 .^itp(log10.(alph)) .* faFac
            
            data = open(readdlm, "rate_sve/411_411_422_BH_NR.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["411_411^422^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
            
            data = open(readdlm, "rate_sve/211_422_433_BH_NR.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["211_422^433^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
            
            data = open(readdlm, "rate_sve/422_322_211_Inf_NR.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["422_322^211^Inf"] = 10 .^itp(log10.(alph)) .* faFac
            
            data = open(readdlm, "rate_sve/433_433_211_Inf_NR.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["433_433^211^Inf"] = 10 .^itp(log10.(alph)) .* faFac
            
            data = open(readdlm, "rate_sve/322_433_211_Inf_NR.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["322_433^211^Inf"] = 10 .^itp(log10.(alph)) .* faFac
            
            data = open(readdlm, "rate_sve/211_322_433_BH_NR.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["211_322^433^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
            
            data = open(readdlm, "rate_sve/322_411_433_BH_NR.dat")
            itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
            Drate["322_411^433^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac

            
            if solve_n5
                data = open(readdlm, "rate_sve/322_322_544_BH_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["322_322^544^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
                
                data = open(readdlm, "rate_sve/411_322_533_BH_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["322_411^533^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
                
                data = open(readdlm, "rate_sve/422_322_544_BH_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["322_422^544^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
            
                data = open(readdlm, "rate_sve/411_411_522_BH_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["411_411^522^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
                
                data = open(readdlm, "rate_sve/522_322_544_BH_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["322_522^544^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
                
                data = open(readdlm, "rate_sve/422_422_544_BH_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["422_422^544^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
                
                data = open(readdlm, "rate_sve/522_422_544_BH_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["422_422^544^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
                
                data = open(readdlm, "rate_sve/522_522_544_BH_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["522_522^544^BH"] = 10 .^itp(log10.(alph)) .* rP_ratio .* faFac
                    
                data = open(readdlm, "rate_sve/422_544_322_Inf_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["422_544^322^Inf"] = 10 .^itp(log10.(alph))  .* faFac

                data = open(readdlm, "rate_sve/433_544_322_Inf_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["433_544^322^Inf"] = 10 .^itp(log10.(alph))  .* faFac
                
                data = open(readdlm, "rate_sve/422_533_322_Inf_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["422_533^322^Inf"] = 10 .^itp(log10.(alph))  .* faFac
                
                data = open(readdlm, "rate_sve/433_533_322_Inf_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["433_533^322^Inf"] = 10 .^itp(log10.(alph))  .* faFac
                
                data = open(readdlm, "rate_sve/433_522_322_Inf_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["433_522^322^Inf"] = 10 .^itp(log10.(alph))  .* faFac
                
                data = open(readdlm, "rate_sve/422_522_322_Inf_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["422_522^322^Inf"] = 10 .^itp(log10.(alph))  .* faFac
                
                
                
                data = open(readdlm, "rate_sve/522_544_322_Inf_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["522_544^322^Inf"] = 10 .^itp(log10.(alph))  .* faFac
                
                data = open(readdlm, "rate_sve/533_544_322_Inf_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["533_544^322^Inf"] = 10 .^itp(log10.(alph))  .* faFac


                data = open(readdlm, "rate_sve/544_544_322_Inf_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["544_544^322^Inf"] = 10 .^itp(log10.(alph))  .* faFac
                
                
                data = open(readdlm, "rate_sve/522_533_322_Inf_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["522_533^322^Inf"] = 10 .^itp(log10.(alph))  .* faFac
                
                data = open(readdlm, "rate_sve/533_533_322_Inf_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["533_533^322^Inf"] = 10 .^itp(log10.(alph))  .* faFac
                
                data = open(readdlm, "rate_sve/522_522_322_Inf_NR.dat")
                itp = LinearInterpolation(log10.(data[:, 1]), log10.(data[:, 2]), extrapolation_bc=Line())
                Drate["522_522^322^Inf"] = 10 .^itp(log10.(alph))  .* faFac
                     
                
            end
        end
    end

    return Drate
end

function key_to_indx(keyN; solve_n4=false, solve_n5=false)
    if keyN[9:10] == "GW"
        totN = 3
    else
        totN = 4
    end
    outPix = zeros(Int, totN)
    sgn = zeros(totN)
    for i in 1:totN
        if i == 1
            nme = keyN[1:3]
            sgn[i] = -1.0
        elseif i == 2
            nme = keyN[5:7]
            sgn[i] = -1.0
        elseif i == 3
            if totN == 4
                nme = keyN[9:11]
                if keyN[8] == "_"
                    sgn[i] = -1.0
                else
                    sgn[i] = 1.0
                end
            else
                nme = "GW"
                sgn[i] = 1.0
            end
        else
            nme = keyN[13:end]
            sgn[i] = 1.0
        end

        
            
        if solve_n4
            if nme == "211"
                outPix[i] = 1
            elseif nme == "322"
                outPix[i] = 2
            elseif nme == "411"
                outPix[i] = 3
            elseif nme == "422"
                outPix[i] = 4
            elseif nme == "433"
                outPix[i] = 5
            elseif nme == "BH"
                outPix[i] = -1
            elseif (nme == "GW")||(nme == "Inf")
                outPix[i] = 0
            else
                if solve_n5
                    if nme == "522"
                        outPix[i] = 6
                    elseif nme == "533"
                        outPix[i] = 7
                    elseif nme == "544"
                        outPix[i] = 8
                    else
                        print("WHAT!? \t ", nme, "\n")
                    end
                else
                    print("WHAT!? \t ", nme, "\n")
                end
            end
        else
            if nme == "211"
                outPix[i] = 1
            elseif nme == "322"
                outPix[i] = 2
            elseif nme == "BH"
                outPix[i] = -1
            elseif (nme == "GW")||(nme == "Inf")
                outPix[i] = 0
            else
                print("WHAT!? \t ", nme, "\n")
            end
        end
    end
    return outPix, sgn
end
