using JLD2, Statistics, LinearAlgebra
using Dierckx, Polynomials
using Plots, LaTeXStrings, Measures

Pb_centrality = load("Pb_centrality_0_30.jld2")["Nch"]
Xe60_centrality = load("Xe60_centrality_0_30.jld2")["Nch"]
Xe27_centrality = load("Xe_centrality_0_30.jld2")["Nch"]
Xe0_centrality = load("Xe0_centrality_0_30.jld2")["Nch"]
Pb_allevents = load("triaxial/Pb_epsilon2345_2023-05-05.jld2")["data"]
Xe60_allevents = load("triaxial/Xe60_epsilon2345_2023-05-06.jld2")["data"]
Xe27_allevents = load("triaxial/Xe_epsilon2345_2023-05-04.jld2")["data"]
Xe0_allevents = load("triaxial/Xe0_epsilon2345_2023-05-06.jld2")["data"]
# Xe60_allevents = load("triaxial/gamma60_2023-05-18.jld2")["data"]
# Xe27_allevents = load("triaxial/gamma27_2023-05-18.jld2")["data"]
# Xe0_allevents = load("triaxial/gamma0_2023-05-18.jld2")["data"]
Pb_Nch = [event[1] for event in Pb_allevents]
Xe60_Nch = [event[1] for event in Xe60_allevents]
Xe27_Nch = [event[1] for event in Xe27_allevents]
Xe0_Nch = [event[1] for event in Xe0_allevents]
events = length(Pb_allevents)

clist = 0.5:1:24.5

function data2dot(centrality, allevents, Nch, clist)
  index = [findfirst(x -> x < centrality[i], Nch) for i in 1:length(clist)]
  data = Vector{Tuple{Int64,Float64,Vararg{ComplexF64,4}}}[]
  # data = Vector{Tuple{Int64,Int64}}[]
  # N = zeros(length(clist))
  for i in 1:length(clist)
    st = (i == 1 ? 1 : index[i-1])
    ed = (i == length(clist) ? length(Nch) : index[i] - 1)
    push!(data, allevents[st:ed])
    # N[i] = mean(Nch[st:ed])
  end
  return data#, N
end

Pb_data = data2dot(Pb_centrality, Pb_allevents, Pb_Nch, clist)
Xe60_data = data2dot(Xe60_centrality, Xe60_allevents, Xe60_Nch, clist)
Xe27_data = data2dot(Xe27_centrality, Xe27_allevents, Xe27_Nch, clist)
Xe0_data = data2dot(Xe0_centrality, Xe0_allevents, Xe0_Nch, clist)

# ## binary collisions
# function calbin(data, clist)
#   bin = zeros(length(clist))
#   for i = 1:length(clist)
#     event = data[i]
#     bin[i] = mean([a[2] for a in event])
#   end
#   return bin
# end

# Xe60_bin = calbin(Xe60_data, clist)
# Xe27_bin = calbin(Xe27_data, clist)
# Xe0_bin = calbin(Xe0_data, clist)

# ## Nch
# N1 = NXe1 ./ NPb
# N2 = NXe2 ./ NPb
# N3 = NXe3 ./ NPb

p0 = plot(size=(800, 600), minorgrid=true, margins=3mm)
plot!(clist, Xe60_bin, lc=2, lw=4, ls=:dot, label="Xe+Xe " * L"\beta=0.207,\gamma=60^\circ")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
plot!(clist, Xe27_bin, lc=3, lw=4, ls=:dash, label="Xe+Xe " * L"\beta=0.207,\gamma=27^\circ")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
plot!(clist, Xe0_bin, lc=4, lw=4, ls=:dashdot, label="Xe+Xe " * L"\beta=0.207,\gamma=0^\circ")
xlims!(0, 25)
ylims!(500, 1800)
xlabel!("centrality(%)")
ylabel!(L"N_{\mathrm{q-coll}}")
annotate!(8, 800, text("Number of binary collisions", 20, "Computer Modern"))
annotate!(8, 650, text("vs " * L"N_{\mathrm{ch}}" * " based centrality ", 20, "Computer Modern"))
# savefig(p0, "figures/beta/Nch_Xe_Pb.svg")

## 计算双关联系数
function rho(a, b)
  δa = a .- mean(a)
  δb = b .- mean(b)
  ρ = mean(δa .* δb) / sqrt(mean(δa .^ 2) * mean(δb .^ 2))
  return ρ
end
function rholist(data, clist)
  ρlist = zeros(length(clist))
  for i in 1:length(clist)
    event = data[i]
    # sort!(event, by=n -> -n[1])
    p̄T = [y[2] for y in event]
    # ε2 = [abs(y[3]) for y in event]
    # ε3 = [abs(y[4]) for y in event]
    # ε4 = [abs(y[5]) for y in event]
    # ε5 = [abs(y[6]) for y in event]
    Φ2 = [angle(-y[3]) / 2 for y in event]
    # Φ3 = [angle(-y[4]) / 3 for y in event]
    # Φ4 = [angle(-y[5]) / 4 for y in event]
    # Φ5 = [angle(-y[6]) / 5 for y in event]
    # observables = [p̄T, ε2, ε3, ε4, ε5, Φ2, Φ3, Φ4, Φ5]
    # ρlist[i] = rho(ε4 .^ 2, Φ4 .^2)
    ρlist[i] = rho(cos.(2*Φ2) ,p̄T)
    # ρlist[i] = mean(ε2 .^2)
  end
  return ρlist
end

Pb_ρlist = rholist(Pb_data, clist)
Xe60_ρlist = rholist(Xe60_data, clist)
Xe27_ρlist = rholist(Xe27_data, clist)
Xe0_ρlist = rholist(Xe0_data, clist)
# Xe1_ρlist = rholist(Xe1_data, clist)
# Xe2_ρlist = rholist(Xe2_data, clist)
# Xe3_ρlist = rholist(Xe3_data, clist)

ρ60list = Xe60_ρlist ./ Pb_ρlist
ρ27list = Xe27_ρlist ./ Pb_ρlist
ρ0list = Xe0_ρlist ./ Pb_ρlist

# p0 = plot(size=(800, 600), minorgrid=true, margins=3mm)
# plot!(clist, Xe60_ρlist./Pb_ρlist, lc=2, lw=4, ls=:dot, label="Xe+Xe " * L"\beta=0.207,\gamma=60^\circ")
# scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
# plot!(clist, Xe27_ρlist./Pb_ρlist, lc=3, lw=4, ls=:dash, label="Xe+Xe " * L"\beta=0.207,\gamma=27^\circ")
# scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
# plot!(clist, Xe0_ρlist./Pb_ρlist, lc=4, lw=4, ls=:dashdot, label="Xe+Xe " * L"\beta=0.207,\gamma=0^\circ")
# xlims!(0, 25)
# ylims!(1, 2.2)
# xlabel!("centrality(%)")
# ylabel!("Xe+Xe/Pb+Pb")
# annotate!(16, 1.7, text(L"\epsilon_2^2"*" ratio vs "*L"N_{\mathrm{ch}}"*" based centrality ", 20, "Computer Modern"))
# savefig(p0, "figures/epsilon2_Xe_Pb.svg")

## 计算三关联系数
function ac(a, b, c)
  numerator = mean(abs.(a) .* abs.(b) .* abs.(c) .* cos.(angle.(-a) + angle.(-b) - angle.(-c)))
  denominator = sqrt(mean((abs.(a)).^2 .* (abs.(b)).^2) * mean(abs.(c).^2))
  return numerator / denominator
end
function aclist(data, clist)
  aclist = zeros(length(clist))
  for i in 1:length(clist)
    event = data[i]
    sort!(event, by=n -> -n[1])
    ϵ2 = [y[3] for y in event]
    # ϵ3 = [y[4] for y in event]
    ϵ4 = [y[5] for y in event]
    # ϵ5 = [y[6] for y in event]
    aclist[i] = ac(ϵ2, ϵ2, ϵ4)
  end
  return aclist
end

Pb_aclist = aclist(Pb_data, clist)
Xe60_aclist = aclist(Xe60_data, clist)
Xe27_aclist = aclist(Xe27_data, clist)
Xe0_aclist = aclist(Xe0_data, clist)

ac60list = Xe60_aclist ./ Pb_aclist
ac27list = Xe27_aclist ./ Pb_aclist
ac0list = Xe0_aclist ./ Pb_aclist

###############################################################
## ρ2
Pb_ATLAS_cen = [0.065036, 0.15521, 0.28418, 0.50332, 0.80539, 1.2389, 1.9339, 2.9556, 4.1093, 5.5531,
                7.3742, 9.3482, 11.525, 13.883, 15.821, 17.174, 18.606, 20.08, 21.646, 23.288, 25.014]
Pb_ATLAS_rho = [0.14228, 0.14366, 0.1632, 0.18992, 0.22461, 0.26368, 0.2905, 0.3218, 0.32663, 0.32589,
                0.32525, 0.32766, 0.32568, 0.32291, 0.31921, 0.31288, 0.31346, 0.31335, 0.30601, 0.30703, 0.30292]
Xe_ATLAS_cen = [0.19408, 1.2063, 2.2081, 3.6334, 5.2434, 7.0539, 9.0113, 11.087, 13.339, 15.758, 18.429, 21.232, 24.397]
Xe_ATLAS_rho = [0.13656, 0.15903, 0.19554, 0.22659, 0.23505, 0.24894, 0.25364, 0.25053, 0.24492, 0.23603, 0.22798, 0.2354, 0.22228]

# 插值后才能画比值
itp_Pb = fit(Polynomial, Pb_ATLAS_cen, Pb_ATLAS_rho, 8)
itp_Xe = fit(Polynomial, Xe_ATLAS_cen, Xe_ATLAS_rho, 8)

## Plot
using Plots, LaTeXStrings, Measures

p1 = plot(size=(800, 600), minorgrid=true, margins=3mm)
scatter!(Pb_ATLAS_cen, Pb_ATLAS_rho, mc=1, ms=9, label = "")#"ATLAS, Pb+Pb " * L"N_{\mathrm{ch}}" * " based")
scatter!(Xe_ATLAS_cen, Xe_ATLAS_rho, markershape=:rect, mc=3, ms=7, label = "")#"ATLAS, Xe+Xe " * L"N_{\mathrm{ch}}" * " based")
plot!(clist, Pb_ρlist, lc=1, lw=4, label="Pb+Pb " * L"\beta=0.062,\gamma=27^\circ")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
plot!(clist, Xe1_ρlist, lc=2, lw=4, ls=:dot, label="Xe+Xe " * L"\beta=0.1,\gamma=60^\circ")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
plot!(clist, Xe2_ρlist, lc=3, lw=4, ls=:dash, label="Xe+Xe " * L"\beta=0.2,\gamma=27^\circ")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
plot!(clist, Xe3_ρlist, lc=4, lw=4, ls=:dashdot, label="Xe+Xe " * L"\beta=0.3,\gamma=0^\circ")
xlims!(0, 25)
ylims!(-0.02, 0.35)
xlabel!("centrality(%)")
ylabel!(L"\rho(\varepsilon_2^2,\langle{p_T}\rangle)")
annotate!(15, 0.275, text(L"\rho_2" * " vs " * L"N_{\mathrm{ch}}" * " based centrality ", 20, "Computer Modern"))
savefig(p1, "figures/rho(epsilon22,pT).svg")

p2 = plot(size=(800, 600), minorgrid=true, margins=3mm)
plot!(clist, ρ1list, lc=2, lw=4, ls=:dot, label="")#L"\beta=0.207,\gamma=60^\circ")
plot!(clist, ρ2list, lc=3, lw=4, ls=:dash, label="")#L"\beta=0.207,\gamma=27^\circ")
plot!(clist, ρ3list, lc=4, lw=4, ls=:dashdot, label="")#L"\beta=0.207,\gamma=0^\circ")
scatter!([30],[1], mc=1, ms=9, label="ATLAS, Pb+Pb at " * L"5.02~\mathrm{TeV}")
scatter!([30],[1], ms=0, mc=:white, msc=:white, label=" ")
scatter!([30],[1], markershape=:rect, mc=3, ms=7, label="ATLAS, Xe+Xe at " * L"5.44~\mathrm{TeV}")
scatter!([30],[1], ms=0, mc=:white, msc=:white, label=" ")
scatter!(LinRange(1, 23, 12), x -> itp_Xe(x) / itp_Pb(x), markershape=:diamond, mc=6, ms=9,
  label = "ATLAS, "*L"\rho_2"*"(Xe+Xe)/"*L"\rho_2"*"(Pb+Pb)")
xlims!(0, 25)
ylims!(0.4, 1.)
xlabel!("centrality(%)")
ylabel!("Xe+Xe/Pb+Pb")
annotate!(15, 0.96, text(L"0.5< p_T <5.0~\mathrm{GeV},~|\eta|<2.5", 20))
annotate!(15, 2, text(L"\epsilon_2^2"*" ratio vs "*L"N_{\mathrm{ch}}"*" based centrality ", 20, "Computer Modern"))
# savefig(p2, "figures/Xe+Xe_Pb+Pb.svg")


######################################################
## ρ3
Pb_cen_3 = [0.032269, 0.2156, 0.37964, 0.65532, 1.0173, 1.4995, 2.4259, 3.8998, 5.5531, 7.3742, 9.3482, 11.525, 13.883, 16.504, 19.329, 22.463, 25.923]
Pb_rho_3 = [0.076142, 0.080688, 0.083653, 0.092948, 0.10197, 0.10186, 0.1034, 0.093274, 0.084449, 0.073862, 0.069925, 0.063607, 0.057259, 0.055398, 0.050814, 0.043765, 0.041603]
Xe_cen_3 = [0.4347, 1.4072, 3.7741, 7.7327, 12.424, 18.025, 24.317]
Xe_rho_3 = [0.079931, 0.10878, 0.10926, 0.087943, 0.064597, 0.049554, 0.042446]

# 插值后才能画比值
itp_Pb = fit(Polynomial, Pb_cen_3, Pb_rho_3, 8)
# itp_Xe = fit(Polynomial, Xe_cen_3, Xe_rho_3, 4)
itp_Xe = Spline1D(Xe_cen_3, Xe_rho_3; k=3, bc="extrapolate")

## Plot
using Plots, LaTeXStrings, Measures

p3 = plot(size=(800, 600), minorgrid=true, margins=3mm)
scatter!(Pb_cen_3, Pb_rho_3, mc=1, ms=9, label="")
scatter!(Xe_cen_3, Xe_rho_3, markershape=:rect, mc=3, ms=7, label="")
plot!(clist, Pb_ρlist, lc=1, lw=4, label="Pb+Pb " * L"\beta=0.062,\gamma=27^\circ")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
plot!(clist, Xe60_ρlist, lc=2, lw=4, ls=:dot, label="Xe+Xe " * L"\beta=0.207,\gamma=60^\circ")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
plot!(clist, Xe27_ρlist, lc=3, lw=4, ls=:dash, label="Xe+Xe " * L"\beta=0.207,\gamma=27^\circ")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
plot!(clist, Xe0_ρlist, lc=4, lw=4, ls=:dashdot, label="Xe+Xe " * L"\beta=0.207,\gamma=0^\circ")
xlims!(0, 25)
ylims!(0.04, 0.12)
xlabel!("centrality(%)")
ylabel!(L"\rho(\varepsilon_3^2,\langle{p_T}\rangle)")
annotate!(8, 0.048, text(L"\rho_3" * " vs " * L"N_{\mathrm{ch}}" * " based centrality ", 20, "Computer Modern"))
# savefig(p3, "figures/rho(epsilon32,pT).svg")

p4 = plot(size=(800, 600), minorgrid=true, margins=3mm)
plot!(clist, ρ60list, lc=2, lw=4, ls=:dot, label="")#L"\beta=0.207,\gamma=60^\circ")
plot!(clist, ρ27list, lc=3, lw=4, ls=:dash, label="")#L"\beta=0.207,\gamma=27^\circ")
plot!(clist, ρ0list, lc=4, lw=4, ls=:dashdot, label="")#L"\beta=0.207,\gamma=0^\circ")
scatter!([30], [1], mc=1, ms=9, label="ATLAS, Pb+Pb at " * L"5.02~\mathrm{TeV}")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
scatter!([30], [1], markershape=:rect, mc=3, ms=7, label="ATLAS, Xe+Xe at " * L"5.44~\mathrm{TeV}")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
scatter!(LinRange(1, 23, 12), x -> itp_Xe(x) / itp_Pb(x), markershape=:diamond, mc=6, ms=9,
  label="ATLAS, " * L"\rho_3" * "(Xe+Xe)/" * L"\rho_3" * "(Pb+Pb)")
xlims!(0, 25)
ylims!(0.6, 1.3)
xlabel!("centrality(%)")
ylabel!("Xe+Xe/Pb+Pb")
annotate!(5, 0.72, text(L"|\eta|<2.5,", 20))
annotate!(5, 0.64, text(L"0.5< p_T <5.0~\mathrm{GeV}", 20))
annotate!(16, 0.86, text(L"\rho_3" * " ratio vs " * L"N_{\mathrm{ch}}" * " based centrality ", 20, "Computer Modern"))
savefig(p4, "figures/Xe+Xe_Pb+Pb_3.svg")

########################################################
## ρ4
Pb_cen_4 = [0.068297, 0.18513, 0.80539, 2.4259, 3.8998, 5.5531, 7.14, 8.6018, 10.14, 11.818, 13.588, 15.5, 17.514, 19.713, 22.041, 24.569]
Pb_rho_4 = [0.044168, 0.070322, 0.082919, 0.093602, 0.087116, 0.081344, 0.069951, 0.067826, 0.067704, 0.074934, 0.070991, 0.079754, 0.097316, 0.1028, 0.11701, 0.12568]
Xe_cen_4 = [0.8232, 2.8789, 6.122, 9.9752, 15.758, 21.979, 26.874]
Xe_rho_4 = [0.0604, 0.092031, 0.086291, 0.069694, 0.09543, 0.11916, 0.1252]

# 插值后才能画比值
# itp_Pb = fit(Polynomial, Pb_cen_4, Pb_rho_4, 8)
# itp_Xe = fit(Polynomial, Xe_cen_4, Xe_rho_4, 4)
itp_Pb = Spline1D(Pb_cen_4, Pb_rho_4; k=3, bc="extrapolate")
itp_Xe = Spline1D(Xe_cen_4, Xe_rho_4; k=3, bc="extrapolate")

## Plot
using Plots, LaTeXStrings, Measures

p5 = plot(size=(800, 600), minorgrid=true, margins=3mm)
scatter!(Pb_cen_4, Pb_rho_4, mc=1, ms=9, label="")
scatter!(Xe_cen_4, Xe_rho_4, markershape=:rect, mc=3, ms=7, label="")
plot!(clist, Pb_ρlist, lc=1, lw=4, label="Pb+Pb " * L"\beta=0.062,\gamma=27^\circ")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
plot!(clist, Xe60_ρlist, lc=2, lw=4, ls=:dot, label="Xe+Xe " * L"\beta=0.207,\gamma=60^\circ")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
plot!(clist, Xe27_ρlist, lc=3, lw=4, ls=:dash, label="Xe+Xe " * L"\beta=0.207,\gamma=27^\circ")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
plot!(clist, Xe0_ρlist, lc=4, lw=4, ls=:dashdot, label="Xe+Xe " * L"\beta=0.207,\gamma=0^\circ")
xlims!(0, 25)
ylims!(0.03, 0.12)
xlabel!("centrality(%)")
ylabel!(L"\rho(\varepsilon_4^2,\langle{p_T}\rangle)")
annotate!(8, 0.04, text(L"\rho_4" * " vs " * L"N_{\mathrm{ch}}" * " based centrality ", 20, "Computer Modern"))
savefig(p5, "figures/rho(epsilon42,pT).svg")

p6 = plot(size=(800, 600), minorgrid=true, margins=3mm)
plot!(clist, ρ60list, lc=2, lw=4, ls=:dot, label="")#L"\beta=0.207,\gamma=60^\circ")
plot!(clist, ρ27list, lc=3, lw=4, ls=:dash, label="")#L"\beta=0.207,\gamma=27^\circ")
plot!(clist, ρ0list, lc=4, lw=4, ls=:dashdot, label="")#L"\beta=0.207,\gamma=0^\circ")
scatter!([30], [1], mc=1, ms=9, label="ATLAS, Pb+Pb at " * L"5.02~\mathrm{TeV}")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
scatter!([30], [1], markershape=:rect, mc=3, ms=7, label="ATLAS, Xe+Xe at " * L"5.44~\mathrm{TeV}")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
scatter!(LinRange(1, 23, 12), x -> itp_Xe(x) / itp_Pb(x), markershape=:diamond, mc=6, ms=9,
  label="ATLAS, " * L"\rho_4" * "(Xe+Xe)/" * L"\rho_4" * "(Pb+Pb)")
xlims!(0, 15)
ylims!(0.6, 1.4)
xlabel!("centrality(%)")
ylabel!("Xe+Xe/Pb+Pb")
annotate!(3.4, 1.3, text(L"0.5< p_T <5.0~\mathrm{GeV},", 20))
annotate!(3.4, 1.22, text(L"|\eta|<2.5", 20))
annotate!(8, 0.7, text(L"\rho_4" * " ratio vs " * L"N_{\mathrm{ch}}" * " based centrality ", 20, "Computer Modern"))
savefig(p6, "figures/Xe+Xe_Pb+Pb_4.svg")

######################################################################\
## ac
using Plots, LaTeXStrings, Measures

# p10 = plot(size=(800, 600), minorgrid=true, margins=3mm)
# plot!(clist, Xe60_aclist, lc=2, lw=4, ls=:dot, label="Xe+Xe " * L"\beta=0.207,\gamma=60^\circ")
# # scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
# plot!(clist, Xe27_aclist, lc=3, lw=4, ls=:dash, label="Xe+Xe " * L"\beta=0.207,\gamma=27^\circ")
# # scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
# plot!(clist, Xe0_aclist, lc=4, lw=4, ls=:dashdot, label="Xe+Xe " * L"\beta=0.207,\gamma=0^\circ")
# xlims!(0, 25)
# ylims!(-0.6, -0.3)
# xlabel!("centrality(%)")
# ylabel!(L"\rho(\epsilon_2,\epsilon_3,\epsilon_5)")

p11 = plot(size=(800, 600), minorgrid=true, margins=3mm)
plot!(clist, ac60list, lc=2, lw=4, ls=:dot, label="Xe+Xe " * L"\beta=0.207,\gamma=60^\circ")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
plot!(clist, ac27list, lc=3, lw=4, ls=:dash, label="Xe+Xe " * L"\beta=0.207,\gamma=27^\circ")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
plot!(clist, ac0list, lc=4, lw=4, ls=:dashdot, label="Xe+Xe " * L"\beta=0.207,\gamma=0^\circ")
xlims!(0, 25)
ylims!(0.8, 1.7)
xlabel!("centrality(%)")
ylabel!("Xe+Xe/Pb+Pb")
annotate!(11, 0.86, text(L"\rho(\mathcal{E}_2,\mathcal{E}_2,\mathcal{E}_4)" * " ratio vs " * L"N_{\mathrm{ch}}" * " based centrality ", 20, "Computer Modern"))
savefig(p11, "figures/nac_224.svg")

###################################################
# test
using Plots, LaTeXStrings, Measures
plot(size=(800, 600), minorgrid=true, margins=3mm)
plot!(clist, Xe60_ρlist, lc=2, lw=4, ls=:dot, label="Xe+Xe " * L"\beta=0.207,\gamma=60^\circ")
plot!(clist, Xe27_ρlist, lc=3, lw=4, ls=:dash, label="Xe+Xe " * L"\beta=0.207,\gamma=27^\circ")
plot!(clist, Xe0_ρlist, lc=4, lw=4, ls=:dashdot, label="Xe+Xe " * L"\beta=0.207,\gamma=0^\circ")
xlims!(0, 25)
ylims!(-0.01, 0.05)
xlabel!("centrality(%)")
ylabel!(L"\rho(\varepsilon_4^2,\Psi_4^2)")

p7 = plot(size=(800, 600), minorgrid=true, margins=3mm)
plot!(clist, Pb_ρlist, lc=1, lw=4, label="Pb+Pb " * L"\beta=0.062,\gamma=27^\circ")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
plot!(clist, Xe60_ρlist, lc=2, lw=4, ls=:dot, label="Xe+Xe " * L"\beta=0.207,\gamma=60^\circ")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
plot!(clist, Xe27_ρlist, lc=3, lw=4, ls=:dash, label="Xe+Xe " * L"\beta=0.207,\gamma=27^\circ")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
plot!(clist, Xe0_ρlist, lc=4, lw=4, ls=:dashdot, label="Xe+Xe " * L"\beta=0.207,\gamma=0^\circ")
xlims!(0, 25)
ylims!(0, 0.16)
xlabel!("centrality(%)")
ylabel!(L"\rho(\cos(2\Psi_2),\langle{p_T}\rangle)")
annotate!(12.5, 0.02, text(L"\rho(\cos(2\Psi_2),\langle{p_T}\rangle)" * " vs " * L"N_{\mathrm{ch}}" * " based centrality ", 20, "Computer Modern"))
savefig(p7, "figures/gamma/rho(cos(2Psi2),pT).svg")

p8 = plot(size=(800, 600), minorgrid=true, margins=3mm)
plot!(clist, ρ60list, lc=2, lw=4, ls=:dot, label="Xe+Xe " * L"\beta=0.207,\gamma=60^\circ")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
plot!(clist, ρ27list, lc=3, lw=4, ls=:dash, label="Xe+Xe " * L"\beta=0.207,\gamma=27^\circ")
scatter!([30], [1], ms=0, mc=:white, msc=:white, label=" ")
plot!(clist, ρ0list, lc=4, lw=4, ls=:dashdot, label="Xe+Xe " * L"\beta=0.207,\gamma=0^\circ")
xlims!(0, 25)
ylims!(0.5, 1.2)
xlabel!("centrality(%)")
ylabel!("Xe+Xe/Pb+Pb")
# annotate!(15, 0.56, text(L"\rho(\cos(2\Psi_2),\langle{p_T}\rangle)" * " ratio vs " * L"N_{\mathrm{ch}}" * " based centrality ", 20, "Computer Modern"))
annotate!(17, 0.63, text(L"\rho(\cos(2\Psi_2),\langle{p_T}\rangle)" * " ratio vs ", 20, "Computer Modern"))
annotate!(17, 0.56, text(L"N_{\mathrm{ch}}" * " based centrality ", 20, "Computer Modern"))
savefig(p8, "figures/gamma/Xe+Xe_Pb+Pb_cos2Psi2.svg")

###############################################################################################
using JLD2, Statistics, LinearAlgebra
using Dierckx, Polynomials
using Plots, LaTeXStrings, Measures

Xe27_centrality = load("Xe_centrality_0_30.jld2")["Nch"]
Xe27_allevents = load("triaxial/Xe_epsilon2345_2023-05-04.jld2")["data"]
Xe27_Nch = [event[1] for event in Xe27_allevents]

index1 = findfirst(x -> x < Xe27_centrality[4], Xe27_Nch)
index2 = findfirst(x -> x < Xe27_centrality[5], Xe27_Nch)
chooseevents = rand(Xe27_allevents[index1:(index2-1)],1000)

pT = [event[2] for event in chooseevents]
Ψ22 = [cos(angle(-event[3])) for event in chooseevents]
# Ψ32 = [(angle(-event[4])/3)^2 for event in chooseevents]
# Ψ52 = [(angle(-event[6])/5)^2 for event in chooseevents]

p10 = plot(size=(800, 600), minorgrid=true, margins=3mm)
scatter!(Ψ22, pT)
xlabel!(L"\Psi_2^2")
ylabel!(L"⟨ p_T ⟩")
ylabel!(L"\Psi_5^2")