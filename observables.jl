using LinearAlgebra, Statistics
using Distributed, SharedArrays
using Plots, LaTeXStrings
using ProgressMeter
events = 6_000_000
# pbar = Progress(events)
# channel = RemoteChannel(() -> Channel{Bool}(), 1)

## centrality对应的Nch
using JLD2, Dates
centrality = load("Xe_centrality_0_30.jld2")["Nch"]

## 计算某中心度对应的各统计量
# versioninfo()
if length(workers())==1
  addprocs(256)
end
workers()
wp = CachingPool(workers())

@everywhere include("MCmodel.jl")
@everywhere using ProgressMeter

# filename = "triaxial/Xe_0-25_$(today()).jld2"
# for i in 1:25
#   println("task $(i) |  centrality: $(i-1)-$(i)%")
# 	name = "d"*string(i)
# 	st = (i == 1 ? 10000 : centrality[i-1])
#   ed = centrality[i]
# 	data = progress_pmap(_->get_observables(st, ed), wp, 1:events)
# 	jldopen(filename, "a+") do f
# 		f[name] = data
# 	end
# end
get_observables(10000, centrality[25])
# data = progress_pmap(_ -> get_charged_particals(10000, centrality[25]), wp, 1:events)

filename = "triaxial/Xe27_$(today()).jld2"
data = progress_pmap(_->get_observables(10000, centrality[25]), wp, 1:events)
sort!(data,by=x->-x[1])
jldsave(filename, data=data)

rmprocs(workers())    # 释放worker
workers()

# allevents = load("triaxial/Pb2TeV_epsilon2345_2023-05-07.jld2")["data"]
allevents = load("triaxial/Pb_2023-05-26.jld2")["data"]
Nch = [event[1] for event in allevents]

clist = 0.5:1:24.5
index = [findfirst(x -> x < centrality[i], Nch) for i in 1:length(clist)]

data = Vector{Tuple{Int64, Float64, Vararg{ComplexF64, 4}}}[]
# data = Vector{Tuple{Int64, Int64}}[]
for i in 1:length(clist)
  st = (i == 1 ? 1 : index[i-1])
  ed = (i == 25 ? events : index[i] - 1)
  push!(data, allevents[st:ed])
end

# for i in 1:25
#   x = load("triaxial/Xe_0-25_2023-05-02.jld2")["d"*string(i)]
#   push!(data, x)
# end
data

## 计算关联系数
function rho(a, b)
  δa = a .- mean(a)
  δb = b .- mean(b)
  ρ = mean(δa .* δb) / sqrt(mean(δa .^ 2) * mean(δb .^ 2))
  return ρ
end
ρlist = zeros(length(data))
for j in 1:length(clist)
  event = data[j]
  sort!(event, by=n->-n[1])
  p̄T = [y[2] for y in event]
  # Nbin = [y[end] for y in event]
  ε2 = [abs(y[3]) for y in event]
  # ε3 = [abs(y[4]) for y in event]
  # ε4 = [abs(y[5]) for y in event]
  # Φ2 = [angle(y[3])/2 for y in event]
  # Φ4 = [angle(y[5]) / 4 for y in event]
  # ρlist[j] = rho(ε4 .^ 2, Φ4 .^ 2)
  ρlist[j] = rho(ε2.^2, p̄T)
end

ρlist

plot(clist, ρlist)
scatter!(clist, ρlist)
xlims!(0, 25)
ylims!(-0.35, 0.15)
xlabel!("centrality(%)")
ylabel!(L"\rho(N_{\mathrm{q-coll}},\langle{p_T}\rangle)")
# ylabel!(L"\rho(\Phi_2^2,\langle{p_T}\rangle)")


#######################################################################
## 算v2和v3的关系
using LinearAlgebra, Statistics
using Distributed, SharedArrays
using Plots, LaTeXStrings
events = 6_000_000

using JLD2
allevents = load("triaxial/Pb2TeV_epsilon2345_2023-05-07.jld2")["data"]
centrality = load("Pb2TeV_centrality_0_30.jld2")["Nch"]
Nch = [event[1] for event in allevents]

# 只要15%-20%的中心度事件
clist = [0.15, 0.20]
st = findfirst(x -> x < centrality[15], Nch)
ed = findfirst(x -> x < centrality[20], Nch)
data = allevents[st:(ed-1)]
events = length(data)

sort!(data, by = x->abs(x[3]))
ε2 = [abs(event[3]) for event in data]
ε3 = [abs(event[4]) for event in data]

xlist = collect(0.01 : 0.02 : 0.57)
index = [findfirst(x -> x > (0.02*i), ε2) for i in 1:(length(xlist)-1)]

ylist = zeros(length(xlist))
for i in 1:length(xlist)
  st = (i == 1 ? 1 : index[i-1])
  ed = (i == length(xlist) ? events : index[i] - 1)
  ylist[i] = mean(ε3[st:ed])
end
ylist

plot(xlist, ylist, label="Pb+Pb at 2.76TeV, centrality 15-20%")
xlims!(0.0, 0.6)
ylims!(0.14, 0.2)
xlabel!(L"\varepsilon_2")
ylabel!(L"\varepsilon_3")

index_plot = collect(1:2000:605666)
ε2_plot = zeros(length(index_plot))
ε3_plot = zeros(length(index_plot))
for i in 1:length(index_plot)
  ε2_plot[i] = ε2[index_plot[i]]
  ε3_plot[i] = ε3[index_plot[i]]
end
plot(ε2_plot, ε3_plot)

using Polynomials
fit(Polynomial, ε2_plot, ε3_plot, 1)