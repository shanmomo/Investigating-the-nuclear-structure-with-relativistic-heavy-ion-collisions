include("TupleAsVector.jl")

## parameters for Xe-129 at √sNN = 5.44 TeV
# woods-saxon profile
N = 129
p = (5.40, 0.59, 0.207, π * 27 / 180)    # p=(r0,d,β,γ)
rmax = 2 * p[1] + 10 * p[2]    # 生成核子的最大范围
dmin = 0.8    # 最小核子距离
Y20(θ) = (3 * cos(θ)^2 - 1) / 4 * √(5 / π)
Y22(θ, ϕ) = cos(2 * ϕ) * sin(θ)^2 / 4 * √(15 / π)    # (Y_2^2+Y_2^-2)/√2
function ws(r, θ, ϕ, para=p)
  r0, d, β, γ = para
  R = r0 * (1 + β * (cos(γ) * Y20(θ) + sin(γ) * Y22(θ, ϕ)))
  f = (1 + exp((r - R) / d))^(-1)
  return f
end
# p = (5.40, 0.59, 0.207)    # p=(r0,d,β2)
# rmax = 2 * p[1] + 10 * p[2]    # 生成核子的最大范围
# dmin = 0.8    # 最小核子距离
# Y20(θ) = (3 * cos(θ)^2 - 1) / 4 * √(5 / π)
# function ws(r, θ, ϕ, para=p)
#   r0, d, β = para
#   R = r0 * (1 + β * Y20(θ))
#   f = (1 + exp((r - R) / d))^(-1)
#   return f
# end
# quark potential
const nq::Int = 3
const r0::Float64 = 0.30    # √sNN=5.02TeV时参数，单位fm
const σqq::Float64 = 1.36   # fm^2, 1mb=0.1fm^2
quark_potential(r) = exp(-r / r0)
# charged partical parameters
# const n̄::Float64 = 1.7036045743335102    # 0.1
# const κ::Float64 = 0.23588371029233218
const n̄::Float64 = 1.706590829306753    # 0.2
const κ::Float64 = 0.23629719175016584
# const n̄::Float64 = 1.7286355370688862    # 0.3
# const κ::Float64 = 0.2393495359018458
# entropy density
const σ::Float64 = 0.55    # fm
const τ::Float64 = 0.4     # fm

# ## parameters for Pb-208 at √sNN = 5.02 TeV
# # woods-saxon profile
# N = 208
# p = (6.62, 0.546, 0.062, π * 27 / 180)    # p=(r0,d,β,γ)
# rmax = 2 * p[1] + 10 * p[2]    # 生成核子的最大范围
# dmin = 0.8    # 最小核子距离
# Y20(θ) = (3 * cos(θ)^2 - 1) / 4 * √(5 / π)
# Y22(θ, ϕ) = cos(2 * ϕ) * sin(θ)^2 / 4 * √(15 / π)    # (Y_2^2+Y_2^-2)/√2
# function ws(r, θ, ϕ, para=p)
#   r0, d, β, γ = para
#   R = r0 * (1 + β * (cos(γ) * Y20(θ) + sin(γ) * Y22(θ, ϕ)))
#   f = (1 + exp((r - R) / d))^(-1)
#   return f
# end
# # quark potential
# const nq::Int = 3
# const r0::Float64 = 0.30    # √sNN=5.02TeV时参数，单位fm
# const σqq::Float64 = 1.36   # fm^2, 1mb=0.1fm^2
# quark_potential(r) = exp(-r / r0)
# # charged partical parameters
# const n̄::Float64 = 1.7437287329172955
# const κ::Float64 = 0.24143936301931784
# # entropy density
# const σ::Float64 = 0.55    # fm
# const τ::Float64 = 0.4     # fm

## parameters for Pb-208 at √sNN = 2.76 TeV
# woods-saxon profile
# N = 208
# p = (6.5, 0.54)    # p=(r0,d) fm
# rmax = 2 * p[1] + 10 * p[2]    # 生成核子的最大范围 fm
# dmin = 0.9    # 最小核子距离 fm
# function ws(r, para=p)
#   r0, d = para
#   f = (1 + exp((r - r0) / d))^(-1)
#   return f
# end
# # quark potential
# const nq::Int = 3
# const r0::Float64 = 0.29    # √sNN=2.76TeV时参数，单位fm
# const σqq::Float64 = 1.19   # fm^2, 1mb=0.1fm^2
# quark_potential(r) = exp(-r / r0)
# # charged partical parameters
# const n̄::Float64 = 1.4440102839407825
# const κ::Float64 = 0.28364487720265374

# N = 129
# p = (5.40, 0.59)    # p=(r0,d) fm
# rmax = 2 * p[1] + 10 * p[2]    # 生成核子的最大范围 fm
# dmin = 0.8    # 最小核子距离 fm
# function ws(r, para=p)
#   r0, d = para
#   f = (1 + exp((r - r0) / d))^(-1)
#   return f
# end
# # quark potential
# const nq::Int = 3
# const r0::Float64 = 0.30    # √sNN=5.02TeV时参数，单位fm
# const σqq::Float64 = 1.36   # fm^2, 1mb=0.1fm^2
# quark_potential(r) = exp(-r / r0)
# # charged partical parameters
# const n̄::Float64 = 1.706590829306753
# const κ::Float64 = 0.23629719175016584
# # entropy density
# const σ::Float64 = 0.55    # fm
# const τ::Float64 = 0.4     # fm

## wounded quark model
include("MCfunc.jl")

# # 看下和binary collision的关系
# function get_charged_particals(Nch_max::Int, Nch_min::Int)
#   Nch = 0
#   Nbin = 0
#   while (Nch >= Nch_max) || (Nch < Nch_min)
#     b = sqrt(rand()) * 2 * rmax    # 碰撞参数
#     nucleusA = get_nucleus(ws, N, dmin, rmax)    # 返回N个核子的(x,y,z)坐标
#     nucleusB = get_nucleus(ws, N, dmin, rmax)    # 使用Xe+Xe碰撞
#     nucleusA, nucleusB = collision_position(nucleusA, nucleusB, b)    # 旋转任意角度，在x轴上距离b碰撞
#     quarkA = get_quarks(nucleusA, nq, quark_potential, 6 * r0)
#     quarkB = get_quarks(nucleusB, nq, quark_potential, 6 * r0)
#     n_wq, Nbin = quark_wounding(quarkA, quarkB, σqq)
#     Nch, nk = charged_partical(n_wq)
#   end
#   return Nch, Nbin
# end

function get_charged_particals()
  b = sqrt(rand()) * 2 * rmax    # 碰撞参数
  nucleusA = get_nucleus(ws, N, dmin, rmax)    # 返回N个核子的(x,y,z)坐标
  nucleusB = get_nucleus(ws, N, dmin, rmax)    # 使用Xe+Xe碰撞
  nucleusA, nucleusB = collision_position(nucleusA, nucleusB, b)    # 旋转任意角度，在x轴上距离b碰撞
  # collision_position!(nucleusA, nucleusB, b)    # 在x轴上距离b碰撞
  quarkA = get_quarks(nucleusA, nq, quark_potential, 6 * r0)
  quarkB = get_quarks(nucleusB, nq, quark_potential, 6 * r0)
  n_wq, quarks = quark_wounding(quarkA, quarkB, σqq)
  Nch, nk = charged_partical(n_wq)
  quarks_nk = [qn for qn ∈ zip(quarks,nk)]
  return quarks_nk, Nch
end

function get_observables(Nch_max::Int, Nch_min::Int)
  quarks_nk, Nch = get_charged_particals()
  while (Nch >= Nch_max) || (Nch < Nch_min)
    quarks_nk, Nch = get_charged_particals()
  end
  filter!(qn->qn[2]≠0, quarks_nk)
  p̄T, ε2, ε3, ε4, ε5 = global_obs(integrate_obs(quarks_nk))
  observables = (Nch, p̄T, ε2, ε3, ε4, ε5)
  return observables
end
