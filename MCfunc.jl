using MKL, LinearAlgebra, Distributions, SpecialFunctions
using HCubature    # 数值积分包

@inline sph2car(r, θ, ϕ) = (r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ))
@inline car2sph(x, y, z) = (hypot(x, y, z), atan(hypot(x, y), z), atan(y, x))
@inline distance(x1, x2) = hypot((x1 - x2)...)

# 将球坐标分布函数转化为直角坐标分布函数
function f_car(f_sph, x)
  r, θ, ϕ = car2sph(x[1], x[2], x[3])
  return f_sph(r, θ, ϕ)
end

# 产生球rmax内部的三维分布点
function ball_even(rmax)
  r = cbrt(rand()) * rmax      # 以正比于r^2概率产生r
  θ = acos(2 * rand() - 1)     # 以正比于sin(θ)概率产生θ
  ϕ = 2π * rand()              # 均匀产生φ
  return r, θ, ϕ
end
# 产生一个服从f_sph分布的样本点
function sample_car(f_sph, rmax)
	r, θ, ϕ = ball_even(rmax)
  ρ = f_sph(r, θ, ϕ)     # 注意分布函数f的最大值不大于1
  while ρ < rand()
    r, θ, ϕ = ball_even(rmax)
    ρ = f_sph(r, θ, ϕ)
  end
  return sph2car(r, θ, ϕ)    # 返回直角坐标
end

# 产生球均匀分布的径向坐标
ball_even_radial(rmax) = cbrt(rand()) * rmax    # 以正比于r^2概率产生r

# 产生一个服从球对称分布f_r的样本点
function sample_sym_car(f_r, rmax)
  r, θ, ϕ = ball_even(rmax)
  ρ = f_r(r)     # 注意分布函数f的最大值不大于1
  while ρ < rand()
    r = ball_even_radial(rmax)
    ρ = f_r(r)
  end
  return sph2car(r, θ, ϕ)    # 返回直角坐标
end

#################################################################
# 以原子核为单位
# 增加核子距离判断
function add_new_nucleon!(f_sph, nucleons, dmin, rmax)
  p = (0.0, 0.0, 0.0)
  acc = false
  while ~acc
    p = sample_car(f_sph, rmax)
    acc = all(q -> distance(p, q) > dmin, nucleons)
    # acc = all([distance(p, collect(nucleons[:, i])) > dmin for i in 1:size(nucleons)[2]])    # 如果用matrix储存就这样写
  end
  push!(nucleons, p)
end

# 调整样本的分布中心为原点
function centralize!(samples)
  c = mean(samples)
  @inbounds for i ∈ eachindex(samples)
    samples[i] = samples[i] - c
  end
end

# 生成N个核子的原子核
function get_nucleus(f_sph, N, dmin, rmax)
  nucleons = NTuple{3,Float64}[]
  sizehint!(nucleons, N)
  for _ = 1:N
    add_new_nucleon!(f_sph, nucleons, dmin, rmax)
  end
  centralize!(nucleons)
  return nucleons
end

# 随机旋转，并将核子坐标投影到截面上
function rotation(samples)
  α = rand() * 2π
  β = acos(2 * rand() - 1)
  γ = rand() * 2π
  R1 = [cos(α)*cos(γ)-cos(β)*sin(α)*sin(γ), sin(α)*cos(γ)+cos(β)*cos(α)*sin(γ), sin(β)*sin(γ)]
  R2 = [-cos(α)*sin(γ)-cos(β)*sin(α)*cos(γ), -sin(α)*sin(γ)+cos(β)*cos(α)*cos(γ), sin(β)*cos(γ)]
  # R3 = [sin(β)*sin(α), -sin(β)*cos(α), cos(β)]
  samples_rot = NTuple{2,Float64}[]
  for sample in samples
    x_rot = dot(R1, sample)
    y_rot = dot(R2, sample)
    # z_rot = dot(R3, sample)
    push!(samples_rot, (x_rot,y_rot))
  end
  return samples_rot
end

# 生成两原子核碰撞的相对坐标
function collision_position(A, B, b)    # 把两原子核中心连线都认为转到x轴
  A_rot = rotation(A)
  B_rot = rotation(B)
  A = [n + (b / 2, 0.0) for n in A_rot]
  B = [n - (b / 2, 0.0) for n in B_rot]
  return A, B
end

# # 增加核子距离判断
# function add_new_nucleon!(f_r, nucleons, dmin, rmax)
#   p = (0.0, 0.0, 0.0)
#   acc = false
#   while ~acc
#     p = sample_sym_car(f_r, rmax)
#     acc = all(q -> distance(p, q) > dmin, nucleons)
#   end
#   push!(nucleons, p)
# end

# # 调整样本的分布中心为原点
# function centralize!(samples)
#   c = mean(samples)
#   @inbounds for i ∈ eachindex(samples)
#     samples[i] = samples[i] - c
#   end
# end

# # 生成N个核子的原子核（由相对论效应，只关心截面投影坐标即可）
# function get_nucleus(f_r, N, dmin, rmax)
#   nucleons = NTuple{3,Float64}[]
#   sizehint!(nucleons, N)
#   for _ = 1:N
#     add_new_nucleon!(f_r, nucleons, dmin, rmax)
#   end
#   centralize!(nucleons)
#   return NTuple{2,Float64}[(n[1], n[2]) for n ∈ nucleons]
# end

# # 生成两原子核碰撞的相对坐标
# function collision_position!(A, B, b)    # 把两原子核中心连线都认为转到x轴
#   bias = (b / 2, 0.0)
#   @inbounds for i ∈ eachindex(A)
#     A[i] += bias
#   end
#   @inbounds for i ∈ eachindex(B)
#     B[i] -= bias
#   end
# end

###############################################################
## 以夸克为单位
# 生成 N*nq个夸克在截面上的投影坐标
function get_quarks(nucleus, nq, quark_pdf, rmax)
  quarks = [n for n in nucleus for _ ∈ 1:nq]
  @inbounds for n ∈ 1:length(nucleus)
    quark = [sample_sym_car(quark_pdf, rmax) for _ ∈ 1:nq]
    centralize!(quark)
    qs = NTuple{2,Float64}[(q[1], q[2]) for q ∈ quark]
    quarks[(3*n-2):3*n] += qs
  end
  return quarks
end

# 判断 wounded quark数目和 binary collision

@inline function bcollision(n1, n2, σqq)
  b = distance(n1, n2)
  p = exp(-π * b^2 / σqq)
  return p > rand()
end

# function quark_wounding(A, B, σqq)
#   Ap = falses(size(A))
#   Bp = falses(size(B))
#   # binary_collisions = Tuple{Int,Int}[]
#   Nbin = 0
#   for i in eachindex(A), j in eachindex(B)
#     if bcollision(A[i], B[j], σqq)
#       # push!(binary_collisions, (i, j))
#       Ap[i] = Bp[j] = true
#       Nbin = Nbin + 1
#     end
#   end
#   participants = sum(Ap) + sum(Bp)
#   return participants, Nbin #, binary_collisions
# end

@inline function get_positions(positions, A, Ap)
  for i in eachindex(A)
    Ap[i] ? push!(positions, A[i]) : positions = positions
  end
  return positions
end

function quark_wounding(A, B, σqq)
  Ap = falses(size(A))
  Bp = falses(size(B))
  # binary_collisions = Tuple{Int,Int}[]
  for i in eachindex(A), j in eachindex(B)
    if bcollision(A[i], B[j], σqq)
      # push!(binary_collisions, (i, j))
      Ap[i] = Bp[j] = true
    end
  end
  participants = sum(Ap) + sum(Bp)
  positions = NTuple{2,Float64}[]
  positions = get_positions(positions, A, Ap)
  positions = get_positions(positions, B, Bp)
  return participants, positions#, binary_collisions
end

####################################################################################
# 计算charged partical数目
function charged_partical(i, para=(n̄, κ))
  n̄_, κ_ = para
  ni = rand(NegativeBinomial(κ_, κ_ / (n̄_ + κ_)), i)
  return sum(ni), ni
end

# 生成熵密度分布函数
# function entropy_density(x, wounded_quarks, nk, para=(σ,τ))
#   σ, τ = para
#   g = [exp(-(distance(x, xi)^2/2/σ^2)) for xi in wounded_quarks]
#   ω = 6 .* nk    # 为每个participants加权
#   s = sum(g.*ω) / (2*π*σ^2) / τ
#   return s    # 单位 fm^-3
# end

function entropy_density(x, quarks_nk)::Float64
  s = sum(n * exp(-(abs2(x[1] - xi[1]) + abs2(x[2] - xi[2])) / (2σ^2)) for (xi, n) ∈ quarks_nk)
  s = s * 6 / (2 * π * σ^2) / τ
  return s
end

# const xTmax::Float64 = rmax * 2 / 3
const zeta3::Float64 = zeta(3)

# 计算被积函数在局域的值
function local_obs(x, quarks_nk)
  # [s, pt, r2real, r2imag, r2abs]
  s = entropy_density(x, quarks_nk)
  pt = s^(4/3)
  rc = complex(x[1],x[2])
  r2 = s * rc^2
  r3 = s * rc^3
  r4 = s * rc^4
  r5 = s * rc^5
  obs = [s, pt,
        real(r2), imag(r2), abs(r2), real(r3), imag(r3), abs(r3), 
        real(r4), imag(r4), abs(r4), real(r5), imag(r5), abs(r5)]
  return obs
end

# 计算全空间积分
function integrate_obs(quarks_nk)
  xrange = extrema(r[1] for (r, _) in quarks_nk) + (-3σ, 3σ)
  yrange = extrema(r[2] for (r, _) in quarks_nk) + (-3σ, 3σ)
  obs_int,_ = hcubature(
    x -> local_obs(x, quarks_nk),
    [xrange[1], yrange[1]], [xrange[2], yrange[2]]; norm = o->abs(o[1]), rtol=1e-3, initdiv=8, maxevals=100000, 
  )
  return obs_int
end

# 计算各物理量
function global_obs(obs_int)
  N, PT, r2_re, r2_im, r2_abs, r3_re, r3_im, r3_abs, r4_re, r4_im, r4_abs, r5_re, r5_im, r5_abs = obs_int
  P̄T = 19 / 1376 / zeta3 * (π^10 / 57)^(1 / 3) * PT / N
  ε2 = complex(r2_re, r2_im) / r2_abs
  ε3 = complex(r3_re, r3_im) / r3_abs
  ε4 = complex(r4_re, r4_im) / r4_abs
  ε5 = complex(r5_re, r5_im) / r5_abs
  return P̄T, ε2, ε3, ε4, ε5
end

# # 计算n阶离心率
# function eccentricity(quarks_nk, n)
#   complexr, _ = hcubature(
#     x -> entropy_density(x, quarks_nk) * (complex(x[1], x[2]))^n,
#     [-xTmax, -xTmax], [xTmax, xTmax]; rtol=1e-3, initdiv=32
#   )
#   absr, _ = hcubature(
#     x -> entropy_density(x, quarks_nk) * (hypot(x...))^n,
#     [-xTmax, -xTmax], [xTmax, xTmax]; rtol=1e-3, initdiv=32
#   )
#   ε = complexr / absr
#   Φ = angle(ε) / n
#   return abs(ε), Φ
# end
