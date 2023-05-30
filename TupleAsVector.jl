using LinearAlgebra
import Base.+, Base.-, Base.*, Base./,Base.hypot,Base.abs2
import LinearAlgebra.dot

@inline function *(R::Matrix, v::NTuple{2})
    @inbounds begin
        return promote(
            R[1] * v[1] + R[3] * v[2],
            R[2] * v[1] + R[4] * v[2],
        )
    end
end

@inline function *(R::Matrix, v::NTuple{3})
    @inbounds begin
        return promote(
            R[1] * v[1] + R[4] * v[2] + R[7] * v[3],
            R[2] * v[1] + R[5] * v[2] + R[8] * v[3],
            R[3] * v[1] + R[6] * v[2] + R[9] * v[3],
        )
    end
end

@inline function *(R::Matrix, v::NTuple{4})
    @inbounds begin
        return promote(
            R[01] * v[1] + R[05] * v[2] + R[09] * v[3] + R[13] * v[4],
            R[02] * v[1] + R[06] * v[2] + R[10] * v[3] + R[14] * v[4],
            R[03] * v[1] + R[07] * v[2] + R[11] * v[3] + R[15] * v[4],
            R[04] * v[1] + R[08] * v[2] + R[12] * v[3] + R[16] * v[4],
        )
    end
end

dot(u::NTuple{2,Real}, v::NTuple{2}) = u[1] * v[1] + u[2] * v[2]
dot(u::NTuple{3,Real}, v::NTuple{3}) = u[1] * v[1] + u[2] * v[2] + u[3] * v[3]
dot(u::NTuple{4,Real}, v::NTuple{4}) = u[1] * v[1] + u[2] * v[2] + u[3] * v[3] + u[4] * v[4]
dot(u::NTuple{2,Complex}, v::NTuple{2}) = u[1]' * v[1] + u[2]' * v[2]
dot(u::NTuple{3,Complex}, v::NTuple{3}) = u[1]' * v[1] + u[2]' * v[2] + u[3]' * v[3]
dot(u::NTuple{4,Complex}, v::NTuple{4}) = u[1]' * v[1] + u[2]' * v[2] + u[3]' * v[3] + u[4]' * v[4]

import Base.hypot
@inline hypot(u::NTuple{2,Real}, v::NTuple{2,Real}) = hypot(u[1] - v[1], u[2] - v[2])
@inline hypot(u::NTuple{3,Real}, v::NTuple{3,Real}) = hypot(u[1] - v[1], u[2] - v[2], u[3] - v[3])
@inline hypot(u::NTuple{4,Real}, v::NTuple{4,Real}) = hypot(u[1] - v[1], u[2] - v[2], u[3] - v[3], u[4] - v[4])
@inline hypot(u::NTuple,v::NTuple) = hypot(u-v)

import Base.abs2
@inline abs2(u::NTuple{2,Real}, v::NTuple{2,Real}) = abs2(u[1] - v[1]) + abs2(u[2] - v[2])
@inline abs2(u::NTuple{3,Real}, v::NTuple{3,Real}) = abs2(u[1] - v[1]) + abs2(u[2] - v[2]) + abs2(u[3] - v[3])
@inline abs2(u::NTuple{4,Real}, v::NTuple{4,Real}) = abs2(u[1] - v[1]) + abs2(u[2] - v[2]) + abs2(u[3] - v[3]) + abs2(u[4] - v[4])
@inline abs2(u::NTuple,v::NTuple) = sum(abs2,u-v)

+(u::NTuple, v::NTuple) = u.+v
-(u::NTuple, v::NTuple) = u.-v
*(α::Number, v::NTuple) = α.*v
/(v::NTuple, β::Number) = v./β

# # [test]
# N = 3
# testtype = Float64
# v1 = rand(testtype, N)
# v2 = rand(testtype, N)
# v1t = (v1...,)
# v2t = (v2...,)
# A = rand(testtype, N,N)

# v1 ⋅ v2 == v1t ⋅ v2t
# A * v1 .- A * v1t

# using BenchmarkTools

# @btime *($A,$v1)
# @btime *($A,$v1t)


