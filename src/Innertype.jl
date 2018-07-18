types=[:Lattice3Vector,:Lattice4Vector,:Lattice3Tensor,:Lattice4Tensor,:Sample,:Crystal]
for x in types
    @eval innertype{T<:Real}(::$x{T})=T
    @eval innertype{T<:Real}(a::Array{$x{T}})=T
end
