"""
A `Lattice3Vector` is a (position or momentum) 3-vector in (`Direct` or `Reciprocal`) lattice units.
"""
# The signature for Lattice3Vector could include the type of its internal Lattice, allowing for
# multiple dispatch calls based on Direct vs Reciprocal, but since this is typically not enough
# information (e.g., dot needs to know if the two Lattices are *the same* not just of the same type)
# dismiss this possibility for the moment.
#type Lattice3Vector{T<:Real} <: LatticeVector
immutable Lattice3Vector{T<:Number} <: LatticeVector
    # Reciprocal or real-space vector in reciprocal lattice (r.l.u.) or direct lattice units (fractional coordinates), respectively.
    # Two vectors with the same (hkl) indicies in the same lattice are identical, ergo type Lattice3Vector is immutable.
    h::T
    k::T
    l::T
    lat::Lattice
end
Lattice3Vector{Q<:Number,R<:Number,S<:Number}(h::Q,k::R,l::S,lat::Lattice)=(T=promote_type(Q,R,S);Lattice3Vector{T}(T(h),T(k),T(l),lat))
Lattice3Vector(lat::Lattice)=Lattice3Vector(0,0,0,lat)
## Add initializer that takes a 3-element vector and a Lattice
Lattice3Vector{T<:Number}(v::AbstractArray{T,1},l::Lattice)=(@assert length(v)==3;Lattice3Vector(v[1],v[2],v[3],l))
function Lattice3Vector{T<:Number}(m::AbstractArray{T,2},l::Lattice=one(Reciprocal{T}))
    @assert size(m,1)==3
    out=Array{Lattice3Vector{T}}(size(m,2))
    for i=1:size(m,2);out[i]=Lattice3Vector{T}(m[1,i],m[2,i],m[3,i],l);end
    out
end
function Lattice3Vector{T<:Number}(m::AbstractArray{T,3},l::Lattice=one(Reciprocal{T}))
    (a,b,c)=size(m)
    @assert size(m,1)==3
    out=Array{Lattice3Vector{T}}(b,c)
    @inbounds for j=1:c
        @inbounds for i=1:b;
            out[i,j]=Lattice3Vector{T}(m[1,i,j],m[2,i,j],m[3,i,j],l);
        end
    end
    return out
end

function showLattice3Vector(io::IO,m::Lattice3Vector,compact::Bool=false)
    Base.print(io,"(")
    compact ? Base.showcompact(io,m.h) : Base.show(io,m.h)
    Base.print(io," ")
    compact ? Base.showcompact(io,m.k) : Base.show(io,m.k)
    Base.print(io," ")
    compact ? Base.showcompact(io,m.l) : Base.show(io,m.l)
    Base.print(io,")")
    isa(m.lat,Reciprocal)? Base.print(io,"ʳ") : Base.print(io,"ᵈ") # reciprocal (ᵣ) or direct (ᵈ) lattice units)
end
Base.show(io::IO,m::Lattice3Vector) = showLattice3Vector(io,m,false)
Base.showcompact(io::IO,m::Lattice3Vector)=showLattice3Vector(io,m,true)
Base.zero{T<:Number}(::Type{Lattice3Vector{T}})=Lattice3Vector{T}(zero(T),zero(T),zero(T),zero(Reciprocal{T}))
(==)(a::Lattice3Vector,b::Lattice3Vector)=(norm(a-b)==0)
#Base.hash{T<:Real}(a::Type{Lattice3Vector{T}},h::UInt=0)=for f in fieldnames(a); h=Base.hash(a.(f),h); end


Base.isapprox(a::Lattice3Vector,b::Lattice3Vector)=isapprox(norm(a-b),0)



isReciprocal(a::Lattice3Vector)=isa(a.lat,Reciprocal)
isDirect(a::Lattice3Vector)=isa(a.lat,Direct)

# Return the Miller indexes of a vector defined by its fractional cell coordinates, or vice versa.
#star(x::Lattice3Vector)=Lattice3Vector(2*pi*contravariantMetricTensor(x.lat)*[x.h,x.k,x.l],star(x.lat))
star(x::Lattice3Vector)=Lattice3Vector(covariantMetricTensor(x.lat)*[x.h,x.k,x.l]/2/pi,star(x.lat))

"""
The scalar product between two 3-vectors in a (`Direct` or `Reciprocal`) `Lattice` basis is defined
if both vectors are described in the same `Lattice` or in one `Lattice` and its `star`.
The function `dot` performs the scalar product or throws an error for any case that the scalar
product is not defined.
"""
# function Base.dot(x::Lattice3Vector,y::Lattice3Vector)
#     # The scalar product between two (reciprocal)lattice vectors
#     if sameLattice(x.lat,y.lat)
#         L=x.lat
#         d=x.h*y.h*L[:a]^2 + x.k*y.k*L[:b]^2 + x.l*y.l*L[:c]^2
#         d+=(x.h*y.k+x.k*y.h)*L[:a]*L[:b]*cos(L[:γ])
#         d+=(x.h*y.l+x.l*y.h)*L[:c]*L[:a]*cos(L[:β])
#         d+=(x.l*y.k+x.k*y.l)*L[:b]*L[:c]*cos(L[:α])
#     elseif invLattices(x.lat,y.lat)
#         typeof(geth(x))(2*pi)*Base.dot([x.h,x.k,x.l],[y.h,y.k,y.l])
#     else
#         println(x.lat," and ",y.lat," are not inverse or same!")
#         sameLattice(x.lat,y.lat)
#         Lattices.sameLattice(x.lat,y.lat)
#         isapprox(x.lat,y.lat)
#         error("Both vectors must be described in the same or inverse lattices.")
#     end
# end
function Base.dot(x::Lattice3Vector,y::Lattice3Vector)
    # The scalar product between two (reciprocal)lattice vectors
    if sameLattice(x,y)
        L=x.lat
        d=x.h*y.h*L[:a]^2 + x.k*y.k*L[:b]^2 + x.l*y.l*L[:c]^2
        d+=(x.h*y.k+x.k*y.h)*L[:a]*L[:b]*cos(L[:γ])
        d+=(x.h*y.l+x.l*y.h)*L[:c]*L[:a]*cos(L[:β])
        d+=(x.l*y.k+x.k*y.l)*L[:b]*L[:c]*cos(L[:α])
    elseif invLattices(x,y)
        typeof(geth(x))(2*pi)*Base.dot([x.h,x.k,x.l],[y.h,y.k,y.l])
    else
        println(x.lat," and ",y.lat," are not inverse or same!")
        error("Both vectors must be described in the same or inverse lattices.")
    end
end
Base.norm(x::Lattice3Vector)=Base.sqrt(Base.dot(x,x)) # The norm of a lattice vector is just the square root of its scalar product with itself
Base.norm{L<:Lattice3Vector}(v::Array{L})=map(Base.norm,v)
"""
The vector product between two vectors in a lattice gives a third vector most easily expressed in
the reciprocal lattice. If v⃑₁= h₁a⃑ + k₁b⃑ + l₁c⃑, and v⃑₂ = h₂a⃑ + k₂b⃑ + l₂c⃑, then v⃑₁×v⃑₂ gives

v⃑₁×v⃑₂=(k₁l₂-k₂l₁) b⃑×c⃑ + (h₂l₁-h₁l₂) c⃑×a⃑ + (h₁k₂-h₂k₁) a⃑×b⃑

v⃑₁×v⃑₂=V/2π * [(k₁l₂-k₂l₁) a⃑⋆ + (h₂l₁-h₁l₂) b⃑⋆ + (h₁k₂-h₂k₁) c⃑⋆]

The function `cross` ensures that two `Lattice3Vector` objects are described in the same `Lattice`,
creates a new `Lattice3Vector` V/2π*[k₁l₂-k₂l₁,h₂l₁-h₁l₂,h₁k₂-h₂k₁] in `star(`[original `Lattice`]`)`
and finaly utilizes `star` to change the basis of the new vector back to the original `Lattice`.
"""
function Base.cross(x::Lattice3Vector,y::Lattice3Vector)
    (x.lat == y.lat) ? (L=x.lat) : (error("Only defined for two vectors within the same lattice"))
    # cross([h1,k1,l1],[h2,k2,l2])*V/2π  in the reciprocal lattice coordinate system
    z=star(Lattice3Vector(Base.cross([x.h,x.k,x.l],[y.h,y.k,y.l])* L[:V]/2/pi , star(L))) # the cross product in the original lattice basis
    Lattice3Vector(z.h,z.k,z.l,L) # in case star(star(L)) introduced rounding errors
end
"""
    orthonormalize(v1::Lattice3Vector, v2::Lattice3Vector)

Given two non-colinear vectors in a `Direct` or `Reciprocal` `Lattice`, `orthonormalize` creates
a set of three vectors that describe a right-handed orthonormal basis expressed as `Latice3Vector`s
in the same `Lattice`. The returned orthonormal basis [`x`,`y`,`z`] has `x` along `v1` (`v1⋅x=|v1|`),
`y` along the component of `v2` perpendicular to `x`, and `z=x×y`.

An optional third `Lattice3Vector`, `v3`, can be passed to `orthonormalize` to force a left-handed
orthonormal basis if `v3⋅(x×y) < 0`.
"""
function orthonormalize(o1::Lattice3Vector,o2::Lattice3Vector,o3::Lattice3Vector=cross(o1,o2))
    # From two (reciprocal) lattice vectors o1 & o2, determine an orthonormal set of three vectors x, y & z
    x=o1/Base.norm(o1)    # ensure x is a unit vector along o1
    y=o2-x*Base.dot(o2,x) # first set y to the component of o2 perpendicular to x
    (Base.norm(y) > 0) ? y=y/Base.norm(y) : error("x & y are colinear.")
    z=Base.cross(x,y)
    (Base.dot(z,o3) < 0) && (z=-z) # allow for left-handed coordinate systems (maybe this is useful?)
    z=z/Base.norm(z)
    (x,y,z)
end

import Base: angle
"""
    angle(x::Lattice3Vector,y::Lattice3Vector)
Returns the angle between two `Lattice3Vectors` which are descirbed in the same `Lattice` or have
lattices which are the `star` of one another.
"""
angle(x::Lattice3Vector,y::Lattice3Vector)= Base.acos(Base.dot(x,y)/(Base.norm(x)*Base.norm(y)))

# rounding a Lattice3Vector returns the near-by vector with its incicies rounded to the nearest integer
Base.round{T<:Integer}(::Type{T},a::Lattice3Vector)=Lattice3Vector(Base.round(T,a.h),Base.round(T,a.k),Base,round(T,a.l),a.lat)
Base.round{T<:AbstractFloat}(::Type{T},a::Lattice3Vector)=Lattice3Vector(Base.round(a.h),Base.round(a.k),Base.round(a.l),a.lat)
Base.round(a::Lattice3Vector)=Lattice3Vector(Base.round(a.h),Base.round(a.k),Base.round(a.l),a.lat)

Base.ceil{T<:Real}(::Type{T},a::Lattice3Vector)=Lattice3Vector(Base.ceil(T,a.h),Base.ceil(T,a.k),Base.ceil(T,a.l),a.lat)
Base.floor{T<:Real}(::Type{T},a::Lattice3Vector)=Lattice3Vector(Base.floor(T,a.h),Base.floor(T,a.k),Base.floor(T,a.l),a.lat)
"""
    roundfinite(T::Type,q::Lattice3Vector)
Returns `Lattice3Vector{T}` that has indicies of `q` rounded to the nearest integer.
If the rounded result is `[0,0,0]` returns the closest of {[100],[010],[001]} and
their Friedel pairs to `q`
"""
function roundfinite{T<:Real}(::Type{T},a::Lattice3Vector)
    b=Base.round(T,a)
    Base.norm(b)>0&&(return b)
    c=Lattice3Vector(T[1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1]',a.lat)
    d=map(Base.norm,c-a) # find the distance to each element in c
    c[findmin(d)[2]] # return the closest element of c to a
end
roundfinite{T<:Real}(a::Lattice3Vector{T})=roundfinite(T,a)
roundfinite{T<:Real}(a::Array{Lattice3Vector{T}})=roundfinite(T,a)
roundfinite{S<:Real,T<:Real}(::Type{S},a::Array{Lattice3Vector{T}})=map(x->roundfinite(T,x),a)

#Base.getindex(a::Lattice3Vector,s::Symbol)= any(s.==fieldnames(Lattice3Vector))?a.(s):a.lat[s]
Base.getindex(a::Lattice3Vector,s::Symbol)= any(s.==fieldnames(a))?a.(s):a.lat[s]
geth(a::Lattice3Vector)=a.h
getk(a::Lattice3Vector)=a.k
getl(a::Lattice3Vector)=a.l
getE(a::Lattice3Vector)=zero(typeof(a.h))
getlat(a::Lattice3Vector)=a.lat

get3vector(a::Lattice3Vector)=a
get4vector(a::Lattice3Vector)=Lattice4Vector(a,0)
