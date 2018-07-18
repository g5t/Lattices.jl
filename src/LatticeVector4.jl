"""
As implemented, `Lattice4Vector` is only properly defined for a `Reciprocal` `Lattice` and
has fields for momentum-transfer, `h`, `k`, and `l`, and energy-transfer in meV, `E`, plus
`lat` to hold a `Reciprocal` or `Direct` type `Lattice`.

The `Lattice4Vector` of a `Direct` `Lattice` is well defined, with the 3-vector part
representing a vector in real space and the fourth component representing a time.
However, since the `Direct` `Lattice4Vector` is not used (and there exists an ambiguity due to
choice of time-units) no conversion between `Direct` and `Reciprocal` `Lattice4Vector`s has
yet been implemented.

| `Lattice4Vector(`...`)` | Initialize |
|:---|:---|
|Q`::Lattice3Vector`,E`::Number`| standard, with promotion of `Q` and `E` to same `Type<:Number`|
|h`::Number`,k`::Number`,l`::Number`,E`::Number`,lat`::Lattice`| individual components specified|
|[h;k;l;E]`::Vector{Number}`,lat`::Lattice` | individual components specified as a `Julia` built-in `Vector{Number}`|
|lat`::Lattice`| the zero-vector in a given `Lattice`|
|Q`::Lattice3Vector`| a zero-energy 4-vector|
|[Q1,...,QM]`::Array{Lattice3Vector,N}`| `M` zero-energy 4-vectors from `M` `Lattice3Vector`s in an `Array` of dimension `N`|
|[Q1,...,QM]`::Array{Lattice3Vector,N}`,[E1,...,EM]`::Array{Number,N}`| `M` 4-vectors from `M` `Lattice3Vector`s and `M` energies in equivalently-shaped `Array`s of dimension `N`.|
|[Q1[:],...,QM[:]]`::Array{Number,2}`,lat`::Lattice`| `M` 4-vectors from a `4×M` `Julia` builtin `Matrix{Number}`.|
"""
#type Lattice4Vector{T<:Number} <: LatticeVector
immutable Lattice4Vector{T<:Number} <: LatticeVector
    h::T
    k::T
    l::T
    E::T
    lat::Lattice
end
Lattice4Vector(a1::Number,a2::Number,a3::Number,a4::Number,a5::Lattice)=((a1,a2,a3,a4)=promote(a1,a2,a3,a4);Lattice4Vector(a1,a2,a3,a4,a5))
Lattice4Vector{S<:Number,R<:Number}(Q::Lattice3Vector{S},E::R)=(T=promote_type(S,R); Lattice4Vector{T}(T(Q.h),T(Q.k),T(Q.l),T(E),Q.lat))

Lattice4Vector{T<:Number}(v::AbstractArray{T,1},l::Lattice)=(@assert length(v)==4;Lattice4Vector(v...,l))
Lattice4Vector(h::Number,k::Number,l::Number,lat::Lattice)=Lattice4Vector(h,k,l,0,lat)
Lattice4Vector(lat::Lattice)=Lattice4Vector(0,0,0,0,lat)
Lattice4Vector(Q::Lattice3Vector)=Lattice4Vector(Q,0)
Lattice4Vector{T<:Lattice3Vector}(Q::Array{T})=map(Lattice4Vector,Q)
Lattice4Vector{T<:Lattice3Vector,R<:Number}(Q::Array{T},E::Array{R})=(@assert compatible(Q,E);map(Lattice4Vector,Q,E))

function Lattice4Vector{T<:Number}(m::AbstractArray{T,2},lat::Lattice=one(Reciprocal{T}))
    @assert size(m,1)==4
    b=size(m,2)
    out=Array{Lattice4Vector{T}}(b)
    for i=1:b; out[i]=Lattice4Vector(m[1,i],m[2,i],m[3,i],m[4,i],lat); end
    out
end
function Lattice4Vector{T<:Number}(m::AbstractArray{T,3},lat::Lattice=one(Reciprocal{T}))
    @assert size(m,1)==4
    b=size(m,2); c=size(m,3)
    out=Array{Lattice4Vector{T}}(b,c)
    @inbounds for j=1:c
        @inbounds for i=1:b
            out[i,j]=Lattice4Vector{T}(m[1,i,j],m[2,i,j],m[3,i,j],m[4,i,j],lat)
            # m[:,i,j]... instead of m[1,i,j],m[2,i,j],m[3,i,j],m[4,i,j] is 17x slower and takes 10x the memory!
        end
    end
    return out
end

function showLattice4Vector(io::IO,a::Lattice4Vector,compact::Bool=false)
    Base.showcompact(io,get3vector(a)) # Always show the compact version
    Base.showcompact(io,a.E)
    compact? Base.print(io,"meV") : Base.print(io," meV")
end
Base.show(io::IO,QE::Lattice4Vector)=showLattice4Vector(io,QE,false)
Base.showcompact(io::IO,QE::Lattice4Vector)=showLattice4Vector(io,QE,true)
#Base.zero{T<:Number}(::Type{Lattice4Vector{T}})=Lattice4Vector{T}(zero(Lattice3Vector{T}),zero(T))
Base.zero{T<:Number}(::Type{Lattice4Vector{T}})=Lattice4Vector(zero(T),zero(T),zero(T),zero(T),zero(Reciprocal{T}))
(==)(a::Lattice4Vector,b::Lattice4Vector)=(f=fieldnames(a); all(map(x->a.(x),f).==map(x->b(x),f)))
# get3vector(a::Lattice4Vector) could be replaced by a.lat on the next five lines,
# but since the definition of Lattice4Vector is in flux, this is safer for now

isReciprocal(a::Lattice4Vector)=isReciprocal(a.lat)
isDirect(a::Lattice4Vector)=isDirect(a.lat)

dot(a::Lattice3Vector,b::Lattice4Vector)=dot(a,get3vector(b))
dot(a::Lattice4Vector,b::Lattice3Vector)=dot(get3vector(a),b)
dot(a::Lattice4Vector,b::Lattice4Vector)=dot(get3vector(a),get3vector(b))+a.E*b.E

function Base.norm(x::Lattice4Vector)
    q,E=dot(get3vector(x),get3vector(x)),x.E^2
    (q>0&&E>0)&&(warn("The length of a full 4-vector is ill defined and this result mixes units. ($q,$E)"))
    sqrt(q+E)
end
Base.norm{L<:Lattice4Vector}(v::Array{L})=map(Base.norm,v)

# rounding a Lattice4Vector returns the near-by vector with its indicies rounded to the nearest integer
#round{T<:Number}(::Type{T},a::Lattice4Vector)=Lattice4Vector(round(T,a.h),round(T,a.k),round(T,a.l),round(T,a.E),getlat(a))
#round(a::Lattice4Vector)=Lattice4Vector(round(a.Q),round(a.E))
round{T<:Real}(::Type{T},a::Lattice4Vector)=Lattice4Vector(round(T,a.h),round(T,a.k),round(T,a.l),round(T,a.E),a.lat)
round(a::Lattice4Vector)=Lattice4Vector(round(a.h),round(a.k),round(a.l),round(a.E),a.lat)


#Base.getindex(a::Lattice4Vector,s::Symbol)= any(s.==fieldnames(a))?a.(s):a.Q[s]
#geth(a::Lattice4Vector)=geth(a.Q)
#getk(a::Lattice4Vector)=getk(a.Q)
#getl(a::Lattice4Vector)=getl(a.Q)
#getQ(a::Lattice4Vector)=a.Q
#getE(a::Lattice4Vector)=a.E
#getlat(a::Lattice4Vector)=getlat(a.Q)

#get3vector(a::Lattice4Vector)=a.Q
#get4vector(a::Lattice4Vector)=a

Base.getindex(a::Lattice4Vector,s::Symbol)= any(s.==fieldnames(a))?a.(s):a.lat[s]
geth(a::Lattice4Vector)=a.h
getk(a::Lattice4Vector)=a.k
getl(a::Lattice4Vector)=a.l
getE(a::Lattice4Vector)=a.E
getlat(a::Lattice4Vector)=a.lat
get3vector(a::Lattice4Vector)=Lattice3Vector(a.h,a.k,a.l,a.lat)
get4vector(a::Lattice4Vector)=a
