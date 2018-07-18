# Make life easier with a LatticeVector abstract type
#abstract LatticeVector <: Number # moved to common Abstract.jl file
# watch https://github.com/SimonDanisch/FixedSizeArrays.jl for a possible way to implement general N-D LatticeVectors
# TODO general N-D latticevectors, latticepoints, and latticetensors implemented in Symmetry.jl, making this module obsolete
# TODO Move all dependent modules (just NXS, I think) to Symmetry.jl's new LQVec, LDVec, etc.

# include LatticeVectorN definitions
include("LatticeVector3.jl")
include("LatticeVector4.jl")

Base.hash(a::Type{LatticeVector},h::UInt=0)=for f in fieldnames(a); h=Base.hash(a.(f),h); end

# Define general LatticeVector constructors
LatticeVector(h::Number,k::Number,l::Number,lat::Lattice)=Lattice3Vector(promote(h,k,l,lat)...)
LatticeVector(h::Number,k::Number,l::Number,E::Number,lat::Lattice)=Lattice4Vector(promote(h,k,l,E,lat)...)
LatticeVector{T<:Number}(v::AbstractArray{T,1},lat::Lattice)=LatticeVector(v[:]...,lat) # unwrap the array and let the general constructor handle it.
#function LatticeVector{T<:Number}(v::AbstractArray{T,2},lat::Lattice)
#    (a,b)=size(v)
#    if a==3
#        out=Array(Lattice3Vector{T},b)
#        for i=1:b; out[i]=Lattice3Vector{T}(v[:,i]...,lat); end
#    elseif a==4
#        out=Array(Lattice4Vector{T},b)
#        for i=1:b; out[i]=Lattice4Vector{T}(v[:,i]...,lat); end
#    end
#    out
#end
function LatticeVector{T<:Number}(v::AbstractArray{T,2},lat::Lattice)
    if size(v,1)==3
        out=Lattice3Vector(v,lat)
    elseif size(v,1)==4
        out=Lattice4Vector(v,lat)
    end
    out
end
function LatticeVector{T<:Number}(v::AbstractArray{T},lat::Lattice)
    n=ndims(v)
    out=LatticeVector(slicedim(v,n,1),lat)
    for i=2:size(v,n); out=cat(n-1,out,LatticeVector(slicedim(v,n,i),lat)); end
    out
end
LatticeVector(lat::Lattice)=Lattice3Vector(lat)


direct_or_reciprocal(a::Lattice)=typeof(a)<:Direct?Direct:Reciprocal
direct_or_reciprocal(a::LatticeVector)=direct_or_reciprocal(a.lat)

geth{T<:LatticeVector}(a::Array{T,1})=(z=zeros(typeof(geth(a[1])),length(a));    for i=1:length(a);                      z[i]=geth(a[i]); end; z)
geth{T<:LatticeVector}(a::Array{T,2})=(z=zeros(typeof(geth(a[1])),size(a)...);   for j=1:size(a,2); for i=1:size(a,1); z[i,j]=geth(a[i,j]); end; end; z)
geth{T<:LatticeVector,N}(a::Array{T,N})=(z=zeros(typeof(geth(a[1])),size(a)...); for i=1:length(a);                      z[i]=geth(a[i]); end; z)
getk{T<:LatticeVector}(a::Array{T,1})=(z=zeros(typeof(geth(a[1])),length(a));    for i=1:length(a);                      z[i]=getk(a[i]); end; z)
getk{T<:LatticeVector}(a::Array{T,2})=(z=zeros(typeof(geth(a[1])),size(a)...);   for j=1:size(a,2); for i=1:size(a,1); z[i,j]=getk(a[i,j]); end; end; z)
getk{T<:LatticeVector,N}(a::Array{T,N})=(z=zeros(typeof(geth(a[1])),size(a)...); for i=1:length(a);                      z[i]=getk(a[i]); end; z)
getl{T<:LatticeVector}(a::Array{T,1})=(z=zeros(typeof(geth(a[1])),length(a));    for i=1:length(a);                      z[i]=getl(a[i]); end; z)
getl{T<:LatticeVector}(a::Array{T,2})=(z=zeros(typeof(geth(a[1])),size(a)...);   for j=1:size(a,2); for i=1:size(a,1); z[i,j]=getl(a[i,j]); end; end; z)
getl{T<:LatticeVector,N}(a::Array{T,N})=(z=zeros(typeof(geth(a[1])),size(a)...); for i=1:length(a);                      z[i]=getl(a[i]); end; z)
getE{T<:LatticeVector}(a::Array{T,1})=(z=zeros(typeof(geth(a[1])),length(a));    for i=1:length(a);                      z[i]=getE(a[i]); end; z)
getE{T<:LatticeVector}(a::Array{T,2})=(z=zeros(typeof(geth(a[1])),size(a)...);   for j=1:size(a,2); for i=1:size(a,1); z[i,j]=getE(a[i,j]); end; end; z)
getE{T<:LatticeVector,N}(a::Array{T,N})=(z=zeros(typeof(geth(a[1])),size(a)...); for i=1:length(a);                      z[i]=getE(a[i]); end; z)
getlat{T<:LatticeVector}(v::Array{T,1})=(rl=[getlat(v[1])]; for i=2:length(v); push!(rl,getlat(v[i])); end; rl)
getlat{T<:LatticeVector}(v::Array{T,2})=(rl=getlat(v[:,1]); for i=2:size(v,2); rl=hcat(rl,getlat(v[:,i])); end; rl)
#getlat{T<::LatticeVector}(v::Array{T,3})= #== what's the 3-Array version of hcat? ==#

gethkl(a::LatticeVector)=[geth(a),getk(a),getl(a)]
gethklE(a::LatticeVector)=[geth(a),getk(a),getl(a),getE(a)]
gethkl{T<:LatticeVector,N}(a::Array{T,N})=permutedims(cat(N+1,geth(a),getk(a),getl(a)),circshift(1:N+1,1))
gethklE{T<:LatticeVector,N}(a::Array{T,N})=permutedims(cat(N+1,geth(a),getk(a),getl(a),getE(a)),circshift(1:N+1,1))

round{T<:LatticeVector}(a::Array{T})=map(round,a)
round{T<:Real,R<:LatticeVector}(::Type{T},a::Array{R})=map(x->round(T,x),a)

get3vector{T<:LatticeVector}(a::Array{T})=map(get3vector,a)
get4vector{T<:LatticeVector}(a::Array{T})=map(get4vector,a)

# all *,/,+,- operator overloading needs to be in one place to sort out conflicting definitions

# unitary minus
(-)(a::Lattice3Vector)=Lattice3Vector(-a.h,-a.k,-a.l,a.lat)
(-)(a::Lattice4Vector)=Lattice4Vector(-a.h,-a.k,-a.l,-a.E,a.lat)

# addition
(+)(x::Lattice3Vector,y::Lattice3Vector)=(@assert sameLattice(x,y); Lattice3Vector(x.h+y.h,x.k+y.k,x.l+y.l,x.lat))
(+)(a::Lattice4Vector,b::Lattice4Vector)=(@assert sameLattice(a,b); Lattice4Vector(a.h+b.h,a.k+b.k,a.l+b.l,a.E+b.E,a.lat))

# subtraction
(-)(x::Lattice3Vector,y::Lattice3Vector)=(@assert sameLattice(x,y); Lattice3Vector(x.h-y.h,x.k-y.k,x.l-y.l,x.lat))
(-)(a::Lattice4Vector,b::Lattice4Vector)=(@assert sameLattice(a,b); Lattice4Vector(a.h-b.h,a.k-b.k,a.l-b.l,a.E-b.E,a.lat))

# multiplication
(*)(x::Lattice3Vector,y::Lattice3Vector)= dot(x,y)
(*)(x::Lattice4Vector,y::Lattice4Vector)= dot(x,y)
(*){T<:LatticeVector,R<:LatticeVector}(::T,::R)=throw(InexactError())
(*)(x::Lattice3Vector,y::Number)=Lattice3Vector(x.h*y,x.k*y,x.l*y,x.lat)
(*)(x::Lattice4Vector,y::Number)=Lattice4Vector(x.h*y,x.k*y,x.l*y,x.E*y,x.lat)
(*){T<:LatticeVector}(::Bool,::T)=throw(InexactError())
(*){T<:LatticeVector}(y::Number,x::T)=x*y


# division
(/)(::Lattice3Vector,::Lattice3Vector)=throw(InexactError())
(/)(::Lattice4Vector,::Lattice4Vector)=throw(InexactError())
(/)(::Lattice3Vector,::Lattice4Vector)=throw(InexactError())
(/)(::Lattice4Vector,::Lattice3Vector)=throw(InexactError())
(/)(x::Lattice3Vector,y::Number)=Lattice3Vector(x.h/y,x.k/y,x.l/y,x.lat)
(/)(y::Number,x::Lattice3Vector)=Lattice3Vector(y/x.h,y/x.k,y/x.l,x.lat)
(/)(x::Lattice4Vector,y::Number)=Lattice4Vector(x.h/y,x.k/y,x.l/y,x.E/y,x.lat)
(/)(y::Number,x::Lattice4Vector)=Lattice4Vector(y/x.h,y/x.k,y/x.l,y/x.E,x.lat)
#(/){T<:Number}(x::Lattice3Vector,y::Array{T})=map(z->x/z,y)
#(/){T<:Number}(y::Array{T},x::Lattice3Vector)=map(z->z/x,y)

# FIXME is overloading * of two arrays of LatticeVectors for the dot product between elements necessary?
broadcast{T<:LatticeVector,L<:LatticeVector,N}(::typeof(*),a::Array{T,N},b::Array{L,N})=broadcast(dot,a,b)

(^)(a::LatticeVector,n::Rational)=a^(n.num/n.den)
(^)(a::LatticeVector,n::Integer)=Base.norm(a)^n
(^)(a::LatticeVector,n::Number)= (n==0?1.:(n<0?1/a^abs(n):Base.norm(a)^n))


function discernable{L<:LatticeVector}(v::Array{L},tol=1e-4)
    d=trues(size(v)...)
    for j=length(v):-1:2
        any(norm(v[1:j-1].-v[j]).<tol) && (d[j]=false)
    end
    return d
end
dunique{L<:LatticeVector}(v::Array{L},tol=1e-4)=v[discernable(v)]

# define dot product between vectors of LatticeVectors
dot{T<:LatticeVector}(a::Array{T,1},b::Array{T,1})=map(dot,a,b)
dot{T<:LatticeVector,S<:LatticeVector,N}(a::Array{T,N},b::Array{S,N})=map(dot,a,b)



Base.convert{T<:Number}(::Type{Lattice4Vector{T}},a::Lattice3Vector{T})=get4vector(a)
Base.convert{T<:Number}(::Type{Lattice3Vector{T}},a::Lattice4Vector{T})=get3vector(a)
#convert{T<:Number}(::Type{Lattice4Vector{T}},a::Number)=Lattice4Vector(a) # this won't work because the lattice is unknown
Base.promote_rule{T<:Number}(::Type{Lattice4Vector{T}},::Type{Lattice3Vector{T}})=Lattice4Vector{T}


#Base.promote_rule{T<:Number,R<:Number}(::Lattice4Vector{T},::Lattice4Vector{R})=Lattice4Vector{Base.promote_rule{T,R}}
Base.promote_type{T<:Number,R<:Number}(::Type{Lattice4Vector{T}},::Type{Lattice4Vector{R}})=Lattice4Vector{Base.promote_type(T,R)}
Base.convert{T<:Number}(::Type{Lattice4Vector{T}},a::Lattice4Vector{T})=a # no conversion necessary if R==T
Base.convert{T<:Number,R<:Number}(::Type{Lattice4Vector{T}},a::Lattice4Vector{R})=Lattice4Vector{T}(T(a.h),T(a.k),T(a.l),T(a.E),direct_or_reciprocal(a){T}(a.lat))
#Base.promote_type{T<:Number,R<:Number}(::Type{Lattice3Vector{T}},::Type{Lattice3Vector{R}})=Lattice3Vector{Base.promote_type{T,R}}
#Base.promote_rule{T<:Number,R<:Number}(::Type{Lattice3Vector{T}},::Type{Lattice3Vector{R}})=Lattice3Vector{Base.promote_rule{T,R}}
Base.promote_type{T<:Number,R<:Number}(::Lattice3Vector{T},::Lattice3Vector{R})=Lattice3Vector{Base.promote_type{T,R}}
Base.convert{T<:Number}(::Type{Lattice3Vector{T}},a::Lattice3Vector{T})=a # no conversion necessary if R==T
Base.convert{T<:Number,R<:Number}(::Type{Lattice3Vector{T}},a::Lattice3Vector{R})=Lattice3Vector{T}(T(a.h),T(a.k),T(a.l),direct_or_reciprocal(a){T}(a.lat))

# XXX These should be redundant with the implicit broadcasting in v0.6, e.g. "dot.(a,b)" becomes "broadcast(dot,a,b)
# define special dot functions between arrays of vectors, including an effective .dot scaling-up version
#Base.dot{T<:LatticeVector,N}(a::Array{T,N},b::Array{T,N})=map(Base.dot,a,b)
#Base.dot{T<:LatticeVector,N}(a::Array{T,N},b::T)=map(x->Base.dot(x,b),a)
#Base.dot{T<:LatticeVector,N}(a::T,b::Array{T,N})=map(x->Base.dot(a,x),b)

#Base.promote_op{T<:Number,R<:Number}(::Any,::Lattice3Vector{T},::Lattice3Vector{R})=Lattice3Vector{Base.promote_type{T,R}}


# it's dumb that I have to define these, but something is wrong with the julia reduce.jl variants:
function Base.sum{T<:LatticeVector}(v::Vector{T})
    @assert length(v)>0
    s=0*v[1]
    for i=1:length(v); s+=v[i]; end
    return s
end
function Base.sum{T<:LatticeVector}(v::Matrix{T},n)
    1==n && (v=v'; n=2)
    @assert 2==n "reduction can only take place along 1st or 2nd dimension for a matrix"
    @assert size(v,1)>0 && size(v,2)>0
    s=0*v[:,1]
    for i=1:size(v,1); s[i]=sum(v[i,:]); end
    return s
end
Base.sum{T<:LatticeVector}(v::Matrix{T})=Base.sum(Base.sum(v,2))
Base.mean{T<:LatticeVector}(v::Array{T})=Base.sum(v)/length(v) # good for matricies and vectors
function Base.mean{T<:LatticeVector}(m::Matrix{T},n)
    1==n && (m=m'; n=2)
    @assert 2==n "reduction can only take place along 1st or 2nd dimension for a matrix"
    v=sum(m,n)/size(m,n)
    return v
end
