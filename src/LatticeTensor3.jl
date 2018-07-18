"""
A `Lattice3Tensor` is a (position or momentum) 3×3 rank-2 Tensor in (`Direct` or `Reciprocal`) lattice units.
"""
immutable Lattice3Tensor{T<:Real} <: LatticeTensor
    hh::T;hk::T;hl::T
    kh::T;kk::T;kl::T
    lh::T;lk::T;ll::T
    lat::Lattice
end
Lattice3Tensor(xx,xy,xz,yx,yy,yz,zx,zy,zz,lat::Lattice)=Lattice3Tensor(promote(xx,xy,xz,yx,yy,yz,zx,zy,zz)...,lat)
Lattice3Tensor(lat::Lattice)=Lattice3Tensor(1,0,0,0,1,0,0,0,1.,lat)
## Add initializer that takes a 3×3 Array and a Lattice
Lattice3Tensor{T<:Real}(v::AbstractArray{T,2},l::Lattice=one(Reciprocal{T}))=(@assert size(v)==(3,3);Lattice3Tensor(v...,l))
function Lattice3Tensor{T<:Real}(m::AbstractArray{T,3},l::Lattice=one(Reciprocal{T}))
    @assert size(m,1,2)==(3,3)
    out=Array{Lattice3Tensor{T}}(size(m,3))
    for i=1:size(m,3);out[i]=Lattice3Tensor{T}(m[:,:,i]...,l);end
    out
end
function Lattice3Tensor{T<:Real}(m::AbstractArray{T,4},l::Lattice=one(Reciprocal{T}))
    (a,b,c,d)=size(m)
    @assert size(m,1,2)==(3,3)
    out=Array{Lattice3Tensor{T}}(c,d)
    @inbounds for j=1:d
        @inbounds for i=1:c
            out[i,j]=Lattice3Tensor{T}(m[:,:,i,j]...,l);
        end
    end
    return out
end

Base.zero{T<:Real}(::Type{Lattice3Tensor{T}},l::Lattice=zero(Reciprocal{T}))=Lattice3Tensor{T}(zeros(T,9),l)
Base.zero{T<:Real}(a::Lattice3Tensor{T})=Base.zero(typeof(a),a.lat)
(==)(a::Lattice3Tensor,b::Lattice3Tensor)=(norm(a-b)==0)
#Base.hash{T<:Real}(a::Type{Lattice3Tensor{T}},h::UInt=0)=for f in fieldnames(a); h=Base.hash(a.(f),h); end


Base.isapprox(a::Lattice3Tensor,b::Lattice3Tensor)=isapprox(norm(a-b),0)



isReciprocal(a::Lattice3Tensor)=isa(a.lat,Reciprocal)
isDirect(a::Lattice3Tensor)=isa(a.lat,Direct)

matrixform(a::Lattice3Tensor)=[a.hh a.hk a.hl; a.kh a.kk a.kl; a.lh a.lk a.ll]
# Return the Miller indexes of a Tensor defined by its fractional cell coordinates, or vice versa.
function star(x::Lattice3Tensor) # Tᵢⱼ=gᵢᵧTᵞᵝgᵦⱼ/4π²
    cmt=covariantMetricTensor(x.lat)
    Lattice3Tensor(cmt*(matrixform(x)*cmt)/4/pi/pi,star(x.lat))
end
"""
The matrix product between a `Lattice3Tensor` and a `Lattice3Vector` depends on their order
and the orientation of the vector. Since the `Lattices` implementation of vectors does not
distinguish between row and column vectors the matrix product assumes the appropriate orientation.
"""
function (*)(m::Lattice3Tensor,v::Lattice3Vector)
    @assert sameLattice(m,v)
    lat=m.lat;abc=1./[lat[:a],lat[:b],lat[:c]]; ten=abc*abc' # this probably wrong for non-orthogonal lattices, and specialized to dot(v,m*v)
    h=dot(Lattice3Vector([m.hh,m.hk,m.hl].*ten[:,1],lat),v)
    k=dot(Lattice3Vector([m.kh,m.kk,m.kl].*ten[:,2],lat),v)
    l=dot(Lattice3Vector([m.lh,m.lk,m.ll].*ten[:,3],lat),v)
    Lattice3Vector(h,k,l,lat)
end
function (*)(v::Lattice3Vector,m::Lattice3Tensor)
    @assert sameLattice(m,v)
    lat=m.lat; abc=1./[lat[:a],lat[:b],lat[:c]]; ten=abc*abc'  # this probably wrong for non-orthogonal lattices, and specialized to dot(v*m,v)
    h=dot(v,Lattice3Vector([m.hh,m.kh,m.lh].*ten[1,:],lat))
    k=dot(v,Lattice3Vector([m.hk,m.kk,m.lk].*ten[2,:],lat))
    l=dot(v,Lattice3Vector([m.hl,m.kl,m.ll].*ten[3,:],lat))
    Lattice3Vector(h,k,l,lat)
end

function outer(a::Lattice3Vector,b::Lattice3Vector)
    @assert sameLattice(a,b)
    va=gethkl(a); vb=gethkl(b)
    Lattice3Tensor(hcat(va*vb[1],va*vb[2],va*vb[3]),a.lat)
end

Base.norm(x::Lattice3Tensor)=Base.norm(matrixform(x))
Base.norm{L<:Lattice3Tensor}(v::Array{L})=map(Base.norm,v)

# rounding a Lattice3Tensor returns the near-by Tensor with its incicies rounded to the nearest integer
Base.round{T<:Integer}(::Type{T},x::Lattice3Tensor)=Lattice3Tensor(Base.round(T,matrixform(x)),x.lat)
Base.round{T<:AbstractFloat}(::Type{T},x::Lattice3Tensor)=Lattice3Tensor(Base.round(matrixform(x)),x.lat)
Base.round(x::Lattice3Tensor)=Lattice3Tensor(Base.round(matrixform(x)),x.lat)

Base.ceil{T<:Real}(::Type{T},x::Lattice3Tensor)=Lattice3Tensor(Base.ceil(T,matrixform(x)),x.lat)
Base.floor{T<:Real}(::Type{T},x::Lattice3Tensor)=Lattice3Tensor(Base.floor(T,matrixform(x)),x.lat)

#Base.getindex(a::Lattice3Tensor,s::Symbol)= any(s.==fieldnames(Lattice3Tensor))?a.(s):a.lat[s]
Base.getindex(a::Lattice3Tensor,s::Symbol)= any(s.==fieldnames(a))?a.(s):a.lat[s]
getlat(a::Lattice3Tensor)=a.lat
