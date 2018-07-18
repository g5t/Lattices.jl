"""
A `Lattice4Tensor` is a (position or momentum) 4×4 rank-2 Tensor in (`Direct` or `Reciprocal`) lattice units.
"""
immutable Lattice4Tensor{T<:Real} <: LatticeTensor
    hh::T;hk::T;hl::T;he::T
    kh::T;kk::T;kl::T;ke::T
    lh::T;lk::T;ll::T;le::T
    eh::T;ek::T;el::T;ee::T
    lat::Lattice
end
Lattice4Tensor(xx,xy,xz,xw,yx,yy,yz,yw,zx,zy,zz,zw,wx,wy,wz,ww,lat::Lattice)=Lattice4Tensor(promote(xx,xy,xz,xw,yx,yy,yz,yw,zx,zy,zz,zw,wx,wy,wz,ww)...,lat)
Lattice4Tensor(lat::Lattice)=Lattice4Tensor(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1.,lat)
## Add initializer that takes a 3×3 Array and a Lattice
Lattice4Tensor{T<:Real}(v::AbstractArray{T,2},l::Lattice=one(Reciprocal{T}))=(@assert size(v)==(4,4);Lattice4Tensor(v...,l))
function Lattice4Tensor{T<:Real}(m::AbstractArray{T,3},l::Lattice=one(Reciprocal{T}))
    @assert size(m,1,2)==(4,4)
    out=Array{Lattice4Tensor{T}}(size(m,3))
    for i=1:size(m,3);out[i]=Lattice4Tensor{T}(m[:,:,i]...,l);end
    out
end
function Lattice4Tensor{T<:Real}(m::AbstractArray{T,4},l::Lattice=one(Reciprocal{T}))
    (a,b,c,d)=size(m)
    @assert size(m,1,2)==(4,4)
    out=Array{Lattice4Tensor{T}}(c,d)
    @inbounds for j=1:d
        @inbounds for i=1:c
            out[i,j]=Lattice4Tensor{T}(m[:,:,i,j]...,l);
        end
    end
    return out
end

Base.zero{T<:Real}(::Type{Lattice4Tensor{T}},l::Lattice=zero(Reciprocal{T}))=Lattice4Tensor{T}(zeros(T,16),l)
Base.zero{T<:Real}(a::Lattice4Tensor{T})=Base.zero(typeof(a),a.lat)
(==)(a::Lattice4Tensor,b::Lattice4Tensor)=(norm(a-b)==0)
#Base.hash{T<:Real}(a::Type{Lattice4Tensor{T}},h::UInt=0)=for f in fieldnames(a); h=Base.hash(a.(f),h); end


Base.isapprox(a::Lattice4Tensor,b::Lattice4Tensor)=isapprox(norm(a-b),0)

isReciprocal(a::Lattice4Tensor)=isa(a.lat,Reciprocal)
isDirect(a::Lattice4Tensor)=isa(a.lat,Direct)

"""
The matrix product between a `Lattice4Tensor` and a `Lattice4Vector` depends on their order
and the orientation of the vector. Since the `Lattices` implementation of vectors does not
distinguish between row and column vectors the matrix product assumes the appropriate orientation.
"""
function (*)(m::Lattice4Tensor,v::Lattice4Vector)
    @assert sameLattice(m,v)
    lat=m.lat;
    abc=1./[lat[:a],lat[:b],lat[:c],1];
    ten=abc*abc' # this probably wrong for non-orthogonal lattices, and specialized to dot(v,m*v)
    h=dot(Lattice4Vector([m.hh,m.hk,m.hl,m.he].*ten[1,:],lat),v)
    k=dot(Lattice4Vector([m.kh,m.kk,m.kl,m.ke].*ten[2,:],lat),v)
    l=dot(Lattice4Vector([m.lh,m.lk,m.ll,m.le].*ten[3,:],lat),v)
    w=dot(Lattice4Vector([m.eh,m.ek,m.el,m.ee].*ten[4,:],lat),v)
    Lattice4Vector(h,k,l,w,lat)
end
function (*)(v::Lattice3Vector,m::Lattice4Tensor)
    @assert sameLattice(m,v)
    lat=m.lat; abc=1./[lat[:a],lat[:b],lat[:c],1];
    ten=abc*abc' # this probably wrong for non-orthogonal lattices, and specialized to dot(v,m*v)
    h=dot(v,Lattice4Vector([m.hh,m.kh,m.lh,m.eh].*ten[:,1],lat))
    k=dot(v,Lattice4Vector([m.hk,m.kk,m.lk,m.ek].*ten[:,2],lat))
    l=dot(v,Lattice4Vector([m.hl,m.kl,m.ll,m.el].*ten[:,3],lat))
    w=dot(v,Lattice4Vector([m.he,m.ke,m.le,m.ee].*ten[:,4],lat))
    Lattice4Vector(h,k,l,w,lat)
end

function outer(a::Lattice4Vector,b::Lattice4Vector)
    @assert sameLattice(a,b)
    va=gethklE(a); vb=gethklE(b)
    Lattice4Tensor(hcat(va*vb[1],va*vb[2],va*vb[3],va*vb[4]),a.lat)
end
matrixform(a::Lattice4Tensor)=[a.hh a.hk a.hl a.he; a.kh a.kk a.kl a.ke; a.lh a.lk a.ll a.le; a.eh a.ek a.el a.ee]

Base.norm(x::Lattice4Tensor)=Base.norm(matrixform(x))
Base.norm{L<:Lattice4Tensor}(v::Array{L})=map(Base.norm,v)

# rounding a Lattice4Tensor returns the near-by Tensor with its incicies rounded to the nearest integer
Base.round{T<:Integer}(::Type{T},x::Lattice4Tensor)=Lattice4Tensor(Base.round(T,matrixform(x)),x.lat)
Base.round{T<:AbstractFloat}(::Type{T},x::Lattice4Tensor)=Lattice4Tensor(Base.round(matrixform(x)),x.lat)
Base.round(x::Lattice4Tensor)=Lattice4Tensor(Base.round(matrixform(x)),x.lat)

Base.ceil{T<:Real}(::Type{T},x::Lattice4Tensor)=Lattice4Tensor(Base.ceil(T,matrixform(x)),x.lat)
Base.floor{T<:Real}(::Type{T},x::Lattice4Tensor)=Lattice4Tensor(Base.floor(T,matrixform(x)),x.lat)

#Base.getindex(a::Lattice4Tensor,s::Symbol)= any(s.==fieldnames(Lattice4Tensor))?a.(s):a.lat[s]
Base.getindex(a::Lattice4Tensor,s::Symbol)= any(s.==fieldnames(a))?a.(s):a.lat[s]
getlat(a::Lattice4Tensor)=a.lat
