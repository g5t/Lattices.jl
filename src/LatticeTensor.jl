# Make life easier with a LatticeTensor abstract type
#abstract LatticeTensor <: Real # moved to common Abstract.jl file
# watch https://github.com/SimonDanisch/FixedSizeArrays.jl for a possible way to implement general N-D LatticeTensors

# include LatticeTensorN definitions
include("LatticeTensor3.jl")
include("LatticeTensor4.jl")

Base.hash(a::Type{LatticeTensor},h::UInt=0)=for f in fieldnames(a); h=Base.hash(a.(f),h); end

Base.diagm(a::Lattice3Vector)=Lattice3Tensor(diagm(gethkl(a)),a.lat)
Base.diagm(a::Lattice4Vector)=Lattice4Tensor(diagm(gethklE(a)),a.lat)

# Define general LatticeTensor constructors
function LatticeTensor{T<:Real}(mat::Array{T},lat::Lattice)
    nel=prod([size(mat,1,2)...])
    if nel==9
        out=Lattice3Tensor(mat,lat)
    elseif nel==16
        out=Lattice4Tensor(mat,lat)
    else
        error("Only 3- and 4-Tensors are supported")
    end
    out
end
direct_or_reciprocal(a::LatticeTensor)=direct_or_reciprocal(a.lat)

function showLatticeTensor(io::IO,m::LatticeTensor,compact::Bool=false)
    !compact && (Base.println(io,typeof(m))) 
    Base.showcompact(io,matrixform(m))
    isa(m.lat,Reciprocal)? Base.print(io,"ʳ") : Base.print(io,"ᵈ") # reciprocal (ᵣ) or direct (ᵈ) lattice units)
end
Base.show(io::IO,m::LatticeTensor) = showLatticeTensor(io,m,false)
Base.showcompact(io::IO,m::LatticeTensor)=showLatticeTensor(io,m,true)

round{T<:LatticeTensor}(a::Array{T})=map(round,a)
round{T<:Real,R<:LatticeTensor}(::Type{T},a::Array{R})=map(x->round(T,x),a)

# all *,/,+,- operator overloading needs to be in one place to sort out conflicting definitions

# unitary minus
(-)(a::Lattice3Tensor)=Lattice3Tensor(-matrixform(a),a.lat)
(-)(a::Lattice4Tensor)=Lattice4Tensor(-matrixform(a),a.lat)

# addition
(+)(x::Lattice3Tensor,y::Lattice3Tensor)=(@assert sameLattice(x,y); Lattice3Tensor(matrixform(x)+matrixform(y),x.lat))
(+)(a::Lattice4Tensor,b::Lattice4Tensor)=(@assert sameLattice(a,b); Lattice4Tensor(matrixform(a)+matrixform(b),a.lat))

# subtraction
(-)(x::Lattice3Tensor,y::Lattice3Tensor)=(@assert sameLattice(x,y); Lattice3Tensor(matrixform(x)-matrixform(y),x.lat))
(-)(a::Lattice4Tensor,b::Lattice4Tensor)=(@assert sameLattice(a,b); Lattice4Tensor(matrixform(a)-matrixform(b),a.lat))

# multiplication
(*)(x::Lattice3Tensor,y::Lattice3Tensor)=(@assert sameLattice(x,y); Lattice3Tensor(matrixform(x)*matrixform(y),x.lat))
(*)(x::Lattice4Tensor,y::Lattice4Tensor)=(@assert sameLattice(x,y); Lattice4Tensor(matrixform(x)*matrixform(y),x.lat))
(*){T<:LatticeTensor,R<:LatticeTensor}(::T,::R)=throw(InexactError())
(*)(x::Lattice3Tensor,y::Real)=Lattice3Tensor(matrixform(x)*y,x.lat) 
(*)(x::Lattice4Tensor,y::Real)=Lattice4Tensor(matrixform(x)*y,x.lat)
(*){T<:LatticeTensor}(::Bool,::T)=throw(InexactError())
(*){T<:LatticeTensor}(y::Real,x::T)=x*y


# division
(/)(::Lattice3Tensor,::Lattice3Tensor)=throw(InexactError())
(/)(::Lattice4Tensor,::Lattice4Tensor)=throw(InexactError())
(/)(::Lattice3Tensor,::Lattice4Tensor)=throw(InexactError())
(/)(::Lattice4Tensor,::Lattice3Tensor)=throw(InexactError())

(/)(::Lattice3Tensor,::Lattice3Vector)=throw(InexactError())
(/)(::Lattice4Tensor,::Lattice4Vector)=throw(InexactError())
(/)(::Lattice3Tensor,::Lattice4Vector)=throw(InexactError())
(/)(::Lattice4Tensor,::Lattice3Vector)=throw(InexactError())

(/)(::Lattice3Vector,::Lattice3Tensor)=throw(InexactError())
(/)(::Lattice4Vector,::Lattice4Tensor)=throw(InexactError())
(/)(::Lattice3Vector,::Lattice4Tensor)=throw(InexactError())
(/)(::Lattice4Vector,::Lattice3Tensor)=throw(InexactError())

(/)(x::Lattice3Tensor,y::Real)=Lattice3Tensor(matrixform(x)./y,x.lat)
(/)(x::Lattice4Tensor,y::Real)=Lattice4Tensor(matrixform(x)./y,x.lat)
(/)(y::Real,x::Lattice3Tensor)=Lattice3Tensor(y./matrixform(x),x.lat)
(/)(y::Real,x::Lattice4Tensor)=Lattice4Tensor(y./matrixform(x),x.lat)

#(*){T<:Real}(x::LatticeTensor,y::Array{T})=map(z->x*z,y)
#(*){T<:Real}(y::Array{T},x::LatticeTensor)=map(z->z*x,y)
#(.*){T<:Real}(a::LatticeTensor,b::Array{T})=map(x->(a*x),b)
#(.*){T<:Real}(a::Array{T},b::LatticeTensor)=map(x->(x*b),a)
#(./){T<:Real}(a::LatticeTensor,b::Array{T})=map(x->(a/x),b)
#(./){T<:Real}(a::Array{T},b::LatticeTensor)=map(x->(x/b),a)
##(.*){T<:LatticeTensor,L<:LatticeTensor,N}(a::Array{T,N},b::Array{L,N})=map(Base.(*),a,b)
#(.*){T<:Real,L<:LatticeTensor,N}(a::Array{T,N},b::Array{L,N})=map((x,y)->(x*y),a,b)
#(.*){T<:Real,L<:LatticeTensor,N}(a::Array{L,N},b::Array{T,N})=map((x,y)->(x*y),a,b)
#(./){T<:LatticeTensor,L<:LatticeTensor,N}(a::Array{T,N},b::Array{L,N})=throw(InexactError())
#(./){T<:Real,L<:LatticeTensor,N}(a::Array{T,N},b::Array{L,N})=map((x,y)->(x/y),a,b)
#(./){T<:Real,L<:LatticeTensor,N}(a::Array{L,N},b::Array{T,N})=map((x,y)->(x/y),a,b)

function discernable{L<:LatticeTensor}(v::Array{L},tol=1e-4)
    d=trues(size(v)...)
    for j=length(v):-1:2
        any(norm(v[1:j-1].-v[j]).<tol) && (d[j]=false)
    end
    return d
end
dunique{L<:LatticeTensor}(v::Array{L},tol=1e-4)=v[discernable(v,tol)]
