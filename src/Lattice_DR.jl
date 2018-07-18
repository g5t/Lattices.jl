# Lattice.jl
abstract Lattice
"""
The parameters that describe a `Direct` or `Recriprocal` `Lattice`.
Two lattices with the same lattice parameters are, by definition, the same lattice.
Therefore the `LatticeParameters` type is immutable with fields:

| Fieldname | Unit cell | `Direct` units | `Reciprocal` units |
|:---------:|:------------|:-------------|:-----------------|
| `a`,`b`,`c` | lengths | Å | Å⁻¹ |
| `α`,`β`,`γ` | angles | radian | radian |
| `V`         | volume | Å³ | Å⁻³ |
"""
immutable LatticeParameters{T<:Real}
    a::T # Å (or Å⁻¹)
    b::T # Å (or Å⁻¹)
    c::T # Å (or Å⁻¹)
    α::T # radians
    β::T # radians
    γ::T # radians
    V::T # Å³ (or Å⁻³)
end
LatticeParameters(a::Real,b::Real,c::Real,α::Real,β::Real,γ::Real,V::Real)=LatticeParameters(promote(a,b,c,α,β,γ,V)...)
LatticeParameters(a::Real,b::Real,c::Real,α::Real=pi/2,β::Real=pi/2,γ::Real=pi/2)=LatticeParameters(a,b,c,α,β,γ,latticeVolume(a,b,c,α,β,γ))
""" A type to indicate its contained `LatticeParameters` represent a `Direct` `Lattice` """
immutable Direct{T<:Real}     <: Lattice; lattice::LatticeParameters{T}; end
""" A type to indicate its contained `LatticeParameters` represent a `Reciprocal` `Lattice` """
immutable Reciprocal{T<:Real} <: Lattice; lattice::LatticeParameters{T}; end

Base.getindex(a::LatticeParameters,s::Symbol)=a.(s)
Base.getindex(a::Lattice,s::Symbol)= (s==:lattice)?(a.lattice):(a.lattice.(s)) # works because Direct and Reciprocal have only one field

function showLattice(io::IO,m::Lattice,compact::Bool=false)
        abc=[m.lattice.a,m.lattice.b,m.lattice.c]
        αβγ=[m.lattice.α,m.lattice.β,m.lattice.γ]/pi*180
        compact?
            (all(map(x->(isa(x,Measured)?isapprox(x,Measured(pi,1e-16)):isapprox(x,pi)),abc))?
            Base.print(io,"[2π,2π,2π]") : Base.showcompact(io,abc)) : Base.show(io,abc)
        isa(m,Reciprocal) ? Base.print(io,"Å⁻¹ ") : Base.print(io,"Å ")
        compact? Base.showcompact(io,αβγ) : Base.show(io,αβγ)
        Base.print(io,"° ")
end
Base.show(io::IO,m::Lattice)=showLattice(io,m,false)
Base.showcompact(io::IO,m::Lattice)=showLattice(io,m,true)
Base.one{T<:Real}(::Type{Reciprocal{T}})=Reciprocal(ones(T,3)...)
Base.one{T<:Real}(::Type{Direct{T}})=Direct(ones(T,3)...)
Base.one(::Type{Reciprocal})=Reciprocal(ones(3)...) # default to Float64
Base.one(::Type{Direct})=Direct(ones(3)...)
Base.zero{T<:Real}(::Type{Reciprocal{T}})=Reciprocal(zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T))
Base.zero{T<:Real}(::Type{Direct{T}})=Direct(zero(T),zero(T),zero(T),zero(T),zero(T),zero(T),zero(T))
Base.isapprox(a::LatticeParameters,b::LatticeParameters)=all([isapprox(a.(fn),b.(fn)) for fn in fieldnames(a)])
Base.isapprox{T<:Lattice}(a::T,b::T)=Base.isapprox(a.lattice,b.lattice) # should catch typeof(a,b) both Reciprocal or both Direct
Base.isapprox(::Lattice,::Lattice)=false # should only catch mixed lattice types
sameLattice(a::Lattice,b::Lattice)=Base.isapprox(a,b)
invLattices(d::Direct,r::Reciprocal)=sameLattice(d,star(r))||sameLattice(star(d),r)
invLattices(r::Reciprocal,d::Direct)=sameLattice(r,star(d))||sameLattice(star(r),d)
invLattices(::Lattice,::Lattice)=false # should catch only same lattice types

# Calculate the volume of a lattice defined by a,b,c,α,β,γ. Taken from star.m of ResLib 3.4a by A. Zheludev
latticeVolume{T<:AbstractFloat}(a::T,b::T,c::T,α::T,β::T,γ::T)=2*a*b*c*sqrt(sin((α+β+γ)/2)*sin((-α+β+γ)/2)*sin((α-β+γ)/2)*sin((α+β-γ)/2))
# And use that definition to calculate the volume of a Lattice
latticeVolume(l::LatticeParameters)=latticeVolume(l.a,l.b,l.c,l.α,l.β,l.γ)
latticeVolume(a::Lattice)=latticeVolume(a.lattice)

# Calculate the reciprocal lattice parameters from a,b,c,α,β,γ,V. Taken from star.m of ResLib 3.4a by A. Zheludev
function star(a::Real,b::Real,c::Real,α::Real,β::Real,γ::Real,V::Real)
    Ta=typeof(a);Tb=typeof(b);Tc=typeof(c);Tα=typeof(α);Tβ=typeof(β);Tγ=typeof(γ);TV=typeof(V)
    astar=Ta(2*pi*b*c*sin(α)/V)
    bstar=Tb(2*pi*c*a*sin(β)/V)
    cstar=Tc(2*pi*a*b*sin(γ)/V)
    αstar=Tα(acos((cos(β)*cos(γ)-cos(α))/(sin(β)*sin(γ))))
    βstar=Tβ(acos((cos(γ)*cos(α)-cos(β))/(sin(γ)*sin(α))))
    γstar=Tγ(acos((cos(α)*cos(β)-cos(γ))/(sin(α)*sin(β))))
    Vstar=TV(8*pi^3/V)
    (astar,bstar,cstar,αstar,βstar,γstar,Vstar)
end
# Use that definition to define the star of the lattice, l, which is another lattice, l⋆.
star(l::LatticeParameters)=LatticeParameters(star(l.a,l.b,l.c,l.α,l.β,l.γ,l.V)...)
star(l::Direct)=Reciprocal(star(l.lattice))
star(l::Reciprocal)=Direct(star(l.lattice))

# Calculate the metric tensor for a lattice.
metrictensor(a::Real,b::Real,c::Real,α::Real,β::Real,γ::Real) = [a^2 a*b*cos(γ) a*c*cos(β); b*a*cos(γ) b^2 b*c*cos(α); c*a*cos(β) c*b*cos(α) c^2]
metrictensor(l::LatticeParameters)=metrictensor(l.a,l.b,l.c,l.α,l.β,l.γ)
metrictensor(l::Lattice)=metrictensor(l.lattice)

covariantMetricTensor(l::Lattice)=metrictensor(l) # gᵢⱼ. Can be used to convert a vector in the opposite basis to l into l's basis: aᵢ=2π gᵢⱼ aʲ
contravariantMetricTensor(l::Lattice)=inv(metrictensor(l)) # gⁱʲ. Can be used to convert a vector in the same basis as l into the opposite basis: aⁱ=2π gⁱʲ aⱼ

Base.convert{T<:Real,R<:Real}(::Type{LatticeParameters{T}},a::LatticeParameters{R})=LatticeParameters{T}(T(a.a),T(a.b),T(a.c),T(a.α),T(a.β),T(a.γ),T(a.V))

for x in (:Direct,:Reciprocal)
    @eval $x(a::Real,b::Real,c::Real,α::Real,β::Real,γ::Real,V::Real)=$x(LatticeParameters(a,b,c,α,β,γ,V))
#     Lattice initialization without specification of volume; defaulting to α,β,γ=90° 
    @eval $x(a::Real,b::Real,c::Real,α::Real=pi/2,β::Real=pi/2,γ::Real=pi/2)=$x(a,b,c,α,β,γ,latticeVolume(a,b,c,α,β,γ))
    @eval Base.convert{T<:Real,R<:Real}(::Type{$x{T}},a::$x{R})=$x{T}(T(a.lattice))
end

BMatrix(a::Real,b::Real,c::Real,α::Real,β::Real,γ::Real,V::Real)=BMatrix(Direct(a,b,c,α,β,γ,V))
BMatrix(d::Direct)=(r=star(d); BMatrix(d,r))
BMatrix(r::Reciprocal)=(d=star(r); BMatrix(d,r))
function BMatrix(x::Direct,y::Reciprocal)
    d=x.lattice; r=y.lattice; 
    # B-matix as in Acta Cryst. (1967). 22, 457 http://dx.doi.org/10.1107/S0365110X67000970
    B=[r.a r.b*cos(r.γ) r.c*cos(r.β); 0 r.a*abs(sin(r.γ)) -r.c*abs(sin(r.β))*cos(d.α); 0 0 2*pi/d.c]
end
BMatrix(r::Reciprocal,d::Direct)=BMatrix(d,r)



