# Lattice.jl
"""
An abstract type with subtypes `Direct` and `Reciprocal`.
"""
abstract type Lattice end
"""
The parameters that describe a `Direct` `Lattice`.
Two lattices with the same lattice parameters are, by definition, the same lattice.
Therefore the `Direct` type is immutable with fields:

| Fieldname | Unit cell | units |
|:---------:|:------------|:-------------|
| `a`,`b`,`c` | lengths | Å |
| `α`,`β`,`γ` | angles | radian |
| `V`         | volume | Å³ |
"""
immutable Direct{T<:Real} <: Lattice
    a::T # Å
    b::T # Å
    c::T # Å
    α::T # radians
    β::T # radians
    γ::T # radians
    V::T # Å³
end
"""
The parameters that describe a `Recriprocal` `Lattice`.
Two lattices with the same lattice parameters are, by definition, the same lattice.
Therefore the `Reciprocal` type is immutable with fields:

| Fieldname | Unit cell | units |
|:---------:|:------------|:-----------------|
| `a`,`b`,`c` | lengths | Å⁻¹ |
| `α`,`β`,`γ` | angles | radian |
| `V`         | volume | Å⁻³ |
"""
immutable Reciprocal{T<:Real} <: Lattice
    a::T # Å⁻¹
    b::T # Å⁻¹
    c::T # Å⁻¹
    α::T # radians
    β::T # radians
    γ::T # radians
    V::T # Å⁻³
end
for x in (:Direct,:Reciprocal)
    @eval $x(a::Real,b::Real,c::Real,α::Real,β::Real,γ::Real,V::Real)=$x(promote(a,b,c,α,β,γ,V)...)
#     Lattice initialization without specification of volume; defaulting to α,β,γ=90°
    @eval $x{T<:Real}(a::T,b::Real,c::Real,α::Real=T(pi/2),β::Real=T(pi/2),γ::Real=T(pi/2))=$x(a,b,c,α,β,γ,T(latticeVolume(a,b,c,α,β,γ)))
    @eval Base.convert{T<:Real}(::Type{$x{T}},a::$x{T})=a # no need to convert types if T==T
    @eval Base.convert{T<:Real,R<:Real}(::Type{$x{T}},a::$x{R})=$x{T}(map(z->Base.convert(T,z),[a.a,a.b,a.c,a.α,a.β,a.γ,a.V])...)
    @eval Base.one{T<:Real}(::Type{$x{T}})=$x(ones(T,3)...)
    @eval Base.one(::Type{$x})=$x(ones(3)...)
    @eval Base.zero{T<:Real}(::Type{$x{T}})=$x(zeros(T,3)...)
    @eval Base.zero(::Type{$x})=$x(zeros(3)...)
end
# Direct{R} to Direct{T} and Reciprocal{R} to Reciprocal{T} already defined above
for (x,y) in ((:Direct,:Reciprocal),(:Reciprocal,:Direct))
    @eval Base.convert{T<:Real,R<:Real}(::Type{$x{T}},a::$y{R})=Base.convert($x{T},star(a))
end


Base.getindex(a::Lattice,s::Symbol)= any(fieldnames(a).==s) ? getfield(a,s) : error("Unknown field $s for $(typeof(a))")

isnpi{T<:Real}(x::T)=Base.abs(x/pi-Base.round(x/pi))/(x/pi) < 1e-14
#isnpi{T<:Real}(x::Array{T})=map(isnpi,x)

function showLattice(io::IO,m::Lattice,compact::Bool=false)
        abc=[m.a,m.b,m.c]
        αβγ=[m.α,m.β,m.γ]#/pi*180
        αβγ = all(isnpi.(2αβγ)) ? Base.round.(Integer,αβγ/pi*180) : αβγ/pi*180
        if all(isnpi.(abc))
            n=Base.round.(Integer,abc/pi)
            q=["$(x)π" for x in n]
            Base.print(io,"[",q[1],",",q[2],",",q[3],"]")
        else
            compact ? Base.showcompact(io,abc) : Base.show(io,abc)
        end
        isa(m,Reciprocal) ? Base.print(io,"Å⁻¹ ") : Base.print(io,"Å ")
        compact? Base.showcompact(io,αβγ) : Base.show(io,αβγ)
        Base.print(io,"° ")
end
Base.show(io::IO,m::Lattice)=showLattice(io,m,false)
Base.showcompact(io::IO,m::Lattice)=showLattice(io,m,true)

# used by sameLattice
# Base.isapprox{T<:Lattice}(a::T,b::T)=all([testfield(Base.isapprox,fn,a,b) for fn in fieldnames(a)])
# Base.isapprox(::Lattice,::Lattice)=false # should only catch mixed lattice types

Base.isapprox(l1::Lattice,l2::Lattice)=all(testfield.(Base.isapprox,fieldnames(l1),l1,l2))

# Calculate the volume of a lattice defined by a,b,c,α,β,γ. Taken from star.m of ResLib 3.4a by A. Zheludev
#latticeVolume{T<:AbstractFloat}(a::T,b::T,c::T,α::T,β::T,γ::T)=2*a*b*c*sqrt(sin((α+β+γ)/2)*sin((-α+β+γ)/2)*sin((α-β+γ)/2)*sin((α+β-γ)/2))
latticeVolume(a::Real,b::Real,c::Real,α::Real,β::Real,γ::Real)=2*a*b*c*sqrt(sin((α+β+γ)/2)*sin((-α+β+γ)/2)*sin((α-β+γ)/2)*sin((α+β-γ)/2))
# And use that definition to calculate the volume of a Lattice
latticeVolume(l::Lattice)=latticeVolume(l.a,l.b,l.c,l.α,l.β,l.γ)

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
star(l::Direct)=Reciprocal(star(l.a,l.b,l.c,l.α,l.β,l.γ,l.V)...)
star(l::Reciprocal)=Direct(star(l.a,l.b,l.c,l.α,l.β,l.γ,l.V)...)

# Calculate the metric tensor for a lattice.
metrictensor(a::Real,b::Real,c::Real,α::Real,β::Real,γ::Real) = [a^2 a*b*cos(γ) a*c*cos(β); b*a*cos(γ) b^2 b*c*cos(α); c*a*cos(β) c*b*cos(α) c^2]
metrictensor(l::Lattice)=metrictensor(l.a,l.b,l.c,l.α,l.β,l.γ)

covariantMetricTensor(l::Lattice)=metrictensor(l) # gᵢⱼ. Can be used to convert a vector in the opposite basis to l into l's basis: aᵢ=2π gᵢⱼ aʲ
contravariantMetricTensor(l::Lattice)=inv(metrictensor(l)) # gⁱʲ. Can be used to convert a vector in the same basis as l into the opposite basis: aⁱ=2π gⁱʲ aⱼ




BMatrix(a::Real,b::Real,c::Real,α::Real,β::Real,γ::Real,V::Real)=BMatrix(Direct(a,b,c,α,β,γ,V))
BMatrix(d::Direct)=(r=star(d); BMatrix(d,r))
BMatrix(r::Reciprocal)=(d=star(r); BMatrix(d,r))
function BMatrix(x::Direct,y::Reciprocal)
    d=x.lattice; r=y.lattice;
    # B-matix as in Acta Cryst. (1967). 22, 457 http://dx.doi.org/10.1107/S0365110X67000970
    B=[r.a r.b*cos(r.γ) r.c*cos(r.β); 0 r.a*abs(sin(r.γ)) -r.c*abs(sin(r.β))*cos(d.α); 0 0 2*pi/d.c]
end
BMatrix(r::Reciprocal,d::Direct)=BMatrix(d,r)
