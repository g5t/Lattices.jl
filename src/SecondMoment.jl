module SecondMoments
import Base: *,/,^,+,-
export SecondMoment
"""
The second moment of any object is strictly a Tensor property, with
<xᵢxⱼ>= (∭ xᵢ xⱼ P(r⃗)dr⃗)/(∭P(r⃗)dr⃗)
where P(r⃗) is the probability of finding part of the object at r⃗
macroscopic objects have P(r⃗)={0,1} depending on their shape, which is used, in practice,
to set limits on the infinite integrals (note that the denominator is just the length/area/volume
of the object, depending on the dimensionality of r⃗)
The immutable type `SecondMoment` should hold only <xᵢxⱼ> values; for use with Multiple-Dispatch
function overloading.
"""
immutable SecondMoment{T<:Real}
    sm::T
end
Base.zero{T<:Real}(::Type{SecondMoment{T}})=Base.zero(T)
Base.one{T<:Real}(::Type{SecondMoment{T}})=Base.one(T)
Base.zero(::Type{SecondMoment})=Base.zero(Real)
Base.one(::Type{SecondMoment})=Base.one(Real)
Base.zero(a::SecondMoment)=Base.zero(Base.typeof(a))
Base.one(a::SecondMoment)=Base.one(Base.typeof(a))

Base.convert{T<:Real}(::Type{Bool},x::SecondMoment{T})=Base.convert(Bool,x.sm)
Base.convert{T<:Real}(::Type{Integer},x::SecondMoment{T})=Base.convert(Integer,x.sm)
Base.convert{T<:Real}(::Type{Rational{Base.GMP.BigInt}},x::SecondMoment{T})=Base.convert(Rational{GMP.BigInt},x.sm)
Base.convert{T<:Integer,S<:Real}(::Type{Rational{T}},x::SecondMoment{S})=Base.convert(Rational{T},x.sm)
Base.convert{S<:Real,T<:Real}(::Type{S},x::SecondMoment{T})=Base.convert(S,x.sm)

#Base.convert{T<:Real}(::Type{SecondMoment{T}},x::Real)=SecondMoment(x)
Base.convert{T<:Real}(::Type{SecondMoment{T}},x::Real)=SecondMoment(Base.convert(T,x))
Base.convert(::Type{SecondMoment},x::Real)=SecondMoment(Base.convert(Real,x)) # allows for, e.g., SecondMoment[1,2,3,4]

Base.promote_rule{S<:Real,T<:Real}(::Type{SecondMoment{S}},::Type{SecondMoment{T}})=Base.promote_rule(S,T)
Base.promote_rule{S<:Real,T<:Real}(::Type{SecondMoment{S}},::Type{T})=Base.promote_rule(S,T)
Base.show{T<:Real}(io::IO,x::SecondMoment{T})=Base.show(io,x.sm)
Base.showcompact{T<:Real}(io::IO,x::SecondMoment{T})=Base.showcompact(io,x.sm)

for x in (:*,:/,:+,:-)
    @eval $x(a::SecondMoment,b::SecondMoment)=SecondMoment($x(a.sm,b.sm)) # preserve second-momentness
    @eval $x(a::SecondMoment,b::Real)=$x(a.sm,b)
    @eval $x(a::Real,b::SecondMoment)=$x(a,b.sm)
end
for x in (:Integer,:Real)
    @eval (^)(a::SecondMoment,b::$x)=SecondMoment(a.sm^b)
end
(^)(a::SecondMoment,b::SecondMoment)=SecondMoment(a.sm^b.sm)

end# module
