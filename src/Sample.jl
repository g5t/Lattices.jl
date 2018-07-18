"""
The `Sample` type represents an element of a spectrometer that is the crystalline sample being studied.
The fields of a `Sample` are:

|Fieldname|Description|
|:-------:|:----------|
|`crystal`    | a `Crystal` object representing the direct |
|             |and reciprocal lattice of the sample |
| `x`,`y`,`z` | orthonormal `Lattice3Vector`s (likely described in |
|             |`crystal.rlat`) defining the orientation of the sample |
| `shape`     | a Tensor of type `Matrix{SecondMoment}` defining the extent|
|             | of the sample in directions of `x`, `y`, and `z`|
|             | with elements `<eᵢeⱼ>`|
| `mosaic`    | an Array of horizontal and vertical mosaic FWHM in radians |
| `sense`     | an `Int8` representing the scattering sense of the sample |
|             |`sense≥0` ↔ to the left, `sense<0` ↔ to the right, from above|
|             | This parameter only applies for triple-axis-type instruments|
"""
type Sample{T<:Real}
    # A sample (can be) comprised of a single crystal, its orienting vectors, its shape, and its scattering sense.
    # Two samples with all fields equivalent are not necessarily the same sample, so type Sample is not immutable.
    # It might be useful to include additional sample parameters here
    crystal::Crystal{T}
    x::Lattice3Vector{T}
    y::Lattice3Vector{T}
    z::Lattice3Vector{T}
    shape::Matrix{SecondMoment{T}}
    mosaic::Array{T,1} # horizontal and vertical mosaic FWHM in radians
    sense::Int8 # scattering sense, +1 or -1 -- could easily be a Bool, but as Int8 the sign is maintained if someone accesses Sample.sense in a calculation
end

setdefault(v,d)=isa(v,Void)?d:v
"""
    Sample(c::Crystal,o1::Vector,o2::Vector,shape::Matrix{SecondMoment},mosaic::Vector,sense::Real)
Initialize a `Sample` object after creating an orthonormal basis from `o1` and `o2`, and converting
`sense` to a `Int8` with only its sign preserved.
"""
function Sample{T<:Real}(c::Crystal{T},
                         o1::Vector=T[1,0,0],
                         o2::Vector=T[0,1,0],
                         shape::Matrix{SecondMoment{T}}=SecondMoment{T}[100/12 0 0; 0 100/12 0; 0 0 100/12],
                         mosaic::Vector{T}=[1,1]/10800*pi,
                         sense::Real=1)
    orient1=Lattice3Vector(convert.(T,o1),c.rlat)::Lattice3Vector{T} # Orienting vectors are typically
    orient2=Lattice3Vector(convert.(T,o2),c.rlat)::Lattice3Vector{T} # defined in reciprocal lattice units
    s=Sample{T}(c,orthonormalize(orient1,orient2)...,shape,mosaic,Int8(signbit(sense)?-1:1))#,sG,sGNo)
end
"""
    Sample(c::Crystal,o1::Vector,o2::Vector,o3::Vector,shape::Matrix{SecondMoment},mosaic::Vector,sense::Real)
Initialize a `Sample` object after creating an orthonormal basis from `o1`,`o2` and `o3`, and converting
`sense` to a `Int8` with only its sign preserved.
"""
function Sample{T<:Real}(c::Crystal{T},o1::Vector,o2::Vector,o3::Vector,
                         shape::Matrix{SecondMoment{T}}=SecondMoment{T}[100/12 0 0; 0 100/12 0; 0 0 100/12],
                         mosaic::Vector{T}=T[1,1]/10800*pi,
                         sense::Real=1)
    orient1=Lattice3Vector(convert.(T,o1),c.rlat)::Lattice3Vector{T} # Orienting vectors are typically
    orient2=Lattice3Vector(convert.(T,o2),c.rlat)::Lattice3Vector{T} # defined in reciprocal lattice units
    orient3=Lattice3Vector(convert.(T,o3),c.rlat)::Lattice3Vector{T}
    s=Sample{T}(c,orthonormalize(orient1,orient2,orient3)...,shape,mosaic,Int8(signbit(sense)?-1:1))#,sG,sGNo)
end
Sample()=Sample(Crystal())

function showSample(io::IO,m::Sample,compact::Bool=false)
    compact && (Base.print(io,signbit(m.sense)?"-":"+"))
    compact ? Base.showcompact(io,m.crystal) : (Base.print(io,"  ");Base.showcompact(io,m.crystal);Base.println(io,""))
    if compact
        Base.showcompact(io,m.x)
        Base.print(io,"—")
        Base.showcompact(io,m.y)
    else
        Base.print(io,"  x̂=")
        Base.showcompact(io,m.x)
        Base.print(io,", ŷ=")
        Base.showcompact(io,m.y)
        Base.print(io,", ẑ=")
        Base.showcompact(io,m.z)
        Base.println(io,"")
    end
    !compact && (Base.println(io,"  sense: ",signbit(m.sense)?"-":"+");)
    !compact && (Base.print(io,"  shape: ");Base.showcompact(io,m.shape))
end
Base.show(io::IO,m::Sample)=showSample(io,m,false)
Base.showcompact(io::IO,m::Sample)=showSample(io,m,true)

(==)(a::Sample,b::Sample)=(a.crystal==b.crystal)&(a.x==b.x)&(a.y==b.y)&(a.z==b.z)&all(a.mosaic.==b.mosaic)&(a.sense==b.sense)
# "all(a.shape.==b.shape)" is causing problems. all(isapprox(a.shape,b.shape)) *might* be more approprite, but leave it out for now

orthonormalize!(a::Sample)= ((a.x,a.y,a.z)=orthonormalize(a.x,a.y))
reshape!{T<:AbstractFloat}(a::Sample{T},s::Matrix{T})=(a.shape=s;)
setshape!{T<:AbstractFloat}(a::Sample{T},s::Matrix{T})=reshape!(a,s)
resense!(a::Sample,s::Real)=(a.sense=signbit(s)?Int8(-1):Int8(1))
setsense!(a::Sample,s::Real)=resense!(a,s)
getsense(a::Sample)=(signbit(a.sense)?-1:1) # return machine-specific Int32 or Int64

"""
    lab2sample!(Q,sample,lv)
Changes the basis of the absolute-unit 3- or 4-vectors described by `Q` in
the laboratory basis to the reciprocal lattice basis of the `sample` storing the resultant
`Lattice3Vector`s or `Lattice4Vectors` in `lv`.

Takes as input an `Array{Real,N}`, `Q`, a triple-axis `Sample` object, and an `Array{LatticeVector,N-1}`, `lv`.
The first dimension of `Q` is taken to represent 3- or 4-vectors and therefore `3 ≤ size(Q,1) ≤ 4`.
The remaining dimensions of `Q` must match the dimensions of `lv` with equivalent sizes.
"""
function lab2sample!{T<:Real,N,L<:LatticeVector,O}(Q::Array{T,N},a::Sample,lv::Array{L,O})
    sQ=[size(Q)...]
    @assert sQ[1]>=3
    @assert ndims(Q)-ndims(lv)==1
    @assert all(sQ[2:end].==[size(lv)...])
    M=eye(typeof(a.x.h),sQ[1])
    M[1:3,1:3]=[ a.x.h a.y.h a.z.h; a.x.k a.y.k a.z.k; a.x.l a.y.l a.z.l ]
    _multMv_LV!(M,Q,lv,a.crystal.rlat)
end
"""
    lab2sample(Q,sample)
Take a vector `Q`, in absolute units (Å⁻¹) defined in the laboratory space,
and return the equivalent `LatticeVector` in the `Sample` reciprocal space
If Q is a 3-vector, assume it is (Qx,Qy,Qz) and return a `Lattice3Vector`
If Q is a 4-vector, assume it is (Qx,Qy,Qz,E) and return a `Lattice4Vector`

If Q is a `Matrix` or higher-order `Array`, take the first dimension to represent
the 3- or 4-vectors and output an Array of `Lattice3Vector`s or `Lattice4Vector`s
with shape equivalent to the remaining dimension(s) of Q.
`Array{Real,N} → Array{LatticeVector,N-1}`
"""
function lab2sample{T<:Real}(Q::Array{T,1},a::Sample)
    M=eye(typeof(a.x.h),size(Q,1))
    M[1:3,1:3]=[ a.x.h a.y.h a.z.h; a.x.k a.y.k a.z.k; a.x.l a.y.l a.z.l ]
    LatticeVector(M*Q,a.crystal.rlat) # pick 3- or 4- vector automatically
end
function lab2sample{T<:Real}(Q::Array{T},a::Sample)
    sz=[size(Q)...]
    @assert (sz[1]==3)||(sz[1]==4)
    lt=sz[1]==3?Lattice3Vector:Lattice4Vector
    lv=Array{lt{T}}(sz[2:end]...)
    lab2sample!(Q,a,lv)
    return lv
end

"""
Take a reciprocal lattice vector, Q, (in rlu) and return the equivalent laboratory space vector (in Å⁻¹)
"""
sample2lab(Q::Lattice3Vector,a::Sample)=[dot(Q,a.x),dot(Q,a.y),dot(Q,a.z)]
sample2lab(Q::Lattice4Vector,a::Sample)=sample2lab(get3vector(Q),a)


# these getx, gety, getz vectorized functions should be combined with those for tripleaxis objects, i.e., T<:Union{Sample,TripleAxis}
# especially if it becomes necessary to ensure the returned Lattic(3)Vectors are promoted to the same type.
getx(a::Sample)=a.x
getx{T<:Sample}(v::Array{T,1})=(rx=[getx(v[1])];  for i=2:length(v);   push!(rx,getx(v[i]));  end; rx)
getx{T<:Sample}(m::Array{T,2})=(rx= getx(v[:,1]); for i=2:size(v,2); rx=hcat(rx,getx(v[:,i])); end; rx)
gety(a::Sample)=a.y
gety{T<:Sample}(v::Array{T,1})=(ry=[gety(v[1])];  for i=2:length(v);   push!(ry,gety(v[i]));  end; ry)
gety{T<:Sample}(m::Array{T,2})=(ry= gety(v[:,1]); for i=2:size(v,2); ry=hcat(ry,gety(v[:,i])); end; ry)
getz(a::Sample)=a.z
getz{T<:Sample}(v::Array{T,1})=(rz=[getz(v[1])];  for i=2:length(v);   push!(rz,getz(v[i]));  end; rz)
getz{T<:Sample}(m::Array{T,2})=(rz= getz(v[:,1]); for i=2:size(v,2); rz=hcat(rz,getz(v[:,i])); end; rz)


getwidth(s::Sample)=sqrt(12*sqrt(det(getshape(s)[1:2,1:2])))
getheight(s::Sample)=sqrt(12*getshape(s)[3,3])
getshape(s::Sample)=deepcopy(s.shape) # second moments too
