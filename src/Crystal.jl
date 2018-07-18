"""
A type for any crystaline sample, comprised of its `Direct` and `Reciprocal` lattice .

Additional optional parameters are:

| keyword |          |
|:-------:|:---------|
| `group`      | the symmetry group name to which this crystal belongs                                         |
| `groupno`    | the symmetry group number to which this crystal belongs                                       |
| `positions`  | a `Direct` lattice `Lattice3Vector` array of the positions of all  atoms in the unit cell     |
| `ions`       | a list of ion names, e.g., 'Fe2+', for form factor calculations                               |
| `types`      | a list of integers to indicate equivalent ions for exchanges                                  |
| `moments`    | a `Direct` lattice `Lattice3Vector` array of the local-moments at each of the atom positions  |
| `mu`         | a unit `Lattice3Vector` defining the local z-direction for spinwave calculations              |
| `anisotropy` | a list of single-ion anisotropy for each local-moment                                         |
| `jions`      | an Nx2 list of atoms involved in each exchange, where the values are indicies into `positions`|
| `jvectors`   | `Direct` lattice `Lattice3Vector` array of extra vectors for long-range exchange interactions |
| `jexchanges` | a list of the values of the exchange energies (in meV)                                        |

all fields with integer elements are `UInt16` (with max value 0xffff=65535)
except for `groupno` which is `UInt8` (with max value 0xff=255 > max group number 230)
"""
type Crystal{T<:Real}
    dlat::Direct{T}
    rlat::Reciprocal{T}
# It might be useful to include additional crystal parameters here
    group::AbstractString
    groupno::UInt8
    pos::Array{Lattice3Vector{T},1}
    occ::Array{Float64,1}  # occupancy of the site, (potentially) usefull for substituted samples
    ion::Array{AbstractString,1} # ion labels for form factor calculations
    typ::Array{UInt16,1}        # ion type, for identifying equivalent sites
    spin::Array{Lattice3Vector{T},1}
    mu::Lattice3Vector{T} # the local z-direction for spinwave calculations
    sia::Array{T,1} # single ion anisotropy
    sig::Array{T,1} # single ion g-factor
    jion::Array{UInt16,2} # Nx2 list of atoms involved in each exchange interaction
    jvec::Array{Lattice3Vector{T},1} # extra vector between atoms with the given exchange
    jexc::Array{T,1}
    jexn::Array{UInt16,1} # equal-distance-exchange number, closest-exchanges all have 1,
    jcel::UInt8 # Number of unit cells to include in exchange determination
    mindist::Float64
    maxdist::Float64
    maxbonds::UInt16
    function Crystal{R}(d::Direct{R},r::Reciprocal{R},
                        g::AbstractString,gn::UInt8,
                        p::Array{Lattice3Vector{R},1},
                        oc::Array{Float64,1},
                        i::Array{AbstractString,1},t::Array{UInt16,1},
                        m::Array{Lattice3Vector{R},1},µ::Lattice3Vector{R},
                        a::Array{R,1},b::Array{R,1},o::Array{UInt16,2},v::Array{Lattice3Vector{R},1},
                        x::Array{R,1},n::Array{UInt16,1},c::UInt8,
                        mind::Float64,maxd::Float64,maxb::UInt16) where R<:Real
        @assert length(p)==length(oc)==length(i)==length(t)==length(m)==length(a)==length(b)
        @assert size(o,1)==length(v)==length(x)==length(n)
        @assert size(o,2)==2
        for z in p; @assert sameLattice(d,z); end
        for z in m; @assert sameLattice(d,z); end
        for z in v; @assert sameLattice(d,z); end
        µ/=norm(µ) # ensure µ is a unit vector
        new(d,r,g,gn,p,oc,i,t,m,µ,a,b,o,v,x,n,c,mind,maxd,maxb)
    end
end
function Crystal{R<:Real}(d::Direct{R},r::Reciprocal{R};
                          groupstr::AbstractString="",
                          groupno::UInt8=0x00,
                          positions::Array{Lattice3Vector{R},1}=zeros(Lattice3Vector{R},0),
                          occupancy::Array{Float64,1}=ones(Float64,length(positions)),
                          ions::Array{AbstractString,1}=Array{AbstractString}(length(positions)),
                          types::Array{UInt16,1}=zeros(UInt16,length(positions)),
                          moments::Array{Lattice3Vector{R},1}=length(positions)>0?positions*0:zeros(eltype(positions),0),
                          µ::Lattice3Vector{R}=Lattice3Vector{R}(0,0,1.,d),
                          anisotropy::Array{R,1}=zeros(R,length(positions)),
                          gfactor::Array{R,1}=2ones(R,length(positions)),
                          jions::Array{UInt16,2}=zeros(UInt16,0,2),
                          jvectors::Array{Lattice3Vector{R},1}=zeros(eltype(positions),size(jions,1)),
                          jexchanges::Array{R,1}=ones(R,size(jions,1)),
                          jexnumber::Array{UInt16,1}=zeros(UInt16,size(jions,1)),
                          jcells::UInt8=0x01,
                          mindist::Float64=0.,
                          maxdist::Float64=Inf,
                          maxbonds::UInt16=0xffff
                          )
    (groupstr,groupno)=isempty(groupstr)?group_string_number(groupno):group_string_number(groupstr)
    @assert length(positions)==length(occupancy)==length(ions)==length(types)==length(moments)==length(anisotropy)==length(gfactor)
    @assert size(jions,1)==length(jvectors)==length(jexchanges)==length(jexnumber)
    @assert size(jions,2)==2
    Crystal{R}(d,r,groupstr,groupno,positions,occupancy,ions,types,moments,µ,anisotropy,gfactor,jions,jvectors,jexchanges,jexnumber,jcells,mindist,maxdist,maxbonds)
end
"""
    Crystal(a,b,c,α,β,γ)
Initialize a `Crystal` object by specifying the lattice parameters of its direct lattice.
All inputs have default values, with lattice angles expected in degrees (then converted to radian).
The unit cell lengths default to 2π, and the angles default to 90° (π/2); such that the components
of vectors in the reduced lattice units of default direct and reciprocal lattices are in absolute
units (Å and Å⁻¹, respectively).
"""
function Crystal(a::Real=2pi,b::Real=2pi,c::Real=2pi,α::Real=90,β::Real=90,γ::Real=90)
    # it's nice to enter degrees, but calculations are easier if angles are in radians
    α,β,γ=α/180*pi,β/180*pi,γ/180*pi
    lat=Direct(promote(a,b,c,α,β,γ)...)
    Crystal(lat,star(lat))
end

Crystal{T<:Real}(r::Reciprocal{T})=Crystal{T}(star(r),r)
Crystal{T<:Real}(d::Direct{T})=Crystal{T}(d,star(d))
Base.show(io::IO,m::Crystal)=Base.show(io,m.dlat);
Base.showcompact(io::IO,m::Crystal)=Base.showcompact(io,m.dlat);

(==)(a::Crystal,b::Crystal)=all(map(x->getfield(a,x)==getfield(b,x),fieldnames(a))) # two crystals with identical elements are equivalent

Base.getindex(a::Crystal,s::Symbol)=any(fieldnames(a).==s)?a.(s):error("Unknown fieldname $s in $(typeof(a))")
