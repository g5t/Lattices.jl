# Types that have a .lat field: Lattice3Vector, Lattice4Vector, Lattice3Tensor, Lattice4Tensor
export sameLattice

# Two lattices are the same if isapprox returns true:
sameLattice(a::Lattice,b::Lattice)=Base.isapprox(a,b)

# Some objects directly contain a lattice field named "lat"
const WithLatticeField = Union{LatticeVector,LatticeTensor}

# An object that contains a lattice can be checked against a given lattice:
sameLattice(a::WithLatticeField,b::Lattice)=sameLattice(a.lat,b)
sameLattice(a::Lattice,b::WithLatticeField)=sameLattice(a,b.lat)
# Two objects that contain lattices can be checked against each other
sameLattice(a::WithLatticeField,b::WithLatticeField)=sameLattice(a.lat,b.lat)

# The crystal object contains two lattices (which are the reciprocal of each other)
# Two crystals are of the sameLattice if both of their direct and reciprocal lattices are the same
sameLattice(a::Crystal,b::Crystal)=(sameLattice(a.dlat,b.dlat)&sameLattice(a.rlat,b.rlat))
# Checking for sameness with other Lattices/Lattice-containing-objects is ambiguous due to the duality
sameLattice(a::Crystal,b::Union{Lattice,WithLatticeField})=(sameLattice(a.dlat,b)|sameLattice(a.rlat,b))
sameLattice(b::Union{Lattice,WithLatticeField},a::Crystal)=(sameLattice(b,a.dlat)|sameLattice(b,a.rlat))

# Sample objects contain a Crystal object and two Sample objects are of the same lattice if their crystals are of the same lattice
sameLattice(a::Sample,b::Sample)=sameLattice(a.crystal,b.crystal)
#
sameLattice(a::Sample,b::Union{Lattice,WithLatticeField,Crystal})=sameLattice(a.crystal,b)
sameLattice(a::Union{Lattice,WithLatticeField,Crystal},b::Sample)=sameLattice(a,b.crystal)



# also define invLattices
export invLattices
invLattices(d::Direct,r::Reciprocal)=sameLattice(d,star(r))||sameLattice(star(d),r)
invLattices(r::Reciprocal,d::Direct)=sameLattice(r,star(d))||sameLattice(star(r),d)
invLattices(::Lattice,::Lattice)=false # should catch only same lattice types

invLattices(a::WithLatticeField,b::Lattice)=invLattices(a.lat,b)
invLattices(a::Lattice,b::WithLatticeField)=invLattices(a,b.lat)
invLattices(a::WithLatticeField,b::WithLatticeField)=invLattices(a.lat,b.lat)
