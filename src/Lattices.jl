module Lattices
include("SecondMoment.jl")
using .SecondMoments # relative submodule
using StaticArrays, FastGaussQuadrature # used in Geometry.jl

importall Measureds # to avoid overwriting discernable

import Base: (==),(+),(-),(/),(*),(./),(.*),(^)
import Base: norm,dot,cross

export Lattice,
       LatticeVector,Lattice3Vector,Lattice4Vector,
       LatticeTensor,Lattice3Tensor,Lattice4Tensor,
       star,roundfinite,
       Crystal,Sample,orthonormalize,reshape!,setshape,resense!,setsense,
       getsense,getshape,getwidth,getheight,
       lab2sample,lab2sample!,sample2lab,getx,gety,getz,getE,geth,getk,getl,
       getlat,get3vector,get4vector, gethkl,gethklE, atompositions
export norm,dot,cross,angle,dunique,discernable
export Direct,Reciprocal,isDirect,isReciprocal,direct_or_reciprocal
export innertype
export _multMv,_multMV!

"""LatticeVector is an abtract type for all vectors given in units of a direct or reciprocal lattice."""
abstract type LatticeVector <: Real end
"""LatticeTensor is an abtract type for all tensors given in units of a direct or reciprocal lattice."""
abstract type LatticeTensor <: Real end

# include("StatusMessages.jl")
global DEBUG=false
global DEBUG_LEVEL = 1
debug(newdebug::Bool=DEBUG)=(global DEBUG=newdebug)
debug_level(newdebuglevel::Integer=DEBUG_LEVEL)=(global DEBUG_LEVEL=newdebuglevel)
function testfield(test,field,a,b)
       res = test(getfield(a,field),getfield(b,field))
       if DEBUG
              if DEBUG_LEVEL > 1
                     info(field," ",res?"ok":"mismatch between $(getfield(a,field)) and $(getfield(b,field))")
              elseif !res
                     info("$field mismatch")
              end
       end
       return res
end
include("Geometry.jl")
include("Lattice_Split.jl")
include("LatticeVector.jl") # internally includes LatticeNVector definitions
include("LatticeTensor.jl") # internally includes LatticeNTensor definitions
include("Crystal.jl") # defines the Crystal type
#include("Sample.jl") # defines the Sample type // uses SecondMoment
include("Sample_geometry.jl")
include("SameLattice.jl") # a single point to define all instances of sameLattice
include("Symmetry.jl")
include("Couplings.jl")
include("Utilities.jl")
include("StructureFactor.jl")
include("FormFactor.jl")
include("Innertype.jl") # defines innertype

include("Atom_Ion_Radii.jl")


end # module
