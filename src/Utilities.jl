# getBytes modified from http://stackoverflow.com/questions/28402989/check-size-in-bytes-of-variable-using-julia
getBytes(x::DataType) = sizeof(x);
function getBytes(x)
   total = 0;
   fieldNames = fieldnames(typeof(x));
   if fieldNames == []
      return sizeof(x);
   else
     for fieldName in fieldNames
        if !isa(getfield(x,fieldName),Function)
            total += getBytes(getfield(x,fieldName));
        end
     end
     return total;
   end
end

_multMv{T<:Real}(M::Array{T,2},v::Array{T})=(Mv=similar(v);_multMv!(M,v,Mv);Mv)
_multMv!{T<:Real}(M::Array{T,2},v::Array{T})=_multMv!(M,v,v)
# Matrix*Vector multiplication is already defined, but this variant of _multMv performs the same error checking and allows for less complex code elsewhere
function _multMv!{T<:Real}(M::Array{T,2},v::Array{T,1},Mv::Array{T,1})
    a=length(v)
    @assert all([size(M)...].==[a,a]) "M must be square and v must be consistently sized [size(v,1)=size(M,1)=size(M,2)]"
    @assert length(Mv)==a "Mv must be a vector the same size as v"
    Mv[:]=M*v # [:] is important to assign result of M*v to passed-in Mv
end
function _multMv!{T<:Real}(M::Array{T,2},v::Array{T,2},Mv::Array{T,2})
    (a,b)=size(v)
    @assert all([size(M)...].==[a,a]) "M must be square and v must be consistently sized [size(v,1)=size(M,1)=size(M,2)]"
    @assert all([size(Mv)...].==[a,b]) "Mv must be an array of vectors the same size as v"
    for i=1:b; Mv[:,i]=M*v[:,i]; end
end
function _multMv!{T<:Real}(M::Array{T,2},v::Array{T,3},Mv::Array{T,3})
    (a,b,c)=size(v)
    @assert all([size(M)...].==[a,a]) "M must be square and v must be consistently sized [size(v,1)=size(M,1)=size(M,2)]"
    @assert all([size(Mv)...].==[a,b,c]) "Mv must be an array of vectors the same size as v"
    for j=1:c; for i=1:b; Mv[:,i,j]=M*v[:,i,j]; end; end
end

function _multMv{T<:Real,R<:Real}(M::Array{T,3},v::Array{R,1})
    (a,b,c)=size(M)
    @assert length(v)==b
    Mv=zeros(Base.promote_type(T,R),a,c)
    _unsafe_multMv!(M,v,Mv,a=a,c=c)
    return Mv
end
function _multMv!{T<:Real,R<:Real,S<:Real}(M::Array{T,3},v::Array{R,1},Mv::Array{S,2})
    @assert S<:Base.promote_type(T,R)
    (a,b,c)=size(M)
    @assert length(v)==b
    @assert all([size(Mv)...].==[a,c])
    _unsafe_multMv!(M,v,Mv;a=a,c=c)
end
function _unsafe_multMv!{T<:Real,R<:Real,S<:Real}(M::Array{T,3},v::Array{R,1},Mv::Array{S,2};a=size(M,1),c=size(M,3))
    @inbounds for j=1:c;
        @inbounds for i=1:a;
            @inbounds Mv[i,j]=sum(vec(M[i,:,j]).*v);
        end;
    end
end

function _multMv{T<:Real,R<:Real}(M::Array{T,3},v::Array{R,2})
    (a,b,c)=size(M)
    @assert all([size(v)...].==[b,c])
    Mv=zeros(Base.promote_type(T,R),a,c)
    for j=1:c
        for i=1:a
            Mv[i,j]=sum(vec(M[i,:,j]).*v[:,j])
        end
    end
    return Mv
end
function _multMv3{T<:Real,R<:Real}(M::Array{T,3},v::Array{R,2})
    (a,b,c)=size(M)
    (vb,d)=size(v)
    @assert vb==b
    Mv=zeros(Base.promote_type(T,R),a,c,d)
    for k=1:d
        for j=1:c
            for i=1:a
                Mv[i,j,k]=sum(vec(M[i,:,j]).*v[:,k])
            end
        end
    end
    return Mv
end


function _multMv_LV{T<:Real}(M::Array{T,2},v::Array{T,1},l::Lattice)
    @assert all(length(v)*[1,1].==[size(M)...])
    LatticeVector(M*v,l)
end
function _multMv_LV{T<:Real}(M::Array{T,2},v::Array{T},l::Lattice)
    sz=[size(v)...]
    @assert (sz[1]==3)||(sz[1]==4)
    lt=sz[1]==3?Lattice3Vector:Lattice4Vector
    Mv=length(sz)>1?Array{lt{T}}(sz[2:end]...):zero(lt{T})
    _multMv_LV!(M,v,Mv,l)
    Mv
end
function _multMv_LV!{T<:Real,L<:LatticeVector}(M::Array{T,2},v::Array{T,2},Mv::Array{L,1},l::Lattice)
    (a,b)=size(v)
    @assert all([size(M)...].==[a,a]) "M must be square and v must be consistently sized [size(v,1)=size(M,1)=size(M,2)]"
    @assert length(Mv)==b "Mv must be a vector the same length as size(v,2)"
    for i=1:b; Mv[i]=LatticeVector(M*v[:,i],l); end
end
function _multMv_LV!{T<:Real,L<:LatticeVector}(M::Array{T,2},v::Array{T,3},Mv::Array{L,2},l::Lattice)
    (a,b,c)=size(v)
    @assert all([size(M)...].==[a,a]) "M must be square and v must be consistently sized [size(v,1)=size(M,1)=size(M,2)]"
    @assert all([size(Mv)...].==[b,c]) "Mv must be an array of vectors the same size as v (without its first dimension)"
    for j=1:c; for i=1:b; Mv[i,j]=LatticeVector(M*v[:,i,j],l); end; end
end

function _multMv_LV{T<:Real,R<:Real}(M::Array{T,3},v::Array{R,2},l::Lattice)
    (M1,M2,M3)=size(M)
    (v1,v2)=size(v)
    @assert all([M1,M2].==[v1,v1]) "M must be square and v must be consistently sized [size(v,1)=size(M,1)=size(M,2)]"
    lt=M1==3?Lattice3Vector:Lattice4Vector
    Mv=Array(lt{Base.promote_type(T,R)},M3,v2)
    _unsafe_multMv_LV!(M,v,Mv,l;v2=v2,M3=M3)
    return Mv
end
function _multMv_LV!{T<:Real,R<:Real,L<:LatticeVector}(M::Array{T,3},v::Array{R,2},Mv::Array{L,2},l::Lattice)
    (v1,v2)=size(v)
    (M1,M2,M3)=size(M)
    (Mv1,Mv2)=size(Mv)
    @assert all([M1,M2].==[v1,v1]) "M must be square and v must be consistently sized [size(v,1)=size(M,1)=size(M,2)]"
    @assert all([Mv1,Mv2].==[M3,v2]) "Mv must be size(M,3) by size(v,2)"
    _unsafe_multMv_LV!(M,v,Mv,l;v2=v2,M3=M3)
end
function _unsafe_multMv_LV!{T<:Real,L<:LatticeVector}(M::Array{T,3},v::Array{T,2},Mv::Array{L,2},l::Lattice;v2=size(v,2),M3=size(M,3))
    for j=1:v2
        for i=1:M3
            Mv[i,j]=LatticeVector(M[:,:,i]*v[:,j],l)
        end
    end
end


# these rotation matricies are rotations that change the direction of a vector within the x-y plane
vectorrotation3xy(θ::Number)=typeof(θ)[cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]
vectorrotation4xy(θ::Number)=typeof(θ)[cos(θ) -sin(θ) 0 0; sin(θ) cos(θ) 0 0; 0 0 1 0; 0 0 0 1]
# these rotation matricies are basis-changes that keep the direction of the vector unchanged while the coordinate system is rotated
basisrotation3xy(θ::Number)=typeof(θ)[cos(θ) sin(θ) 0; -sin(θ) cos(θ) 0; 0 0 1]
basisrotation4xy(θ::Number)=typeof(θ)[cos(θ) sin(θ) 0 0; -sin(θ) cos(θ) 0 0; 0 0 1 0; 0 0 0 1]

function cubepoints(n::Integer)
   cl=2n+1
   np=cl^3
   p1=zeros(typeof(n),np); p2=copy(p1); p3=copy(p1);
   for i=1:cl
   for j=1:cl
   for k=1:cl
       idx=(i-1)*cl^2+(j-1)*cl+k;
       p1[idx]=i-1-n; p2[idx]=j-1-n; p3[idx]=k-1-n;
   end
   end
   end
    (p1,p2,p3)
end
halfcubepoints(n::Unsigned)=halfcubepoints(signed(n))
function halfcubepoints{T<:Signed}(n::T)
    np=fld((2n+1)^3-1,2) # ((2n+1)^3-1)/2 preserving integer-ness
    p1=zeros(T,np); p2=copy(p1); p3=copy(p1);
    # for i=0, j=0, k>0
    p3[1:n]=1:n;
    idx=n+one(T);
    for j=one(T):n # i=0, j>0
    for k=-n:n
        p2[idx]=j; p3[idx]=k; idx+=one(T)
    end
    end
    for i=one(T):n # i>0
    for j=-n:n
    for k=-n:n
        p1[idx]=i; p2[idx]=j; p3[idx]=k; idx+=one(T);
    end
    end
    end
    (p1,p2,p3)
end



function elementsplit(aei::String)
    relement=r"[a-z,A-Z]+"
    s=split(aei,relement)
    (s[1],aei[search(aei,relement)],s[end]) #(Atomic Number, Element symbol, Ion specification)
end
