generators(g)=decodegenerators(generatorstring(g))
function decodegenerators(g::AbstractString)
    gs=split(g,";"); ndim=length(split(gs[1],","))
    ng=length(gs); r=zeros(Int8,ndim,ndim,ng); t=zeros(Rational{Int8},ndim,ng);
    for i=1:ng; (r[:,:,i],t[:,i])=decodegenerator(gs[i]); end
    return (r,t)
end
function decodegenerator(gs::AbstractString)
    @assert length(findin(gs,";"))==0 "decodegenerator only decodes a single generator, not a list of generators"
    gsi=split(gs,","); ndim=length(gsi);
    @assert ndim==3 "This function can be extended to higher dimensions."
    t=zeros(Rational{Int8},ndim); r=zeros(Int8,ndim,ndim);
    x="[1,0,0]"; y="[0,1,0]"; z="[0,0,1]";
    for j=1:ndim
        if length(findin(gsi[j],"+"))>0
            (gsi[j],rem)=split(gsi[j],"+")
            t[j]=rationalize(eval(parse(rem)))
        end
        r[j,:]=eval(parse(replace(replace(replace(gsi[j],"x",x),"y",y),"z",z)))
    end
    return (r,t)
end
function encodegenerators{T<:Real,R<:Real}(r::Array{T,3},t::Array{R,2})
    S=""
    for i=1:size(r,3); S*=encodegenerator(r[:,:,i],t[:,i])*";"; end
    return S[1:end-1]
end
function encodegenerator{T<:Real,R<:Real}(r::Matrix{R},v::Vector{T}=zeros(R,size(r,1)))
    @assert size(r,1)==length(v)
    S=""
    for i=1:size(r,1)
        xyz=vec(r[i,:])
        s="";
        abs(xyz[1])>0 && (s*=(signbit(xyz[1])?"-":"+")*(abs(xyz[1])==1?"":"$(abs(xyz[1]))")*"x")
        abs(xyz[2])>0 && (s*=(signbit(xyz[2])?"-":"+")*(abs(xyz[2])==1?"":"$(abs(xyz[2]))")*"y")
        abs(xyz[3])>0 && (s*=(signbit(xyz[3])?"-":"+")*(abs(xyz[3])==1?"":"$(abs(xyz[3]))")*"z")
        abs(v[i])>0 && (s*=(signbit(v[i])?"-":"+")*(isa(v[i],Rational)?"$(abs(v[i].num))/$(v[i].den)":"$(abs(v[i]))"))
        s*=","
        S*= s[1]=='+' ? s[2:end] : s
    end
    return S[1:end-1]
end
#function encodegenerator_new{T<:Real,R<:Real}(r::Matrix{R},v::Vector{T}=zeros(R,size(r,1)))
#    (r1,r2)=size(r)
#    @assert r1==r2==length(v)==3
#    sr=map(x->signbit(x)?"-":"+",r).*map(x->x==1?"":"$x",abs(r)).*["x";"y";"z"]'
#    sv=map(x->signbit(x)?"-":"+",v).*(eltype(v)<:Rational?map(x->"$(x.num)/$(x.den)",abs(v)):map(x->"$x",abs(v)))
#    sr[r.==0]=""; sv[v.==0]=""; so=similar(sv)
#    for i=1:r1; so[i]=lstrip(prod(sr[i,:])*sv[i],'+'); end
#    return join(so,",")
#end
Base.one(::String)="" # necessary for taking the product of a Matrix of Strings along one direction
Base.one(::Type{String})=""
function encodegenerator_new{T<:Real,R<:Real}(r::Matrix{R},v::Vector{T}=zeros(R,size(r,1)))
    (r1,r2)=size(r)
    @assert r1==r2==length(v)==3
    sr=map(x->signbit(x)?"-":"+",r).*map(x->x==1?"":"$x",abs(r)).*["x";"y";"z"]'
    sv=map(x->signbit(x)?"-":"+",v).*(eltype(v)<:Rational?map(x->"$(x.num)/$(x.den)",abs(v)):map(x->"$x",abs(v)))
    sr[r.==0]=""; sv[v.==0]="";
    return join(map(x->lstrip(x,'+'),prod(sr,2).*sv),",")
end

generateall(group::Union{AbstractString,Integer,Symbol};k...)=generateall(generators(group)...;k...)
function generateall{T<:Real,R<:Real}(inrot::Array{R,3},intrans::Array{T,2};tol::Real=1e-4)
    ng=size(inrot,3)
    tyr=eltype(inrot); tyt=eltype(intrans)
    nfold=determineNfold(inrot,intrans)
    expected=prod(nfold); found=1 # no rotation, no translation is first "found"
    outrot=zeros(tyr,3,3,expected)
    outtrans=zeros(tyt,3,expected)
    outrot[[1,5,9]]=one(tyr)
    ri=zeros(tyr,3,3); ti=zeros(tyt,3);
    for i=1:ng # for each input symmetry operation
        r0=inrot[:,:,i];
        t0=intrans[:,i]
        ri[[2,3,4,6,7,8]]=zero(tyr);
        ri[[1,5,9]]=one(tyr);
        ti[:]=zero(tyt) # reset ri and ti
        for j=1:(nfold[i]-1) # go through each of its n-fold equivalent versions
            ri=r0*ri    # calculate each version's rotation
            ti=r0*ti+t0 # and translation
            for k=1:found # go through each of the already-proven-unique symmetry operations
                # and add the trial symmetry operation to the already-proven-unique operation
                rs=ri*outrot[:,:,k];
                ts=mod(ri*outtrans[:,k]+ti,1)
                # if the combined operations are in the output lists, add the combined result
                if all((vec(sum(abs(outrot.-rs),[1,2])).>tol)|(vec(sum(abs(outtrans.-ts),1)).>tol))
                    found+=1; outrot[:,:,found]=rs; outtrans[:,found]=ts
                end
            end
        end
    end
    return(outrot[:,:,1:found],outtrans[:,1:found])
end

atompositions(group::Union{AbstractString,Integer,Symbol},o...;k...)=atompositions(generateall(generators(group)...;k...)...,o...;k...)
function atompositions{T<:Real,R<:Real,S<:AbstractFloat}(rot::Array{R,3},trans::Array{T,2},inpos::Array{S,2},lat::Lattice=one(Direct{S});tol::Real=1e-4)
    # all possible atom positions in first unit cell, {PositionNumber}x{AtomNumber} Lattice3Vectors
    apap = Lattice3Vector(mod(_multMv3(rot,inpos).+trans,1),lat)
    pos = Array{eltype(apap)}(0)
    typ = Array{UInt16}(0)
    flg = Array{Bool}(0)
    gen = Array{R}(3,3,0)
    # check for repeated coordinate sets for each atom
    for i=1:size(apap,2)
        dap=discernable(apap[:,i],tol) # discernable atom positions
        # check for near-matches across the unit cell boundary by mapping points within 1-tol < {h,k,l} < 1 to -tol < {h,k,l} < 0
        dap&=discernable(Lattice3Vector(mod(gethkl(apap[:,i])+tol,1)-tol,lat),tol)
        # pick out just the discernable atom postions and their types
        pos=vcat(pos, apap[dap,i])
        typ=vcat(typ,convert(eltype(typ),i)*ones(eltype(typ),sum(dap)))
        # it's important to keep track of which atom positions have been generated, and their rotations
        flg=vcat(flg,[false;trues(sum(dap)-1)])
        gen=cat(3,gen, rot[:,:,dap])
    end
    return (pos,typ,flg,gen)
end


function determineNfold{T<:Real,R<:Real}(r::Array{T,3},t::Array{R,2})
    @assert size(r,3)==size(t,2)
    N=zeros(T,size(r,3))
    for i=1:length(N); N[i]=determineNfold(r[:,:,i],t[:,i]); end
    return N
end
function determineNfold{T<:Real,R<:Real}(rotate::Matrix{T},translate::Vector{R}=T[0,0,0])
    N=one(T);
    tol = 1e-5;
    r = rotate;
    t = translate;
    while (norm(r-eye(T,3))>tol || norm(t)>tol) && N<10
        r=rotate*r
        t=mod(rotate*t+translate,1)
        N+=1
    end
    return N
end

#convertspecialcharacters(a::ASCIIString)=a
function convertspecialcharacters(a::String)
    subscripts=["₀","₁","₂","₃","₄","₅","₆","₇","₈","₉"]
    for i=1:length(subscripts)
        while contains(a,subscripts[i])
            idx=search(a,subscripts[i])
            mni=minimum(idx)
            mxi=maximum(idx)
            r1=1:(mni-1)
            r2=nextind(a,mxi):endof(a)
            a=a[r1]*"$(i-1)"*a[r2]
        end
    end
    bar="̄"
    while contains(a,bar)
        idx=search(a,bar)
        mni=minimum(idx)
        mxi=maximum(idx)
        r1=mni>1 ? collect(1:(mni-2)) : collect(1:(mni-1))
        r2=mni>1 ? [mni-1;collect(nextind(a,mxi):endof(a))]:collect(nextind(a,mxi):endof(a))
        a=a[r1]*"-"*a[r2]
    end
    return a #::ASCIIString # a should now be an ASCII string, or there are leftover symbols to remove
end
function standardize(g::AbstractString)
    sg=convertspecialcharacters(g)
    return replace(sg," ","")
end

generatorstring(g::AbstractString)=generatorstring(get(symbol2no,standardize(g),-1))
#generatorstring(g::Symbol)=generatorstring(get(symbol2no,g,-1))
function generatorstring(g::Integer)
    (g<1||g>length(no2gens)) && (g=1; warn("Unkown space group, returning generators for P 1."))
    no2gens[g]
end

function group_string_number(a::Integer)
    n= 0<a<=230 ? UInt8(a) : 0x01
    s=no2symbol[n]
    return (s,n::UInt8)
end
function group_string_number(a::AbstractString)
    n=UInt8(get(symbol2no,standardize(a),1))
    s=no2symbol[n]
    return (s,n)
end
include("Symmetry/symmetry.dict")
#"""
#The `const` dictionary `symbol2no` contains space group symbols as keys and their
#group numbers as values. The symbols have been stripped of spaces and UTF-8 characters.
#"""
#symbol2no
#"""
#The `const` `Vector{ASCIIString}` `no2gens` contains the generators of space group symmetry
#opperators where the ith string contains the generators for the ith group. The function
#`decodegenerator` converts the string form of a single generator to a rotation matrix and
#translation vector. The function `decodegenerators` does the same for a semi-colon separated list
#of string generators. The function `generators` takes either a space group name or number,
#pulls the encoded generators from `no2gens` and then uses `decodegenerators` to return the
#rotation matricies as a 3x3xN `Int8` array and the translation vectors as a 3xN `Rational{Int8}` array
#"""
#no2gens
