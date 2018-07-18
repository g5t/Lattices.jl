
function sortcouplings(crystal::Crystal,tol=1e-5)
    lm=find(norm(crystal.spin) .> 0) #vector of moment-bearing atom indicies
    nm=length(lm)
    (p1,p2,p3)=halfcubepoints(crystal.jcel)
    np=length(p1)
    npm=np*nm
    atrn=Array{eltype(crystal.pos)}(npm) # unit cell translation
    apos=Array{eltype(crystal.pos)}(npm) # atom position
    aidx=Array{eltype(crystal.jion)}(npm)# atom index
    # calculate all magnetic ion positions as Lattice3Vectors
    for i=1:np
        for j=1:nm
            k=(i-1)*nm+j
            atrn[k]=Lattice3Vector(p1[i],p2[i],p3[i],crystal.dlat)
            apos[k]=crystal.pos[lm[j]] + atrn[k]
            aidx[k]=lm[j]
        end
    end
    # number of couplings to keep track of
    nc=fld((2np+1)*nm^2-nm,2)# 2np+1 is odd, if nm is even(odd) nm^2 is e(o), so (2np+1)*nm^2 is e(o) and (2np+1)*nm^2-nm is e(e)
    cdst=Array{eltype(norm(crystal.pos[1]))}(nc)
    ctrn=Array{eltype(atrn)}(nc)
    cidx=Array{eltype(aidx)}(2,nc)
    for i=1:nm
        for j=1:npm
            k=(i-1)*npm+j
            cdst[k]=norm(crystal.pos[lm[i]]-apos[j])
            ctrn[k]=atrn[j]
            cidx[1,k]=lm[i]
            cidx[2,k]=aidx[j]
        end
    end
    k=npm*nm
    # add on the upper/lower elements ( there are (N^2-N)/2 upper triangle elements in an NxN matrix)
    # these were left-out by halfcubepoints because it does not return (0,0,0)
    for i=1:(nm-1)
        for j=(i+1):nm
            k+=1
            cdst[k]=norm(crystal.pos[lm[i]]-crystal.pos[lm[j]])
            ctrn[k]=0*ctrn[1]
            cidx[1,k]=lm[i]
            cidx[2,k]=lm[j]
        end
    end
    @assert nc==k "Running index should be $nc but is $k instead."
    # verify all distances above threshhold
    @assert all(cdst.>=crystal.mindist)
    # Remove too-large distances
    bad= cdst .> crystal.maxdist
    cdst=cdst[!bad]
    ctrn=ctrn[!bad]
    cidx=cidx[:,!bad]
    # Find the sorting permutation for the remaining distances
    p=sortperm(cdst)
    cdst=cdst[p]
    ctrn=ctrn[p]
    cidx=cidx[:,p]
    # Label equivalent distances by increasing integers
    edst=cumsum([true;diff(cdst).>tol]) # only assigns a new integer for ith distance if distance to (i-1)th [smaller] distance is larger than `tol`

    if true #false
        # Equal-distance couplings are not necessarily symmetry-equivalent couplings.
        # For each set of equal-distance couplings determine which are symmetry-equivalent.
        # Note that this process should not produce new couplings within the search-volume
        # but it is possible that large-distance couplings could generate new out-of-search-volume
        # couplings.
        (rot,trn)=generateall(crystal.group)
        ndst=Array{eltype(cdst)}(0)
        ntrn=Array{eltype(ctrn)}(0)
        nidx=Array{eltype(cidx)}(2,0)
        neqd=Array{eltype(cidx)}(0)
        j=1
        for i=1:maximum(edst)
            thisi= find(edst .== i)
            while length(thisi)>0 #&& j<10^3
                (genctrn,gencidx,gencdst)=generatecouplings(crystal.pos,ctrn[thisi[1]],cidx[:,thisi[1]],rot,trn;tol=tol)
                # find the unique newly-generated couplings by comparing Gi[v,a1,a2]==Gj[v,a1,a2] and Gi[-v,a2,a1]==Gj[v,a1,a2]
                unc=uniquecoupling(genctrn,gencidx)
                # the couplings generated by ctrn[this[1]]... are likely symmetry-equivalent to some of
                # the couplings described by ctrn[this[2:end]]...
                # find the logical vector of symmetry-equivalent couplings to thisi[1] in thisi
                # by comparing ctrn[thisi] and cidx[:,thisi] to symmetry-generated couplings
                # `se` is a logical vector the same length as ctrn[thisi] with value `true` if that
                # equal-distance coupling is symmetry-equivalent with ctrn[thisi[1]]
                se=symmetryequivalent(genctrn,gencidx,ctrn[thisi],cidx[:,thisi],tol)
                #remove the list of already generated coupling indicies from cdst, ctrn, cidx
                nse=sum(se) # number of already-generated couplings
                nuq=sum(unc) # number of unique newly-generated couplings
                npg=length(se) # number of pre-generated couplings
                ngn=length(genctrn) # number of generated couplings
                if nse<nuq
                    # It's entirely possible that a symmetry-opperation could generate
                    # a coupling that was missed by the simple half-cube equal-distance search.
                    # If this is the case then accept the addition, otherwise warn the user.
                    #
                    # `extrainsearch` returns a logical vector indicating which extra generated
                    # elements of genctrn were in the brute-force search
                    wis=extrainsearch(hcat(p1,p2,p3)',ctrn[thisi[se]],genctrn,tol)
                    if sum(wis)>0
                        info("Symmetry generation of couplings seems to have produced extra couplings that the brute-force search missed! You should investigate")
                        println(" ",sum(wis)," extra in-search couplings ",genctrn[wis])
                    end
                elseif nse>nuq
                    error("Symmetry error with i=$i, j=$j, d=$(cdst[thisi[1]]); $nse of $npg equal-distance couplings are symmetry equivalent but $nuq of $ngn generated couplings are unique")
                end
                # remove the symmetry-equivalent couplings from the list of equal-distance couplings
                thisi=thisi[!se]
                # store the unique-generated symmetry-equivalent coupling information in n*
                ndst=vcat(ndst,gencdst[unc])
                ntrn=vcat(ntrn,genctrn[unc])
                nidx=hcat(nidx,gencidx[:,unc])
                neqd=vcat(neqd,ones(eltype(neqd),nuq)*j)
                j+=1
            end
        end
    else
        ndst=cdst
        ntrn=ctrn
        nidx=cidx
        neqd=edst
    end

    p=crystal.maxbonds>0 ? neqd.<=crystal.maxbonds : trues(size(neqd)...)
    return (ntrn[p],nidx[:,p],neqd[p],ndst[p])
end
sortcouplings!(crystal::Crystal,tol=1e-5)=((crystal.jvec,crystal.jion,crystal.jexn)=sortcouplings(crystal,tol);)

function extrainsearch{R<:Real,L<:Lattice3Vector}(p::Array{R,2},ed::Array{L,1},sg::Array{L,1},tol=1e-5)
    # `ed` -- equal-distance coupling unit-cell translation vectors
    # `sg` -- symmetry-generated coupling unit-cell translation vectors
    # each of `ed` should be in `sg`
    edinsg=falses(size(sg))
    for i=1:length(sg)
        any(norm(sg[i]-ed).<tol) && (edinsg[i]=true)
    end
    extras=find(!edinsg) # indicies into sg that are not represented in eg
    insearch=falses(size(sg))
    for i=1:length(extras)
        any(sumabs2(p.-gethkl(sg[extras[i]]),1).<tol) && (insearch[extras[i]]=true)
    end
    return insearch
end


function generatecouplings{L<:Lattice3Vector,I<:Integer,S<:Real,R<:Real}(pos::Array{L,1},dl::L,idx::Array{I,1},rot::Array{R,3},trans::Array{S,2};tol=1e-5)
    # pull-out position vectors for the two atoms, and create the unit-cell lattice vector
    r1=pos[idx[1]]
    r2=pos[idx[2]]
    # calculate the symmetry equivalent indicies of the two atoms and the lattice vector
    r1i=_multMv(rot,gethkl(r1)).+trans
    r2i=_multMv(rot,gethkl(r2)).+trans
    dli=_multMv(rot,gethkl(dl))-floor(r1i+tol)+floor(r2i+tol)
    # create vectors from the indicies for the symmetry equivalent positions
    r1n=Lattice3Vector(mod(r1i,1),r1.lat)
    r2n=Lattice3Vector(mod(r2i,1),r2.lat)
    dln=Lattice3Vector(dli,dl.lat)
    # none of the generated positions should be new
    # all should be represented in the first unit cell already
    # use findfirstmod1 to determine the first-unit-cell equivalent indicies of the generated positions
    atom1=findfirstmod1(pos,r1n,tol)
    atom2=findfirstmod1(pos,r2n,tol)
    @assert !(any(atom1.==0)|any(atom2.==0)) "Symmetry operations have produced new positions!"
    dist=norm(pos[atom2]-pos[atom1]+dln)
    ok=abs(dist-dist[1]).<tol
    any(!ok) && (info("Some symmetry generated couplings had wrong distances and have been dropped."))
    return (dln[ok],vcat(atom1[ok]',atom2[ok]'),dist[ok])
end
#findfirst{L<:Lattice3Vector}(a::Array{L,1},b::Array{L,1},tol=1e-5)=map(x->Base.findfirst(norm(a-x).<tol),b)
findfirstmod1{L<:Lattice3Vector}(a::Array{L,1},b::Array{L,1},tol=1e-5)=map(x->Base.findfirst(norm(Lattice3Vector(mod(gethkl(a-x)+tol,1)-tol,x.lat)).<tol),b)

#uniquecoupling{L<:LatticeVector,R<:Real}(v::Array{L,1},a::Array{R,2})=uniquecoupling(v,vec(a[1,:]),vec(a[2,:]))
function uniquecoupling{L<:LatticeVector,R<:Real}(v::Array{L,1},a::Array{R,2},tol=1e-5)
    nc=size(a,2)
    uc=trues(nc)
    for i=1:(nc-1)
        for j=(i+1):nc
            de=(norm(v[i]-v[j])<tol)&(a[1,i]==a[1,j])&(a[2,i]==a[2,j])# both couplings are directly equivalent
            re=(norm(v[i]+v[j])<tol)&(a[1,i]==a[2,j])&(a[2,i]==a[1,j])# both couplings have reverse equivalency
            (de|re) && (uc[j]=false) # we want to keep one unique copy of each coupling. old version got rid of all copies #(uc[i]=false;uc[j]=false)
        end
    end
    return uc
end
function symmetryequivalent{L<:LatticeVector,R<:Real,S<:Real}(gv::Array{L,1},ga::Array{R,2},cv::Array{L,1},ca::Array{S,2},tol=1e-5)
    ag=falses(size(cv)...)
    ng=length(gv)
    nc=length(cv)
    @assert size(ga,1)==size(ca,1)==2
    @assert size(ga,2)==ng
    @assert size(ca,2)==nc
    for i=1:nc
        j=0
        while j<ng && !ag[i]
            j+=1 # check for G[v,a1,a2]==C[v,a1,a2] and G[-v,a2,a1]==C[v,a1,a2]
            ((norm(gv[j]-cv[i])<tol)&(ga[1,j]==ca[1,i])&(ga[2,j]==ca[2,i])) && (ag[i]=true)
            ((norm(gv[j]+cv[i])<tol)&(ga[1,j]==ca[2,i])&(ga[2,j]==ca[1,i])) && (ag[i]=true)
        end
    end
    return ag
end
