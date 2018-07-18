# Instead of loading in all GeometryTypes junk, make our own simple Mesh type:
struct SampleMesh{T,N,M,O}
    vertices :: SVector{O,SVector{3,T}}
    indicies :: SMatrix{N,M,UInt}
    function SampleMesh{T,N,M,O}(verts::SVector{O,SVector{3,T}},inds::SMatrix{N,M,UInt}) where {O,T,N,M}
        @assert all(inds.<=O) "The indicies must index into the vertices vector"
        new{T,N,M,O}(verts,inds)
    end
end

function Base.show(io::IO,p::SampleMesh{T,N,M,O}) where {T,N,M,O}
    Base.print(io,"$O vertices with $N $M-vertex faces")
end
Base.show(io::IO,::MIME"text/plain",p::SampleMesh{T}) where T=Base.print(io,"SampleMesh{$T}:\n ",p)

function SampleMesh(verts::Vector{SVector{3,T}},inds::Matrix{I}) where {T,I<:Integer}
    O=length(verts)
    @assert all(0.<inds.<=O) "The indicies must index into the vertices vector"
    N,M=size(inds)
    SampleMesh{T,N,M,O}(SVector{O,SVector{3,T}}(verts),SMatrix{N,M,UInt}(inds))
end
function SampleMesh(verts::SVector{O,SVector{3,T}},inds::SMatrix{N,M,I}) where {O,T,N,M,I<:Integer}
    @assert all(0.<inds.<=O) "The indicies must index into the vertices vector"
    SampleMesh{T,N,M,O}(verts,convert.(UInt,inds))
end

# utility functions for Geometric calculations:
#const no_intersection= (false, 0., SVector{3}(0.,0.,0.))
function no_intersection(T)
    return (false,zero(T),SVector{3,T}(zero(T),zero(T),zero(T)))
end
function ray_face_intersect(ri::SVector{3,T},rf::SVector{3,T},a::SVector{3,R},b::SVector{3,R},c::SVector{3,R}) where {T,R}
    S=promote_type(T,R)
    # ri,rf -- coordinates for the start and end point of a ray
    # a,b,c -- coordinates for the three corners of the facet
    # n     -- normal direction for the facet (if a,b,c defined in the right order this is (b-a)×(c-a) )
    v1=b-a
    v2=c-a
    n=cross(v1,v2) # calculate the normal instead of relying on the mesh-stored normal
    n̂=n/norm(n)
    r=rf-ri # the ray vector
    x=( dot(a,n̂) - dot(ri,n̂) )/dot(r,n̂)
    0<=x<=1 || return no_intersection(S)
    ip=ri+x*r # vector to the intersection point (from the origin)
    p=ip-a-dot(ip-a,n̂)*n̂ # in-facet-plane vector pointing from a to the intersection

    v̂1=v1/norm(v1)
    v̂2=v2/norm(v2)
    v̂1v̂2=dot(v̂1,v̂2)
    den=1-v̂1v̂2^2
    pv̂1=dot(p,v̂1)
    pv̂2=dot(p,v̂2)
    #sv1= ( pv̂1 - v̂1v̂2*pv̂2 )/den * v̂1
    #tv2= ( pv̂2 - v̂1v̂2*pv̂1 )/den * v̂2
    #s=norm(sv1)/norm(v1)
    #t=norm(tv2)/norm(v2)
    s= ( pv̂1 - v̂1v̂2*pv̂2 )/den/norm(v1)
    t= ( pv̂2 - v̂1v̂2*pv̂1 )/den/norm(v2)
    0<s<1 && 0<t<1 && s+t<1 || return no_intersection(S)
    #return (true, x*norm(r)*sign(dot(r,n)), x*r+ri, a+sv1+tv2) # (flag, ray-length before hitting facet + out/in == (+,-), intersection point calculated two ways)
    return (true, x*norm(r)*sign(dot(r,n)), ip)
end
function ray_mesh_intersect(ri::SVector{3,T},rf::SVector{3,T},mesh::SampleMesh{R,N}) where {T,R,N}
    # a SampelMesh has fields :vertices, :indicies
    all_intersections = Array{typeof(no_intersection(promote_type(T,R)))}(N)
    for i=1:N
        all_intersections[i]=ray_face_intersect(ri,rf,mesh.vertices[mesh.indicies[i,:]]...)
    end
    return all_intersections
end
function ray_mesh_pathlength(ri::SVector{3,T},rf::SVector{3,T},mesh::SampleMesh{R}) where {T,R}
    all_intersections=ray_mesh_intersect(ri,rf,mesh)
    does_intersect = [x[1] for x in all_intersections]
    any(does_intersect) ? sum(x[2] for x in all_intersections[does_intersect]) : zero(promote_type(T,R))
end

function transform(r::SMatrix{3,3,R},mesh::SampleMesh{T,N,M,O}) where {R,T,N,M,O}
    P=promote_type(R,T)
    tverts = Array{SVector{3,P}}(O) # O vertices in the mesh
    for i=1:O; tverts[i]=r*mesh.vertices[i]; end
    SampleMesh(SVector{O,SVector{3,P}}(tverts),mesh.indicies)
end
transform(r::Matrix,f)=transform(SMatrix{size(r)...}(r),f)

function axis_limits(m::SampleMesh,v::SVector{3})
    extrema(dot.([v/norm(v)],m.vertices))
end
function axis_limits(m::SampleMesh{T}) where T
    v=m.vertices
    [ extrema(dot.([q],v)) for q in ( SVector{3,T}(1,0,0), SVector{3,T}(0,1,0), SVector{3,T}(0,0,1) ) ]
end
function line_integration_limits(ri::SVector{3,T},rf::SVector{3,T},mesh::SampleMesh{R}) where {T,R}
    all_intersections=ray_mesh_intersect(ri,rf,mesh)
    does_intersect = [x[1] for x in all_intersections]
    any(does_intersect) ? maximum(x[2] for x in all_intersections[does_intersect]) : zero(promote_type(T,R))
end

function determine_second_moments(mesh::SampleMesh{T};method::Symbol=:gausslegendre,points::Integer=5) where T
    o=SVector{3,T}(0,0,0)
    x=SVector{3,T}(1,0,0)
    y=SVector{3,T}(0,1,0)
    z=SVector{3,T}(0,0,1)
    (p,w)=FastGaussQuadrature.gausslegendre(points)
    result=zeros(T,3,3)
    volume=zero(T)
    mesh_limits=axis_limits(mesh)
    xmin,xmax=mesh_limits[1]
    ylim=mesh_limits[2]; zlim=mesh_limits[3]
    for i=1:points
        xi=(p[i]*(xmax-xmin)+(xmax+xmin))*x/2
        ymin=line_integration_limits(o+xi,o+xi+2*ylim[1]*y,mesh)*sign(ylim[1]) # since sign(ylim) determines the search direction
        ymax=line_integration_limits(o+xi,o+xi+2*ylim[2]*y,mesh)*sign(ylim[2])
        for j=1:points
            yj= (p[j]*(ymax-ymin) + (ymax+ymin))*y/2
            zmin=line_integration_limits(o+xi+yj,o+xi+yj+2*zlim[1]*z,mesh)*sign(zlim[1])
            zmax=line_integration_limits(o+xi+yj,o+xi+yj+2*zlim[2]*z,mesh)*sign(zlim[2])
            for k=1:points
                zk=(p[k]*(zmax-zmin)+(zmax+zmin))*z/2
                pijk=o+xi+yj+zk
                vijk = w[i]*(xmax-xmin)/2 * w[j]*(ymax-ymin)/2 * w[k]*(zmax-zmin)/2
                result+= vijk*(pijk*pijk.')
                volume+= vijk
            end
        end
    end
    sm=result/volume
    # set absolute values smaller than eps(T) to +/-0 (and convert to SecondMoments)
    convert(Matrix{SecondMoment{T}},sm.*(abs.(sm).>eps(T)))
end


pperr(x) = @sprintf("%5.3f",x)

# This absorption correction fuction requires that the facets defining the crystal
# have already been rotated to the correct orientation for the required Q point
# The incident beam is along the x-axis, and the outgoing beam is scattered by
# the angle tθ.
function absorption_correction(tθ::Real,μ::Real,mesh::SampleMesh{T};method::Symbol=:gausslegendre,points::Integer=8) where T
    #DEBUG && status(:debug,"Absorption correction using $method integration")
    o=SVector{3,T}(0,0,0)
    x=SVector{3,T}(1,0,0)
    y=SVector{3,T}(0,1,0)
    z=SVector{3,T}(0,0,1)
    result=zero(T)
    volume=zero(T)
    mesh_limits=axis_limits(mesh)
    biggest_dimension = maximum(maximum.([abs.(x) for x in mesh_limits]))
    din= 2*biggest_dimension*SVector{3,T}(-1,0,0)
    dout=2*biggest_dimension*SVector{3,T}(cos(tθ),sin(tθ),0)
    if method==:gausslegendre
        (p,w)=FastGaussQuadrature.gausslegendre(points)
        xmin,xmax=mesh_limits[1]
        ylim=mesh_limits[2]; zlim=mesh_limits[3]
        for i=1:points
            xi=(p[i]*(xmax-xmin)+(xmax+xmin))*x/2
            ymin=line_integration_limits(o+xi,o+xi+2*ylim[1]*y,mesh)*sign(ylim[1]) # since sign(ylim) determines the search direction
            ymax=line_integration_limits(o+xi,o+xi+2*ylim[2]*y,mesh)*sign(ylim[2])
            for j=1:points
                yj= (p[j]*(ymax-ymin) + (ymax+ymin))*y/2
                zmin=line_integration_limits(o+xi+yj,o+xi+yj+2*zlim[1]*z,mesh)*sign(zlim[1])
                zmax=line_integration_limits(o+xi+yj,o+xi+yj+2*zlim[2]*z,mesh)*sign(zlim[2])
                for k=1:points
                    zk=(p[k]*(zmax-zmin)+(zmax+zmin))*z/2
                    v=o+xi+yj+zk
                    Ti=ray_mesh_pathlength(v,v+din,mesh)
                    To=ray_mesh_pathlength(v,v+dout,mesh)
                    if Ti>=0 && To>=0 && Ti+To>0 # otherwise this point isn't in the crystal
                        #DEBUG && status(:debug,"The path length is $(pperr(Ti+To)) giving exponential factor $(pperr(exp(-μ*(Ti+To))))")
                        vijk = w[i]*(xmax-xmin)/2 * w[j]*(ymax-ymin)/2 * w[k]*(zmax-zmin)/2
                        result += vijk*exp(-μ*(Ti+To))
                        volume += vijk # the volume is the integral of f(x)=1
                    end
                end
            end
        end
    elseif method==:montecarlo
        volume = points^3 # the monte-carlo method needs to be normalized by the number of trials
        meshmins   = minimum.(mesh_limits)
        meshranges = maximum.(mesh_limits)-meshmins
        eamat = rand(3,volume).*meshranges .+ meshmins
        evalat = [SVector{3,T}(eamat[:,i]) for i=1:volume]
        for i=1:volume
            v=o+evalat[i]
            Ti=ray_mesh_pathlength(v,v+din,mesh)
            To=ray_mesh_pathlength(v,v+dout,mesh)
            if Ti>=0 && To>=0 && Ti+To>0 # otherwise this point isn't in the crystal
                #DEBUG && status(:debug,"The path length is $(pperr(Ti+To)) giving exponential factor $(pperr(exp(-μ*(Ti+To))))")
                result += exp(-μ*(Ti+To))
                #volume += 1 # the volume is the integral of f(x)=1
            end
        end
    end

    #return (result,volume)
    invA= result>0 ? volume/result : 0 # A = (∭e^{-μT}dV)/V. if V≡0 A=∞, which we might want to protect againts.
    return invA # since we care about 1/A, V≡0 is OK, but we need to protect against result≡0
end

function absorption_correction(tθ::Real,μ::Measured,mesh::SampleMesh{T};method::Symbol=:gausslegendre,points::Integer=8) where T
    o=SVector{3,T}(0,0,0)
    x=SVector{3,T}(1,0,0)
    y=SVector{3,T}(0,1,0)
    z=SVector{3,T}(0,0,1)
    (p,w)=FastGaussQuadrature.gausslegendre(points)
    din=SVector{3,T}(-1,0,0)
    dout=SVector{3}(cos(tθ),sin(tθ),0)
    result=zero(T)
    volume=zero(T)
    mesh_limits=axis_limits(mesh)
    xmin,xmax=mesh_limits[1]
    ylim=mesh_limits[2]; zlim=mesh_limits[3]
    # uncertainty stuff:
    μval = value(μ)
    μvar = variance(μ)
    VδAδμ = zero(T)
    for i=1:points
        xi=(p[i]*(xmax-xmin)+(xmax+xmin))*x/2
        ymin=line_integration_limits(o+xi,o+xi+2*ylim[1]*y,mesh)*sign(ylim[1]) # since sign(ylim) determines the search direction
        ymax=line_integration_limits(o+xi,o+xi+2*ylim[2]*y,mesh)*sign(ylim[2])
        for j=1:points
            yj= (p[j]*(ymax-ymin) + (ymax+ymin))*y/2
            zmin=line_integration_limits(o+xi+yj,o+xi+yj+2*zlim[1]*z,mesh)*sign(zlim[1])
            zmax=line_integration_limits(o+xi+yj,o+xi+yj+2*zlim[2]*z,mesh)*sign(zlim[2])
            for k=1:points
                zk=(p[k]*(zmax-zmin)+(zmax+zmin))*z/2
                v=o+xi+yj+zk
                Ti=ray_mesh_pathlength(v,v+din,mesh)
                To=ray_mesh_pathlength(v,v+dout,mesh)
                if Ti>=0 && To>=0 && Ti+To>0 # otherwise this point isn't in the crystal
                    vijk = w[i]*(xmax-xmin)/2 * w[j]*(ymax-ymin)/2 * w[k]*(zmax-zmin)/2
                    result += vijk*exp(-μval*(Ti+To))
                    VδAδμ  += vijk*(Ti+To)exp(-μval*(Ti+To))
                    volume += vijk # the volume is the integral of f(x)=1
                end
            end
        end
    end
    #return (result,volume)
    invA= result>0 ? volume/result : 0 # A = (∭e^{-μT}dV)/V. if V≡0 A=∞, which we might want to protect againts.
    δinvAδμ = (VδAδμ/volume) * invA^2 # chain rule
    return Measured(invA,δinvAδμ^2*μvar) # + δinvAδT^2*Tvar + ...
end


function make_circle_z(d::T,N::Integer=25;center::SVector{3,T}=SVector{3,T}(0,0,0)) where T
    @assert N>2
    θ=linspace(0,2π,N+1)
    r=d/2
    ring=[SVector{3,T}(r*cos(x),r*sin(x),0) for x in θ[1:end-1]] # 1:end-1 to avoid the doubled point
    verts=vcat([center],[center].+ring)
    idxs=1+hcat( zeros(Int,N), 1:N, circshift(1:N,-1) ) # 1-based indexing
    return (verts,idxs)
end

function make_cylinder_z(l::T,d::R,N::Integer=25;center::SVector{3,T}=SVector{3,T}(0,0,0)) where {T,R}
    P=promote_type(T,R)
    tvrt,tidx=make_circle_z(convert(P,d),N;center=SVector{3,P}(0,0, l/2))
    bvrt,bidx=make_circle_z(convert(P,d),N;center=SVector{3,P}(0,0,-l/2))
    verts=vcat(tvrt,bvrt).+[center] # 2(N+1) vertices spanning 0:2N+1
    bidx = bidx[:,[1,3,2]] + N+1 # the bottom face indicies must be offset and flipped
    #bnrm *= -1 # and their normals must be reversed
    i1 = 1:N; i2 = N+2:2N+1; i3=circshift(i2,-1); i4=circshift(i1,-1);
    midx = 1+reshape(hcat( i1,i2,i3,  i1,i3,i4 ).',(3,2N)).' # 1-based indexing
    return (verts, vcat(tidx,midx,bidx))
end

function make_cylinder_z_better(l::T,d::R,N::Integer=25;center::SVector{3,T}=SVector{3,T}(0,0,0)) where {T,R}
    P=promote_type(T,R)
    tvrt,tidx=make_circle_z(convert(P,d),N;center=SVector{3,P}(0,0, l/2))
    rot=SMatrix{3,3,P}([cos(pi/N) -sin(pi/N) 0; sin(pi/N) cos(pi/N) 0; 0 0 1])
    tvrt=[rot].*tvrt
    bvrt,bidx=make_circle_z(convert(P,d),N;center=SVector{3,P}(0,0,-l/2))
    verts=vcat(tvrt,bvrt).+[center] # 2(N+1) vertices spanning 0:2N+1
    bidx = bidx[:,[1,3,2]] + N+1 # the bottom face indicies must be offset and flipped
    #bnrm *= -1 # and their normals must be reversed
    i1 = 1:N; i2 = N+2:2N+1; i3=circshift(i2,-1); i4=circshift(i1,-1);
    midx = 1+reshape(hcat( i1,i2,i3,  i1,i3,i4 ).',(3,2N)).' # 1-based indexing
    return (verts, vcat(tidx,midx,bidx))
end

mesh_cylinder_z(o...;k...)=SampleMesh(make_cylinder_z(o...;k...)...)

#
# # Sometimes it's useful to visualize what your mesh looks like:
# function plotvecs3D(vs::Vector{T},o...;k...) where T<:StaticArrays.SVector
#     PyPlot.plot3D([v[1] for v in vs],[v[2] for v in vs],[v[3] for v in vs],o...;k...)
# end
# function plotwireframe(vs::AbstractVector{T},vi::AbstractArray{I,2},o...;k...) where {T<:StaticArrays.SVector,I<:Integer}
#     allv=vs[vi]
#     plotvecs3D(vec(hcat( allv, NaN*allv[:,end] ).'),o...;k...)
# end
# function plotmesh(a::SampleMesh,o...;k...)
#     plotwireframe(a.vertices,a.indicies,o...;k...)
# end
