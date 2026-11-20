"""
    NDLCCBasis

Lowest-order Nédélec (edge) functions on a tetrahedral mesh, see `nedelecc3d`.
"""
mutable struct NDLCCBasis{T,M,P} <: Space{T}
    geo::M
    fns::Vector{Vector{Shape{T}}}
    pos::Vector{P}
end

NDLCCBasis(geo, fns) = NDLCCBasis(geo, fns, Vector{vertextype(geo)}(undef,length(fns)))

refspace(space::NDLCCBasis{T}) where {T} = NDLCCRefSpace{T}()
Base.similar(space::NDLCCBasis{T,M,P} where {T,M,P}, geo, fns, pos) = NDLCCBasis(geo, fns, pos)

"""
    nedelecc3d(mesh, edges=skeleton(mesh, 1))

The lowest-order Nédélec functions of the first kind (curl-conforming edge functions) on
a tetrahedral `mesh`, as an `NDLCCBasis`: one function per cell of the mesh of
segments `edges`, which must share the vertices of `mesh`. Each function is supported
on the tetrahedra that share its edge, has a continuous tangential component, and has
circulation ∫ f ⋅ t ds = +1 along its edge, oriented from the first to the second
vertex of the cell in `edges`, and 0 along all other edges. Leaving out edges (for
example those on a perfectly conducting boundary) sets the tangential component to zero
there.

The integrals over the tetrahedra use the signed volume of their charts, so the
tetrahedra of `mesh` should be positively oriented; the mass matrix of a mesh with
negatively oriented tetrahedra is negative definite.

`curl` of the basis is an `NDLCDBasis`, and `ttrace` gives the
tangential trace on the boundary as an `RTBasis`.
"""
function nedelecc3d(mesh, edges)
    T = coordtype(mesh)
    P = vertextype(mesh)
    num_edges = numcells(edges)

    C = connectivity(edges, mesh, identity)
    rows = rowvals(C)
    vals = nonzeros(C)

    fns = Vector{Vector{Shape{T}}}(undef,num_edges)
    pos = Vector{P}(undef,num_edges)
    for (i,e) in enumerate(edges)

        fns[i] = Vector{Shape{T}}()
        pos[i] = cartesian(center(chart(edges,e)))

        for k in nzrange(C,i)

            j = rows[k]
            s = vals[k]
            push!(fns[i], Shape{T}(j, abs(s), sign(s)))
        end
    end

    NDLCCBasis(mesh, fns, pos)
end

function nedelecc3d(mesh)
    edges = skeleton(mesh,1)
    nedelecc3d(mesh, edges)
end


"""
    nedelec(mesh::AbstractMesh{U,4}, edges=skeleton(mesh, 1))

The edge functions of a tetrahedral mesh, see `nedelecc3d`.
"""
function nedelec(mesh::CompScienceMeshes.AbstractMesh{U,4} where {U},
    edges=skeleton(mesh,1))

    nedelecc3d(mesh, edges)
end

curl(space::NDLCCBasis, geo, fns) = NDLCDBasis(geo, fns, space.pos)

ttrace(X::NDLCCBasis, geo, fns) = RTBasis(geo, fns, deepcopy(X.pos))
