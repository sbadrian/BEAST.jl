"""
    NDBasis

Lowest-order Nédélec (edge) functions on a triangle mesh, see `nedelec`.
"""
struct NDBasis{T,M,P} <: Space{T}
    geo::M
    fns::Vector{Vector{Shape{T}}}
    pos::Vector{P}
end

NDBasis(geo, fns) = NDBasis(geo, fns, Vector{vertextype(geo)}(undef,length(fns)))

refspace(s::NDBasis) = NDRefSpace{scalartype(s)}()


"""
    nedelec(surface, edges=skeleton(surface, 1))

The lowest-order Nédélec functions of the first kind (curl-conforming edge functions) on
a triangle mesh `surface`, one per cell of the mesh of segments `edges`, which must
share the vertices of `surface`. Each function is supported on the triangles adjacent to
its edge, has a continuous tangential component across the edge, and has circulation
∫ f ⋅ t ds = +1 along its edge, oriented from the first to the second vertex of the
cell in `edges`, and 0 along all other edges.

On triangles in 3D, the functions are tangential to the triangles, and their `curl` is
the normal component n ⋅ ∇ × f, with the normal of each triangle. On triangles with 2D
coordinates, the functions are those of the same triangles in the plane z = 0, oriented
with normal ẑ, and the curl is ẑ ⋅ ∇ × f.

For a tetrahedral mesh, `nedelec(mesh, edges)` returns the volume functions of
`nedelecc3d`.
"""
function nedelec(surface, edges=skeleton(surface,1))

    T = coordtype(surface)
    # P = eltype(surface.vertices)
    P = vertextype(surface)
    num_edges = numcells(edges)

    C = connectivity(edges, surface, identity)
    rows = rowvals(C)
    vals = nonzeros(C)

    fns = Vector{Vector{Shape{T}}}(undef,num_edges)
    pos = Vector{P}(undef,num_edges)
    for (i,edge) in enumerate(edges)

        fns[i] = Vector{Shape{T}}()
        pos[i] = cartesian(center(chart(edges,edge)))

        for k in nzrange(C,i)

            j = rows[k] # j is the index of a cell adjacent to edge
            s = vals[k] # s contains the oriented (signed) local index of edge[i] in cell[j]

            # i == 3 && @show s
            push!(fns[i], Shape{T}(j, abs(s), sign(s)))
        end
    end

    NDBasis(surface, fns, pos)
end

# function LinearAlgebra.cross(::NormalVector, s::NDBasis)
#     # @assert CompScienceMeshes.isoriented(s.geo)
#     RTBasis(s.geo, s.fns, s.pos)
# end

function LinearAlgebra.cross(::NormalVector, s::NDBasis)
    @assert CompScienceMeshes.isoriented(s.geo)
    fns = similar(s.fns)
    for (i,fn) in pairs(s.fns)
        fns[i] = [Shape(sh.cellid, sh.refid, -sh.coeff) for sh in fn]
    end
    RTBasis(s.geo, fns, s.pos)
end
function LinearAlgebra.cross(s::NDBasis, ::NormalVector)
    @assert CompScienceMeshes.isoriented(s.geo)
    fns = similar(s.fns)
    for (i,fn) in pairs(s.fns)
        fns[i] = [Shape(sh.cellid, sh.refid, sh.coeff) for sh in fn]
    end
    RTBasis(s.geo, fns, s.pos)
end
function curl(space::NDBasis)
    divergence(n × (n × (n × space)))
end
