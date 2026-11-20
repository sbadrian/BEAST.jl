@testitem "Nédélec functions on triangles in 2D agree with their curl" begin
    using CompScienceMeshes
    using LinearAlgebra
    using StaticArrays

    # By Stokes' theorem, the integral of the curl of a function over a triangle is its
    # circulation counter-clockwise around the boundary. The tangential component of a
    # lowest-order Nédélec function is constant along each edge, so the circulation
    # along an edge is its value at the midpoint, dotted with the edge vector.
    corners = [(0.0, 0.0), (1.13, 0.27), (0.19, 0.82)]
    for order in ((1, 2, 3), (1, 3, 2))  # counter-clockwise and clockwise cells
        for udim in (2, 3)
            point(c) = SVector(ntuple(i -> i <= 2 ? c[i] : 0.0, udim))
            mesh = Mesh(point.(corners), [CompScienceMeshes.SimplexGraph(order...)])
            ch = chart(mesh, first(mesh))
            refs = BEAST.NDRefSpace{Float64}()
            # The parameters of the local vertices of the chart, and the local vertices
            # in counter-clockwise order (corners 1, 2, 3)
            uv(k) = SVector(k == 1 ? (1.0, 0.0) : k == 2 ? (0.0, 1.0) : (0.0, 0.0))
            local_ccw = [findfirst(==(i), order) for i in (1, 2, 3)]
            for r in 1:3
                curl = refs(neighborhood(ch, SVector(1 / 3, 1 / 3)))[r].curl
                circulation = sum(1:3) do e
                    ka, kb = local_ccw[e], local_ccw[mod1(e + 1, 3)]
                    value = refs(neighborhood(ch, (uv(ka) + uv(kb)) / 2))[r].value
                    edge = ch.vertices[kb] - ch.vertices[ka]
                    return dot(value[1:2], edge[1:2])
                end
                # The scalar curl in 2D is ẑ ⋅ ∇ × f; in 3D, n ⋅ ∇ × f with the normal of
                # the cell, which is -ẑ for the clockwise cell
                sgn = udim == 3 && order == (1, 3, 2) ? -1 : 1
                @test sgn * curl * volume(ch) ≈ circulation
            end
            # The functions in 2D equal those in 3D, in the plane z = 0
            if udim == 2
                mesh3 = Mesh(
                    [SVector(c[1], c[2], 0.0) for c in corners],
                    [CompScienceMeshes.SimplexGraph(order...)],
                )
                ch3 = chart(mesh3, first(mesh3))
                p = SVector(0.2, 0.3)
                v2 = refs(neighborhood(ch, p))
                v3 = refs(neighborhood(ch3, p))
                n = normal(ch3)
                for r in 1:3
                    # 3D functions are n × (...), 2D ones ẑ × (...)
                    @test v2[r].value ≈ dot(n, SVector(0, 0, 1)) * v3[r].value[1:2]
                end
            end
        end
    end
end
