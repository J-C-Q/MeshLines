@MakieCore.recipe(MeshLines, points) do scene
    MakieCore.Attributes(
        around=10,
        between = 0,
        width = 0.1
    )
end

# https://math.stackexchange.com/questions/133177/finding-a-unit-vector-perpendicular-to-another-vector
function normal(x::MakieCore.Vec3f)
    y = zeros(eltype(x),3)
    m = findfirst(m->m!=0,x)
    n = mod1(m+1,3)
    y[n] = x[m]
    y[m] = -x[n]
    # norm = sqrt(y[n]^2+y[m]^2)
    # y ./= norm
    return MakieCore.Vec3f(normalize(y))
end

function intersectPlane(x::MakieCore.Vec3f,y::MakieCore.Vec3f)
    z = zeros(eltype(x),3)
    if iszero(x)
        z .= normalize(y)
    elseif iszero(y)
        z .= normalize(x)
    else
        z .= normalize(x) .+ normalize(y)
    end

    # norm = sqrt(z[1]^2+z[2]^2+z[3]^2)
    # z ./= norm
    zvec = MakieCore.Vec3f(normalize(z))


    return zvec,normal(zvec)
end

function intersectPlaneNotZeroAndLinear(x::MakieCore.Vec3f,y::MakieCore.Vec3f)
    z = zeros(eltype(x),3)
    normal = zeros(eltype(x),3)
    z .= normalize(x) .+ normalize(y)
    normal = -1.0.*normalize(x) .+ normalize(y)


    # norm = sqrt(z[1]^2+z[2]^2+z[3]^2)
    # z ./= norm
    zvec = MakieCore.Vec3f(normalize(z))
    normalvec = MakieCore.Vec3f(normalize(normal))


    return zvec,normalvec
end

function cos_angle_between_vectors(a, b)
    return abs(dot(a, b)) / (norm(a) * norm(b))
end

function r_ellipse(ϕ;a=1,b=1)
    return a*b/sqrt(a^2*sin(ϕ)^2+b^2*cos(ϕ)^2)
end

function align_ring(prev_ring::Vector{MakieCore.Point3f}, cur_ring::Vector{MakieCore.Point3f})
    n = length(prev_ring)
    best_shift = 0
    best_sum = Inf
    # Try every possible circular shift
    for shift in 0:(n-1)
        # Use circular shift to reorder cur_ring
        shifted_ring = circshift(cur_ring, shift)
        total_dist = sum(norm(prev_ring[i] - shifted_ring[i]) for i in 1:n)
        if total_dist < best_sum
            best_sum = total_dist
            best_shift = shift
        end
    end
    return circshift(cur_ring, best_shift)
end

function MakieCore.plot!(ml::MeshLines{<:Tuple{<:AbstractVector{<:MakieCore.Point3f}}})

    points = ml[1]
    around = ml.around

    vertices_ob = MakieCore.Observable(Matrix{Float32}(undef, length(points[]) * around[] + 2, 3))
    faces_ob = MakieCore.Observable(Matrix{UInt32}(undef, 2 * around[] * (length(points[]) - 1) + 2around[], 3))





    function update_plot(points, n_around)
        vertices = Matrix{Float32}(undef, length(points) * n_around + 2, 3)
        faces = ones(UInt32, 2 * n_around * (length(points) - 1) + 2n_around, 3)

        linespoints = points

        pointsaround = Matrix{MakieCore.Point3f}(undef, n_around, length(linespoints))
        angle = 2π / n_around


        #  beginning and end are not elliptic
        vector = MakieCore.Vec3f(linespoints[2] .- linespoints[1])
        normalvec = normal(vector)
        R = AngleAxis(angle, vector...)
        for j in 1:n_around
            pointsaround[j, 1] = R^(j - 1) * normalvec .* ml.width[] .+ linespoints[1]
        end
        vector = MakieCore.Vec3f(linespoints[end] .- linespoints[end-1])
        normalvec = normal(vector)
        R = AngleAxis(angle, vector...)
        for j in 1:n_around
            pointsaround[j, end] = R^(j - 1) * normalvec .* ml.width[] .+ linespoints[end]
        end



        for i in 2:length(linespoints)-1
            p = linespoints[i]
            prevPoint = linespoints[i-1]
            nextPoint = linespoints[i+1]

            v1 = MakieCore.Vec3f(p .- prevPoint)
            v2 = MakieCore.Vec3f(nextPoint .- p)

            # linearly dependent
            if rank(hcat(v1, v2)) == 1
                normalvec = normal(v1)
                R = AngleAxis(angle, v1...)
                for j in 1:n_around
                    pointsaround[j, i] = R^(j - 1) * normalvec .* ml.width[] .+ p
                end
            else
                #elliptic
                vec, normalvec = intersectPlaneNotZeroAndLinear(v1, v2)
                R = AngleAxis(angle, vec...)
                cosang = cos_angle_between_vectors(v1, vec)
                for j in 1:n_around
                    pointsaround[j, i] = R^(j - 1) * normalvec .* r_ellipse(angle * (j - 1); b=ml.width[], a=ml.width[] / cosang) .+ p
                end
            end
        end










        pointsaroundVector = vec(reshape(pointsaround, (length(pointsaround), 1)))

        rings = [[pointsaround[j, i] for j in 1:n_around] for i in 1:length(linespoints)]
        # Align each ring (starting from the second) with the previous one:
        for i in 2:length(rings)
            rings[i] = align_ring(rings[i-1], rings[i])
        end

        # Flatten the rings into your vertices array (ensure the same order is maintained):
        for (i, ring) in enumerate(rings)
            for j in 1:n_around
                idx = (i - 1) * n_around + j
                vertices[idx, :] .= ring[j]
            end
        end
        vertices[end-1, :] .= points[1]
        vertices[end, :] .= points[end]

        num_rings = length(linespoints)
        face_idx = 1
        for i in 1:(num_rings-1)
            for j in 1:n_around
                next_j = mod1(j + 1, n_around)
                # First triangle of the quad
                faces[face_idx, :] = [(i - 1) * n_around + j,
                    (i - 1) * n_around + next_j,
                    i * n_around + j]
                face_idx += 1
                # Second triangle of the quad
                faces[face_idx, :] = [(i - 1) * n_around + next_j,
                    i * n_around + next_j,
                    i * n_around + j]
                face_idx += 1
            end
        end

        for j in 1:n_around
            next_j = mod1(j + 1, n_around)
            faces[end-2n_around+j, :] = [next_j, j, size(vertices, 1) - 1]
        end

        for j in 1:n_around
            next_j = mod1(j + 1, n_around)
            faces[end-n_around+j, :] = [
                size(vertices, 1) - 2 - n_around + j,
                size(vertices, 1) - 2 - n_around + next_j,
                size(vertices, 1)]
        end

        vertices_ob[] = vertices
        faces_ob[] = faces
    end


    MakieCore.Observables.onany(update_plot, points, around)
    update_plot(points[], around[])

    mesh!(ml, vertices_ob, faces_ob; transparency=false, overdraw=false, color=:red)


    ml
end
