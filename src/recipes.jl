@MakieCore.recipe(MeshLines, points) do scene
    MakieCore.Attributes(
        around = 3,
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

function sin_angle_between_vectors(a, b)
    return norm(cross(a, b)) / (norm(a) * norm(b))
end

function r_ellipse(ϕ;a=1,b=1)
    return a*b/sqrt(a^2*sin(ϕ)^2+b^2*cos(ϕ)^2)
end


function MakieCore.plot!(ml::MeshLines{<:Tuple{<:AbstractVector{<:MakieCore.Point3f}}})
    linespoints = ml[1][]
    n_around = ml.around[]
    pointsaround = Matrix{MakieCore.Point3f}(undef,n_around,length(linespoints))
    angle = 2π/n_around


    #  beginning and end are not elliptic
    vector = MakieCore.Vec3f(linespoints[2] .- linespoints[1])
    normalvec = normal(vector)
    R = AngleAxis(angle, vector...)
    for j in 1:n_around
        pointsaround[j,1] = R^(j-1)*normalvec.*ml.width[] .+ linespoints[1]
    end
    vector = MakieCore.Vec3f(linespoints[end] .- linespoints[end-1])
    normalvec = normal(vector)
    R = AngleAxis(angle, vector...)
    for j in 1:n_around
        pointsaround[j,end] = R^(j-1)*normalvec.*ml.width[] .+ linespoints[end]
    end



    for i in 2:length(linespoints)-1
        p = linespoints[i]
        prevPoint = linespoints[i-1]
        nextPoint = linespoints[i+1]

        v1 = MakieCore.Vec3f(p .- prevPoint)
        v2 = MakieCore.Vec3f(nextPoint .- p)

        # linearly dependent
        if rank(hcat(v1,v2)) == 1
            normalvec = normal(v1)
            R = AngleAxis(angle, v1...)
            for j in 1:n_around
                pointsaround[j,1] = R^(j-1)*normalvec.*ml.width[] .+ p
            end
        else
            vec,normalvec = intersectPlaneNotZeroAndLinear(v1,v2)
            R = AngleAxis(angle, vec...)
            sinang = sin_angle_between_vectors(v1, vec)
            println(sinang)
            for j in 1:n_around
                pointsaround[j,i] = R^(j-1)*normalvec.*r_ellipse(angle*(j-1);b=ml.width[],a=ml.width[]/sinang) .+ p
            end
        end
    end

    display(pointsaround)
    pointsaroundVector = vec(reshape(pointsaround,(length(pointsaround),1)))
    vertices = Matrix{Float32}(undef,length(pointsaround),3)
    for (i,p) in enumerate(pointsaroundVector)
        vertices[i,:] .= p
    end
    # faces = Matrix{Int32}(undef,length(linespoints)*(n_around-1),3)
    faces = Matrix{Int32}(undef,4n_around,3)
    for j in 1:n_around
        faces[j,1] = j
        faces[j,2] = mod1(j+1,n_around)
        faces[j,3] = mod1(j+1,n_around)+n_around
    end
    # for i in 1:length(linespoints)-1
    for j in 1:n_around
        faces[j+n_around,1] = j + n_around
        faces[j+n_around,2] = mod1(j+1,n_around) + n_around
        faces[j+n_around,3] = mod1(j,n_around)

        faces[j+2n_around,1] = j + n_around
        faces[j+2n_around,2] = mod1(j+1,n_around) + n_around
        faces[j+2n_around,3] = mod1(j,n_around) + 2n_around
    end
    for j in 1:n_around
        faces[j+3n_around,1] = j+ 2n_around
        faces[j+3n_around,2] = mod1(j+1,n_around) + 2n_around
        faces[j+3n_around,3] = mod1(j+1,n_around)+n_around
    end

    mesh!(ml, vertices, faces; transparency=false, overdraw=false)

    display(vertices)
    meshscatter!(ml,vec(reshape(pointsaround,(length(pointsaround),1))), markersize = 0.01)
    meshscatter!(ml,vec(pointsaround[1,:]), markersize = 0.01,color=:red)
    lines!(ml,ml[1])
    ml
end
