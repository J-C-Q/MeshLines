using GLMakie
using MeshLine

points = Observable(Point3f[(x, sin(x), cos(x)) for x in range(0, 4pi, 1000)])

fig = Figure()

scene = LScene(fig[1, 1])

slider = Slider(fig[2, 1], value=1, range=range(0.1, 10, 1000))

on(slider.value) do value
    points[] = Point3f[(x, sin(x * value), cos(x * value)) for x in range(0, 4pi, 1000)]
end
# points[] = Point3f[(x, sin(x * 5), cos(x * 5)) for x in range(0, 4pi, 100)]
meshlines!(scene, points)

display(fig)
