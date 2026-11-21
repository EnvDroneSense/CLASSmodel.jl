# Lapse rate that changes with height.

"""
    Segments(z, gamma)

A lapse rate that changes with height: `gamma[i]` holds below `z[i]` [m] and `gamma[end]`
above `z[end]`, so `gamma` has one entry more than the ascending `z`. `gammatheta`, `gammaq`,
`gammaCO2`, `gammau` and `gammav` take a number or `Segments`.
"""
struct Segments{T <: Real}
    z::Vector{T}
    gamma::Vector{T}
    function Segments{T}(z, gamma) where {T}
        issorted(z) && length(gamma) == length(z) + 1 ||
            throw(ArgumentError("Segments needs ascending z and one more gamma than z"))
        return new{T}(z, gamma)
    end
end
Segments(z, gamma) = Segments{Float64}(z, gamma)

"Lapse rate just above the height `h`."
lapse(g::Real, h) = g
lapse(g::Segments, h) = g.gamma[searchsortedlast(g.z, h) + 1]

# Increase of a quantity with the lapse rate `g` from the height `h` to `ztop`.
rise(g::Real, h, ztop) = g * (ztop - h)
function rise(g::Segments, h, ztop)
    edges = [h; filter(z -> h < z < ztop, g.z); ztop]
    return sum(lapse(g, edges[i]) * (edges[i + 1] - edges[i]) for i in 1:(length(edges) - 1))
end
