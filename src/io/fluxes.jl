# Prescribed surface fluxes: a number, or a shape that gives the flux as a function of time.

"""
    FluxShape

A surface flux that changes in time. `wtheta`, `wq`, `wCO2` and each entry of `wC` take a
number (constant flux) or one of [`Sine`](@ref), [`Cosine`](@ref), [`Block`](@ref),
[`Table`](@ref), [`Sun`](@ref), [`Deposition`](@ref). An amplitude has the unit of the flux.
A new shape needs one method `CLASSmodel.flux(shape, x::ModelInput, t, c)`, with `t` in
seconds from the start of the run and `c` the mixed-layer value of the quantity.
"""
abstract type FluxShape end

"""
    Sine(amplitude; period = 43200, start = 0)

`amplitude * sin(pi (t - start) / period)` between `start` and `start + period` [s], else zero.
"""
struct Sine{T <: Real} <: FluxShape
    amplitude::T
    period::T
    start::T
end

"""
    Cosine(amplitude; period = 43200, start = 0)

`amplitude / 2 * (1 - cos(2 pi (t - start) / period))` between `start` and `start + period` [s],
else zero.
"""
struct Cosine{T <: Real} <: FluxShape
    amplitude::T
    period::T
    start::T
end

"""
    Block(amplitude; period = 43200, start = 0)

`amplitude` between `start` and `start + period` [s], else zero.
"""
struct Block{T <: Real} <: FluxShape
    amplitude::T
    period::T
    start::T
end

"""
    Table(values; interval = 3600)

One value every `interval` [s] from the start of the run, with linear interpolation. After the
last value the flux stays at that value.
"""
struct Table{T <: Real} <: FluxShape
    values::Vector{T}
    interval::Float64
end

"""
    Sun(amplitude)

`amplitude` times the sine of the solar elevation divided by its noon value. Zero at night.
"""
struct Sun{T <: Real} <: FluxShape
    amplitude::T
end

"""
    Deposition(vd)

Dry deposition: the flux is `-vd * c`, with the deposition velocity `vd` [m s-1].
"""
struct Deposition{T <: Real} <: FluxShape
    vd::T
end

for S in (:Sine, :Cosine, :Block)
    @eval $S(amplitude; period = 43200, start = 0) = $S(promote(float(amplitude), period, start)...)
end
function Table(values; interval = 3600)
    interval > 0 || throw(ArgumentError("Table needs interval > 0"))
    return Table(float.(collect(values)), float(interval))
end

# Throw if a number inside a field of an input is not finite (a Bool and other types pass).
check_finite(name, v::Real) = v isa Bool || isfinite(v) ||
    throw(ArgumentError("$name must be finite, got $v"))
check_finite(name, v::Union{Tuple, AbstractVector}) = foreach(w -> check_finite(name, w), v)
check_finite(name, v::FluxShape) = foreach(n -> check_finite(name, getfield(v, n)), fieldnames(typeof(v)))
check_finite(name, v) = nothing

# The value of a flux or lapse-rate field of an input with number type `T`: a number becomes a `T`.
fieldvalue(::Type{T}, v::Real) where {T} = T(v)
fieldvalue(::Type{T}, v::Union{FluxShape, Segments{T}}) where {T} = v
fieldvalue(::Type{T}, v::Segments) where {T} = Segments{T}(v.z, v.gamma)
# The same for a field of type `F`: only a field that takes `T` or a shape is converted.
function fieldvalue(::Type{T}, ::Type{F}, v) where {T, F}
    F !== T && T <: F && return fieldvalue(T, v)
    F <: Tuple && v isa Tuple && return map(w -> fieldvalue(T, eltype(F), w), v)
    F <: AbstractVector && eltype(F) <: Union{Real, FluxShape} && v isa AbstractVector &&
        return eltype(F)[fieldvalue(T, eltype(F), w) for w in v]
    return v
end

# Fraction of the window that has passed, or nothing outside the window.
function window(f, t)
    p = (t - f.start) / f.period
    return 0 <= p <= 1 ? p : nothing
end

"""
    flux(shape, x::ModelInput, t, c)

Flux of `shape` at `t` [s] from the start of the run, for the mixed-layer value `c`.
A number is a constant flux.
"""
flux(v::Real, x, t, c) = v
function flux(f::Sine, x, t, c)
    window(f, t) === nothing && return zero(c)
    return f.amplitude * sin(pi / f.period * (t - f.start))
end
function flux(f::Cosine, x, t, c)
    p = window(f, t)
    return p === nothing ? zero(c) : f.amplitude / 2 * (1 - cospi(2p))
end
flux(f::Block, x, t, c) = window(f, t) === nothing ? zero(c) : f.amplitude
function flux(f::Table, x, t, c)
    i = floor(Int, t / f.interval) + 1
    i >= length(f.values) && return f.values[end]
    p = t / f.interval - (i - 1)
    return f.values[i] + p * (f.values[i + 1] - f.values[i])
end
function flux(f::Sun, x, t, c)
    return f.amplitude * max(solar_elevation(x, utc_seconds(x, t)), 0) / solar_elevation_max(x)
end
flux(f::Deposition, x, t, c) = -f.vd * c
