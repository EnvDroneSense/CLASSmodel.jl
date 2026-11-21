# Time stepping and the run loop.

"""
    compute!(s::ModelState; first = false)

Run the components for the current time: the diagnostics, then radiation, surface layer,
land surface, cumulus and mixed layer, as far as they are switched on. With `first` (the
initial state) the surface layer iterates 10 times and the cumulus gets entrainment fluxes.
"""
function compute!(s::ModelState; first::Bool = false)
    x = s.input
    prescribed_fluxes!(s)
    statistics!(s)
    x.sw_rad && radiation!(s)
    if x.sw_sl
        for _ in 1:(first ? 10 : 1)
            surface_layer!(s)
        end
    end
    x.sw_ls && land_surface!(s)
    if x.sw_cu
        first && mixed_layer!(s)   # cumulus! needs the entrainment fluxes of mixed_layer!
        cumulus!(s)
    end
    x.sw_ml && mixed_layer!(s)
    return s
end

"Advance the land surface and the mixed layer by one time step, then the chemistry."
function integrate!(s::ModelState)
    x = s.input
    x.sw_ls && integrate_land_surface!(s)
    x.sw_ml && integrate_mixed_layer!(s)
    s.t += 1
    x.sw_chem && chemistry!(s, x.dt)   # for the new state, at the new time
    return s
end

"""
    step!(s::ModelState)

Advance the model by one time step: compute the components, then integrate the
land surface, the mixed layer and the chemistry. The diagnostics of the step are left in `s`.
"""
function step!(s::ModelState)
    compute!(s)
    integrate!(s)
    return s
end

"""
    run_model(input; constants = PhysicalConstants{T}(), plant = ..., progress = nothing,
        progress_every = 100) -> ModelOutput

Run the model for `floor(runtime / dt)` steps and return the time series. Each row holds
the state before the time step is integrated. The columns are the base variables, then one
mixing ratio [ppb] for each name in `input.species` (see [`ModelOutput`](@ref)). `constants`
and `plant` go to [`init_model`](@ref). With `ls_type = :ags` it throws a `DomainError`
when the vapor pressure deficit is above `f0 / ad` (C3 12.7 kPa, C4 5.7 kPa).

`progress(i, n, out)` is called after every `progress_every` steps and after the last step
`n`. The rows after `i` of `out` are still zero. If it returns `false`, the run stops and the
output holds the rows `1:i`.
"""
function run_model(x::ModelInput{T}; progress = nothing, progress_every::Integer = 100,
        kwargs...) where {T}
    s = init_model(x; kwargs...)
    n = floor(Int, x.runtime / x.dt)
    return run_loop(s, empty_output(T, n, x.species), n, progress, progress_every)
end

# The loop is a separate function so that it is compiled for the concrete type of `out`.
function run_loop(s::ModelState, out, n, progress, progress_every)
    for i in 1:n
        compute!(s)
        store!(out, s)
        integrate!(s)
        if progress !== nothing && (i % progress_every == 0 || i == n) &&
           progress(i, n, out) === false
            foreach(v -> resize!(v, i), out)
            break
        end
    end
    return out
end
