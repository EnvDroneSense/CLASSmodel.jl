# Time stepping and the run loop.

"""
    compute!(s::ModelState)

Run the components for the current time: the diagnostics, then radiation, surface layer,
land surface, cumulus and mixed layer, as far as they are switched on.
"""
function compute!(s::ModelState)
    x = s.input
    prescribed_fluxes!(s)
    statistics!(s)
    x.sw_rad && radiation!(s)
    x.sw_sl && surface_layer!(s)
    x.sw_ls && land_surface!(s)
    x.sw_cu && cumulus!(s)
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
land surface and the mixed layer. The diagnostics of the step are left in `s`.
"""
function step!(s::ModelState)
    compute!(s)
    integrate!(s)
    return s
end

"""
    run_model(input; constants = PhysicalConstants{T}(), progress = nothing,
        progress_every = 100) -> ModelOutput

Run the model for `floor(runtime / dt)` steps and return the time series. Each row holds
the state before the time step is integrated. With `ls_type = :ags` it throws a `DomainError`
when the vapour pressure deficit is above `f0 / ad` (C3 12.7 kPa, C4 5.7 kPa).

`progress(i, n, out)` is called after every `progress_every` steps and after the last step
`n`. The rows after `i` of `out` are still zero. If it returns `false`, the run stops and the
output holds the rows `1:i`.
"""
function run_model(x::ModelInput{T}; progress = nothing, progress_every::Integer = 100,
        kwargs...) where {T}
    s = init_model(x; kwargs...)
    n = floor(Int, x.runtime / x.dt)
    out = empty_output(T, n)
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
