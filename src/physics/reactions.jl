# The gas-phase reactions of the chemistry scheme and the production and loss terms built
# from them.

"""
    Reaction(name, form, inputs, outputs, k)

One reaction of the chemistry. `inputs` and `outputs` map a species to its coefficient (at
most two inputs). `k` holds the constants A to G of the rate, which `form` selects:

| `form` | rate constant |
|:-------|:--------------|
| `:photolysis` | `cloud * A exp(B / cos(zenith))`, zero at night |
| `:constant` | `A` |
| `:arrhenius` | `A exp(B / T)` |
| `:falloff`, `:decomposition` | pressure fall-off from `A` to `G` (two-body, and the decomposition of N2O5) |
| `:ho2` | HO2 + HO2, with a water term |
| `:power` | `A (T / B)^C exp(D / T)` |

The input `reactions` of `ModelInput` is
a vector of these; the default is `CLASSmodel.REACTIONS`. `phi` and `k_r05` in the output use
the reactions `:r05` and `:r21`.
"""
struct Reaction
    name::Symbol
    form::Symbol
    inputs::Vector{Pair{Symbol, Int}}
    outputs::Vector{Pair{Symbol, Float64}}
    k::NTuple{7, Float64}
end

const FORMS = (:photolysis, :constant, :arrhenius, :falloff, :decomposition, :ho2, :power)

# Throw an ArgumentError for a reaction that the scheme cannot build or evaluate.
function check_reaction(rx::Reaction, species)
    rx.form in FORMS || throw(ArgumentError("reaction $(rx.name): form must be one of $FORMS"))
    length(rx.inputs) <= 2 ||
        throw(ArgumentError("reaction $(rx.name) has more than two input species"))
    for (n, _) in Iterators.flatten((rx.inputs, rx.outputs))
        n in species || throw(ArgumentError("reaction $(rx.name): unknown species $n"))
    end
end

# The table of the C++ model. The coefficient of OH in r19 is its yield, zero by default.
const REACTIONS = [
    Reaction(:r01, :photolysis, [:O3 => 1], [:O1D => 1.0, :O2 => 1.0],
        (3.171e-4, -1.809, 1.2, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r02, :arrhenius, [:O1D => 1, :H2O => 1], [:OH => 2.0],
        (1.63e-10, 60.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r03, :arrhenius, [:O1D => 1, :N2 => 1], [:O3 => 1.0],
        (2.15e-11, 110.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r04, :arrhenius, [:O1D => 1, :O2 => 1], [:O3 => 1.0],
        (3.3e-11, 55.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r05, :photolysis, [:NO2 => 1], [:NO => 1.0, :O3 => 1.0],
        (0.0193, -0.575, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r06, :photolysis, [:CH2O => 1], [:HO2 => 1.0],
        (2.16e-4, -0.778, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r07, :constant, [:OH => 1, :CO => 1], [:HO2 => 1.0],
        (2.4e-13, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r08, :arrhenius, [:OH => 1, :CH4 => 1], [:CH3O2 => 1.0],
        (2.45e-12, -1775.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r09, :constant, [:OH => 1, :ISO => 1], [:RO2 => 1.0],
        (1e-10, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r10, :constant, [:OH => 1, :MVK => 1], [:HO2 => 1.0, :CH2O => 1.0],
        (2.4e-11, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r11, :arrhenius, [:OH => 1, :HO2 => 1], [:H2O => 1.0, :O2 => 1.0],
        (4.8e-11, 250.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r12, :arrhenius, [:OH => 1, :H2O2 => 1], [:H2O => 1.0, :HO2 => 1.0],
        (2.9e-12, -160.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r13, :arrhenius, [:HO2 => 1, :NO => 1], [:OH => 1.0, :NO2 => 1.0],
        (3.5e-12, 250.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(
        :r14, :arrhenius, [:CH3O2 => 1, :NO => 1], [:HO2 => 1.0, :NO2 => 1.0, :CH2O => 1.0],
        (2.8e-12, 300.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r15, :constant, [:RO2 => 1, :NO => 1],
        [:HO2 => 1.0, :NO2 => 1.0, :MVK => 1.0, :CH2O => 1.0],
        (1e-11, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r16, :arrhenius, [:OH => 1, :CH2O => 1], [:HO2 => 1.0, :CO => 1.0],
        (5.5e-12, 125.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r17, :ho2, [:HO2 => 2], [:H2O2 => 1.0],
        (2.2e-13, 600.0, 1.9e-33, 980.0, 1.4e-21, 2200.0, 1.0)),
    Reaction(:r18, :arrhenius, [:CH3O2 => 1, :HO2 => 1], [:Product => 1.0],
        (4.1e-13, 750.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r19, :constant, [:RO2 => 1, :HO2 => 1], [:OH => 0.0, :Product => 1.0],
        (1.5e-11, 750.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r20, :arrhenius, [:OH => 1, :NO2 => 1], [:HNO3 => 1.0],
        (3.5e-12, 340.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r21, :arrhenius, [:NO => 1, :O3 => 1], [:NO2 => 1.0],
        (3e-12, -1500.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r22, :arrhenius, [:NO => 1, :NO3 => 1], [:NO2 => 2.0],
        (1.8e-11, 110.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r23, :arrhenius, [:NO2 => 1, :O3 => 1], [:NO3 => 1.0, :O2 => 1.0],
        (1.4e-13, -2470.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r24, :falloff, [:NO2 => 1, :NO3 => 1], [:N2O5 => 1.0],
        (3.6e-30, -4.1, 0.0, 1.9e-12, 0.2, 0.0, 0.35)),
    Reaction(:r25, :decomposition, [:N2O5 => 1], [:NO3 => 1.0, :NO2 => 1.0],
        (1.3e-3, -3.5, -11000.0, 9.7e14, 0.1, -11080.0, 0.35)),
    Reaction(:r26, :arrhenius, [:N2O5 => 1, :H2O => 1], [:HNO3 => 2.0],
        (2.5e-22, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0)),
    Reaction(:r27, :power, [:N2O5 => 1, :H2O => 2], [:HNO3 => 2.0, :H2O => 1.0],
        (1.8e-39, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0)),
    Reaction(:r28, :arrhenius, [:OH => 1, :O3 => 1], [:HO2 => 1.0, :O2 => 1.0],
        (1.7e-12, -940.0, 1.0, 1.0, 1.0, 1.0, 1.0))
]

# Species that the solver leaves alone: the inert gases, H2O (set from q) and a pure product.
const NOT_SOLVED = (:N2, :O2, :H2O, :Product)

# One term of the production or the loss of a species: `coef * K[reaction] * y[i1]^e1 * y[i2]^e2`.
struct Term{T}
    production::Bool
    reaction::Int
    coef::T
    i1::Int
    e1::Int
    i2::Int
    e2::Int
end

# The reactions that run and, for each species, the terms of its production and loss.
struct ChemistryScheme{T}
    species::Vector{Symbol}
    reactions::Vector{Reaction}
    terms::Vector{Vector{Term{T}}}
    solved::Vector{Int}
end

"""
    ChemistryScheme{T}(species, reactions)

Build the production and loss terms of `reactions` for the list of `species`. The loss of an input species is the
reaction rate divided by that species, so that the solver can treat it implicitly.
"""
function ChemistryScheme{T}(species::Vector{Symbol}, reactions::Vector{Reaction}) where {T}
    index(n) = something(findfirst(==(n), species))
    terms = [Term{T}[] for _ in species]
    function term(production, r, coef, factors)
        f = filter(>(0) ∘ last, factors)
        f = [f; fill((1, 0), 2 - length(f))]
        return Term{T}(production, r, coef, f[1][1], f[1][2], f[2][1], f[2][2])
    end
    for (r, rx) in enumerate(reactions)
        inputs = [(index(n), c) for (n, c) in rx.inputs]
        for (k, (i, c)) in enumerate(inputs)
            factors = [(j, kk == k ? e - 1 : e) for (kk, (j, e)) in enumerate(inputs)]
            push!(terms[i], term(false, r, T(c), factors))
        end
        for (n, c) in rx.outputs
            push!(terms[index(n)], term(true, r, T(c), inputs))
        end
    end
    solved = [i for i in eachindex(species) if !isempty(terms[i]) && !(species[i] in NOT_SOLVED)]
    return ChemistryScheme{T}(species, reactions, terms, solved)
end

# Index of species `name` in the scheme, or 0 when the scheme does not have it.
function species_index(scheme::ChemistryScheme, name)
    i = findfirst(==(name), scheme.species)
    return i === nothing ? 0 : i
end
