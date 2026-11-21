# API

```@meta
CurrentModule = CLASSmodel
```

| Name | Purpose |
|:-----|:--------|
| [`ModelInput`](@ref) | all input of a run, checked by [`validate`](@ref) |
| [`run_model`](@ref) | run the model, return a [`ModelOutput`](@ref) table |
| [`init_model`](@ref), [`step!`](@ref) | run the model one step at a time on a [`ModelState`](@ref) |
| [`input_info`](@ref), [`variable_info`](@ref) | title, unit and symbol of an input or an output |
| [`mixed_layer_profile`](@ref) | theta and q profiles at one output row |
| [`FluxShape`](@ref) | surface flux that changes in time: [`Sine`](@ref), [`Cosine`](@ref), [`Block`](@ref), [`Table`](@ref), [`Sun`](@ref), [`Deposition`](@ref) |
| [`PhysicalConstants`](@ref) | physical constants, passed with `constants` |
| [`AGsParameters`](@ref), [`PlantType`](@ref), [`C3`](@ref), [`C4`](@ref) | plant parameters of the A-Gs scheme |

## Module

```@docs
CLASSmodel
```

## Input

```@docs
ModelInput
validate
input_info
FluxShape
Sine
Cosine
Block
Table
Sun
Deposition
PhysicalConstants
AGsParameters
PlantType
C3
C4
```

## Run

```@docs
run_model
init_model
step!
ModelState
```

## Output

```@docs
ModelOutput
variable_info
mixed_layer_profile
```

## Index

```@index
```
