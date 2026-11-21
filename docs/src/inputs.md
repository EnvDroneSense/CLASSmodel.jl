# Inputs

```@setup quickstart
using CLASSmodel
```

## Input information

`input_info(name)` returns the title, unit, symbol and slider range of a `ModelInput` field.
The slider range is not a validity rule. `validate` checks the input when `ModelInput` is
constructed and throws an `ArgumentError` if the input cannot be run.

```@example quickstart
input_info(:wtheta)
```

## Input fields

The table comes from `input_info` and the defaults of `ModelInput()`. A tuple field holds
one value for each species or reaction. Its row gives the value of one element.

```@eval
using CLASSmodel, Markdown
x = ModelInput()
fmt(v::Tuple) = "$(length(v)) values"
fmt(v) = "`$(repr(v))`"
rows = map(fieldnames(ModelInput)) do name
    info = input_info(name)
    "| `$name` | $(info.title) | $(info.unit) | $(fmt(getfield(x, name))) |"
end
Markdown.parse(join(["| Field | Title | Unit | Default |", "|:--|:--|:--|:--|", rows...], "\n"))
```
