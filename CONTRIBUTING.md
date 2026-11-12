# Contributing

Thank you for your help with CLASSmodel.jl. This note gives the rules for a change.

## Layout

| Folder | Content |
|:-------|:--------|
| `src/io/` | `ModelInput`, `validate`, `input_info`, the output table, `output_info` |
| `src/core/` | `ModelState` and the time loop (`init_model`, `step!`, `run_model`) |
| `src/physics/` | one file for each component: mixed layer, radiation, surface layer, land surface, A-Gs, cumulus, chemistry |
| `test/` | one test file for each component |
| `docs/` | the documentation (Documenter.jl), one page for each topic in `docs/src/` |

## Tests

Run the full test suite before you commit. It includes the Aqua and JET checks.

```bash
julia --project=. -e "using Pkg; Pkg.test()"
```

For a small change, run one test file. This command runs the chemistry tests:

```bash
julia --project=. -e "using CLASSmodel, Test; include(\"test/test_chemistry.jl\")"
```

Write one test for each real behaviour. Do not add a test that repeats another test.

## Code style

- CI formats the code with JuliaFormatter 1. The settings are in `.JuliaFormatter.toml`. Format before you commit:

  ```bash
  julia -e "using JuliaFormatter; format(\".\")"
  ```

- Keep the line endings of the files you edit. Git stores LF. A Windows checkout can use CRLF.
- Write one line for the non-obvious reason in a comment. Use at most two more lines.
- Keep docstrings short. Document a public name with its signature and one short description.

## CHANGELOG

A change of behaviour needs an entry in `CHANGELOG.md`, under `[Unreleased]`. Use the section
names given at the top of that file. A change that only touches the docs or the tests needs no entry.

## Add an input

1. Add the field to `ModelInput` in `src/io/input.jl`, with a default, a unit and a one-line comment.
2. Add a rule to `validate` if the value can be wrong.
3. Add an entry to `INPUT_INFO` in `src/io/input_info.jl`. The keys must follow the order of the `ModelInput` fields.
4. Add the field to the table in the `ModelInput` docstring.
5. Use the field in the physics code.

## Add an output

1. Add the name to `OUTPUT_NAMES` in `src/io/output.jl`.
2. Write the value in `store!`, in the same file.
3. Add an entry to `BASE_INFO` in `src/io/output_info.jl`. The keys must follow the order of `OUTPUT_NAMES`.

## Add a switch or a physics module

1. Add the `sw_` field to `ModelInput`, as for an input. The default is `false`.
2. Put the physics in a file in `src/physics/`. Add it to the `include` list in `src/CLASSmodel.jl`.
3. Add the state it needs to `src/core/state.jl`.
4. Call the new function from `compute!` or `integrate!` in `src/core/model.jl`.
5. Add a test file in `test/` and include it in `test/runtests.jl`.
6. Add the switch to the table in `docs/src/switches.md`.

## Documentation

Build the docs from the repository root:

```bash
julia --project=docs -e "using Pkg; Pkg.develop(path=\".\"); Pkg.instantiate(); include(\"docs/make.jl\")"
```
