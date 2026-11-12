# Running the model

`run_model` runs the whole period. `init_model` and `step!` give control over each time step.

```@setup quickstart
using CLASSmodel
```

## Step by step

```@example quickstart
s = init_model(ModelInput(; runtime = 600))
step!(s)
s.mixed_layer.h
```

## Progress and early stop

`run_model` calls `progress(i, n, out)` after every `progress_every` steps and after the last
step `n`. If the callback returns `false`, the run stops. The output then holds the rows `1:i`.

```@example quickstart
stop_at_one_hour(i, n, out) = out.t[i] - out.t[1] < 1
out = run_model(ModelInput(); progress = stop_at_one_hour, progress_every = 10)
length(out.t)
```
