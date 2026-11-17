# Parameter fitting with PEtab.jl

A [`Schedule`](@ref) can be converted into a [PEtab.jl](https://github.com/sebapersson/PEtab.jl) problem, so that parameters of its regulation model can be estimated from experimental measurements. To this end, the stochastic regulation model is approximated by its deterministic (mean-field) ODE, after which PEtab.jl can be used for classical optimisation, profile likelihood methods or Bayesian inference. The fitted parameters can then be written back into the GRS schedule.

## Setup

The extension is loaded automatically when PEtab.jl is loaded together with GeneRegulatorySystems.jl:

```@example petab
using GeneRegulatorySystems
using PEtab
```

## Schedule

A PEtab problem is built around a single `Catalyst.ReactionSystem`, so the schedule must advance through exactly one regulation model.
Instant adjustments such as `{"{add}": ...}` or `{"{set}": ...}` (see [Plumbing](@ref)) are how experiments perturb the system (in PEtab, they become initial values if they happen at time zero, and events otherwise).
Species names must be valid Julia identifiers once `.` is replaced by `₊`, since PEtab parses them in formulas (e.g. `a.proteins*` is not valid).

As an example, we use a schedule where gene `b` activates gene `a`. After a warm-up period, it branches into a `control` phase and an `induced` phase, the latter introducing `b`-proteins into the system.

```@example petab
source = "../../../examples/toy/induction.schedule.json"
read(source, String) |> print
```

```@example petab
schedule! = Models.load(source, seed = "docs")
nothing # hide
```

## Observables and parameters

PEtab.jl parses observable formulas and parameter ids as Julia expressions, so it refers to species and parameters by their ModelingToolkit names, which use `₊` as a delimiter instead of `.` (e.g. `a₊activation₊b₊at`). Model parameters that are not estimated keep the value from the schedule, so only the parameters in question need to be listed.
Multiple model parameters can be tied to a shared estimated parameter via `shared_parameters` (but the shared parameters need to be listed in `parameters`, also).

```@docs
PEtab.PEtabParameter(::Models.Scheduling.Schedule, ::Symbol)
```

Measurements such as qPCR are typically only proportional to molecule counts, so we observe both genes' mRNAs up to a common `scale` factor that we also infer.
Besides the activation threshold, we estimate a single mRNA `decay` rate for both genes:

```@example petab
observables = [
    PEtabObservable(gene, "scale * $(gene)₊mrnas", "0.1 * scale * $(gene)₊mrnas")
    for gene in ("a", "b")
]
parameters = [
    PEtabParameter(schedule!, :a₊activation₊b₊at, scale = :log10, lb = 1.0, ub = 1e5),
    PEtabParameter(:decay, scale = :log10, lb = 1e-4, ub = 1e-1, value = 0.001),
    PEtabParameter(:scale, scale = :log10, lb = 1e-3, ub = 1e1, value = 0.1),
]
shared_parameters = Dict(:a₊mrna_decay => :decay, :b₊mrna_decay => :decay)
nothing # hide
```

## Measurements

Measurements are given as a PEtab measurement table, where `simulation_id` refers to a label of the schedule and `time` is relative to the start of the phase it marks.
PEtab simulates each condition separately, so each label used in the table becomes its own simulation condition, which simulates the label's branch from time zero, including all phases leading up to it.
A label must therefore lie within a single branch and mark a single contiguous stretch of it, so that its start time is unambiguous.

For this example we simulate the data from the model itself.
We first list which genes to measure at which times, build the problem from that list, and fill in the measurements with simulated values plus noise:

```@example petab
using DataFrames
using Random

measurements = DataFrame([
    (obs_id = gene, simulation_id = label, time, measurement = 0.0)
    for gene in ("a", "b")
    for label in ("control", "induced")
    for time in 1200.0:1200.0:10800.0
])
model = PEtabModel(schedule!, observables, measurements, parameters; shared_parameters)
problem = PEtabODEProblem(model)
truth = get_x(problem)
noise = 1 .+ 0.1 .* randn(Xoshiro(1), nrow(measurements))
measurements.measurement = problem.simulated_values(truth) .* noise
filter(row -> row.time <= 2400.0, measurements)
```

## Building the problem

```@docs
PEtab.PEtabModel(::Models.Scheduling.Schedule, ::Any, ::Any, ::Any)
```

With the synthetic data in place, we build the problem again:

```@example petab
model = PEtabModel(schedule!, observables, measurements, parameters; shared_parameters)
problem = PEtabODEProblem(model)
```

## Fitting

The problem can now be fitted with any of the optimizers PEtab.jl supports.
Since the likelihood may have several local optima, we start from multiple random points within the parameter bounds:

```@example petab
using Optim

result = calibrate_multistart(Xoshiro(2), problem, Optim.LBFGS(), 10)
result.xmin, truth
```

which should recover the values the data was simulated from.

## Writing fitted parameters back

```@docs
Models.remake(::Models.Scheduling.Schedule, ::PEtab.EstimationResult, ::PEtab.PEtabODEProblem)
```

Writing the fitted values back into the schedule allows simulating the fitted model with the stochastic simulation and the CLI tools, for example to check whether the mean-field fit also describes the noisy system:

```@example petab
fitted! = Models.remake(schedule!, result, problem)
nothing # hide
```

## Plotting the fit

Since PEtab simulates the whole branch, its time axis starts with the warm-up, and each phase starts later.
We look up each phase's start in the schedule and plot the fitted simulation around it, in hours since that start, together with the data scaled back to molecule counts:

```@example petab
using Plots

segments = Models.Scheduling.segments(schedule!)
scale = 10^result.xmin[:log10_scale]
panels = map(["control", "induced"]) do label
    start = first(s.from for s in segments if s.label == label)
    solution = get_odesol(result, problem, condition = label)
    plot(solution, idxs = [:a₊mrnas, :b₊mrnas], tspan = (start - 3600.0, start + 10800.0), title = label, label = ["a (fitted)" "b (fitted)"],
        xticks = (start .+ 3600.0 .* (-1:3), -1:3), xlabel = "time (h)")
    vline!([start], linestyle = :dash, color = :gray, label = nothing)
    data = filter(row -> row.simulation_id == label, measurements)
    scatter!(start .+ data.time, data.measurement ./ scale, group = data.obs_id, color = [1 2], label = ["a (measured)" "b (measured)"])
end
plot(panels..., layout = (1, 2), size = (900, 350), bottom_margin = 5Plots.mm)
```

## Identifiability

Several parameter regimes can fit the data equally well, meaning the data cannot tell them apart, so that the parameters are *non-identifiable*.
Likelihood profiles from [LikelihoodProfiler.jl](https://github.com/insysbio/LikelihoodProfiler.jl) show how well each parameter is determined by the data.
Here we compare fitting only the `control` data with fitting both the `control` and `induced` data:

```@example petab
using LikelihoodProfiler
using OptimizationLBFGSB

function profiles(measurements)
    model = PEtabModel(schedule!, observables, measurements, parameters; shared_parameters)
    problem = PEtabODEProblem(model)
    result = calibrate_multistart(Xoshiro(2), problem, Optim.LBFGS(), 10)
    method = OptimizationProfiler(optimizer = OptimizationLBFGSB.LBFGSB(), stepper = FixedStep(initial_step = 5e-2))
    LikelihoodProfiler.solve(ProfileLikelihoodProblem(result, problem), method)
end

control = profiles(filter(row -> row.simulation_id == "control", measurements))
both = profiles(measurements)

names = ["log10 a₊activation₊b₊at", "log10 decay", "log10 scale"]
panels = map(eachindex(names)) do i
    panel = plot(xlabel = names[i], ylabel = "Δ objective", ylims = (0, 6), legend = i == 2 ? :topright : false)
    for (profile, label) in ((control, "control only"), (both, "control + induced"))
        plot!(panel, profile[i].x, profile[i].obj .- minimum(profile[i].obj); label, linewidth = 2)
    end
    hline!(panel, [both[i].obj_level - minimum(both[i].obj)], linestyle = :dash, color = :gray, label = "95% threshold")
end
plot(panels..., layout = (1, 3), size = (900, 300), margin = 4Plots.mm)
```

Without the induction, the system stays at steady state, which many combinations of threshold, decay rate and scale explain equally well.
The response to the induction determines all three.

```@raw html
<!-- TODO: write on limitations -->
```
