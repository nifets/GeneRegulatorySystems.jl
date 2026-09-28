# GeneRegulatorySystems.jl

[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)][docs]

A package to construct models of single-cell gene regulatory dynamics and orchestrate their simulation.

## Overview

This package contains components to construct complex models of single-cell gene regulation from templates and can coordinate reproducible simulation through a flexible scheduling mechanism and output format. The core regulation models are currently defined using [Catalyst.jl](https://github.com/SciML/Catalyst.jl) and simulated via [JumpProcesses.jl](https://github.com/SciML/JumpProcesses.jl), allowing them to be treated in isolation with various tools from the SciML ecosystem. Scheduling supports various forms of instant adjustment and fine-grained output control. See [Overview][docs-overview] for details.

Features include:

* Construction of Catalyst.jl `ReactionSystem`-based gene regulation models using a simple JSON-based template language.
* Simulation using exact stochastic methods from JumpProcesses.jl.
* Experiment scheduling mechanism that supports complex observation models, seed control and simulation branching, as well as persistent and transient interventions (including periodic interventions).
* Export to counts matrices for downstream applications.

## Setup

The package can be used either as a library or as a standalone application via the command line. It requires a working installation of Julia ≥ v1.11.

As a library, it can be installed using Julia's package manager. For example, type
```
]add https://github.com/drostlab/GeneRegulatorySystems.jl
```
in the REPL.

With Git installed, the easiest way to get started using the command line interface is to clone this repository, e.g. like
```sh
git clone https://github.com/drostlab/GeneRegulatorySystems.jl $HOME/src/grs
```
and then run the utilities via the (currently Unix-only) wrapper script. From the cloned directory, run
```sh
tools/grs
```
to see the included utilities.

> [!NOTE]
> This assumes that a recent Julia is available as `julia`; you may alternatively override the used Julia binary by setting the `JULIA` environment variable, for example by prepending it like
> ```sh
> JULIA=/path/to/julia tools/grs
> ```
> for a single invocation.

The wrapper will automatically run the utilities in individual Julia environments with their respective dependencies, which may take a couple of minutes on each first use as they are precompiled.

> [!TIP]
> For convenience, you could optionally create a shortcut to `tools/grs`, for example by defining an `alias` in your `.bashrc` like
> ```sh
> echo "alias grs='\$HOME/src/grs/tools/grs'" >> "$HOME/.bashrc"
> ```
> or by adding `$HOME/src/grs/tools` to your `PATH` environment variable.

On Windows, since the wrapper is not available you need to run the tools manually by setting the correct Julia environment. For example, you can run
```sh
julia --project=tools/experiment tools/experiment/run.jl
```
from the cloned directory.

## Quickstart

Simulation experiments are specified by schedules that define the models to be simulated and how their time evolution should be concatenated along the time dimension. These schedules may be assembled manually or loaded from specification files written in a JSON-based mini-language. See the `examples` directory for a collection of example schedules.

The following assumes that you have made the `tools/grs` wrapper script available as `grs` as described above. From the repository's directory, run
```sh
grs experiment examples/toy/repressilator.schedule.json
```
to run that schedule. Outputs are by default placed in a `results` subdirectory of the current working directory, which can be changed using the `--location` argument when invoking `grs experiment`. Run `grs experiment --help` for more options, and see [Experiment tool][docs-experiment] for details.

The results consist of

* the schedule(s) that were run to produce them (listed in `experiment.schedule.json`),
* the simulated event stream(s) (ending in `.stream.arrow`) in long format and
* the index of all executed simulation segments (`index.arrow`).

You can inspect Arrow files e.g. like
```julia
using DataFrames
import Arrow

DataFrame(Arrow.Table("path/to/results/index.arrow"))
```
both for the streaming and non-streaming variant. See [Results format][docs-experiment-results] for details.

There is an experimental Makie-based `grs inspect` tool that can display the results and supports some interactivity. Run it like
```sh
grs inspect path/to/results/
```
and see [Inspect tool][docs-inspect] for details.

You can export results from their native (long) format to a wide format using the `grs export` tool. Its default behavior is to include the values of all trajectory dimensions at the end of the last simulation segment for each unique label, although both the included dimensions and considered segments can be pre-filtered; run `grs export --help` for options, and see [Export tool][docs-export] for details.

## Updating

When used as a library, the package can just be updated using Julia's package manager. If it was installed by cloning this repository as described in [Setup](#setup), you can update by checking out the new version using Git. The `tools/grs` wrapper script will automatically re-instantiate the changed tools' environments on first use.

[docs]: https://drostlab.github.io/GeneRegulatorySystems.jl/stable/
[docs-overview]: https://drostlab.github.io/GeneRegulatorySystems.jl/stable/overview/
[docs-experiment]: https://drostlab.github.io/GeneRegulatorySystems.jl/stable/tools/experiment/
[docs-experiment-results]: https://drostlab.github.io/GeneRegulatorySystems.jl/stable/tools/experiment/#Results-format
[docs-inspect]: https://drostlab.github.io/GeneRegulatorySystems.jl/stable/tools/inspect/
[docs-export]: https://drostlab.github.io/GeneRegulatorySystems.jl/stable/tools/export/
