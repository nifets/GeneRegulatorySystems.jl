# `experiment` tool

```@eval
Main.include_tool_script(@__MODULE__, "experiment")
Main.tool_help(Script.settings())
```

## Results format

The `experiment` tool collects the simulation results into files at the specified location (as set by the `--location` parameter, or by default in a `results/{TIMESTAMP}/` directory). These files are:
* Verbatim copies of the executed schedules' JSON specifications, named like the original files unless there is a naming conflict.
* An `experiment.schedule.json` that references these original specifications and ties them together into a composite top-level experiment schedule with additional meta information.
* An `index.arrow` table in Arrow format that lists all executed simulation segments.
* Any number of event stream files in Allow stream format (with names ending in `.stream.arrow`) containing the actually realized trajectories.

### Experiment specifications

The `experiment.schedule.json` includes information on the versions of Julia and of this package that were used to run the experiment, as well as the used root seed (as set by the `--seed` parameter, or `"seed"` by default).
This can be used to repeat the experiment, and also the other CLI tools will check if the recorded versions agree with those currently in use and emit a warning otherwise.
If more than one specification was passed to the `experiment` tool, the corresponding schedules will be run in separate branches (see [`Schedule`](@ref)).

### Index of simulation segments

The `index.arrow` table (with one row per simulation segment) contains the following columns:
* `i`: `Int64` index of the simulation segment.
* `path`: `String` identifying the simulation segment within the schedule specified by `experiment.schedule.json`; see [Paths and reification](@ref).
* `from` and `to`: `Float64`s recording the (inclusive) start and end time points covered by the simulation segment. They may be equal (that is, an instant segment), and they may overlap between segments (rows) if the simulation was branched.
* `model`: `String` path identifying the primitive model that produced the simulation segment.
* `label`: Custom `String` attached to the simulation segment if `"label"` was set in the schedule specification. This can be used to filter simulation segments, for example to [export](@ref "`export` tool") them.
* `count`: The `Int64` count of events recorded in the event stream referenced by `into` related to this simulation segment. If output was disabled for this segment (because `"into": null` was specified in the schedule), `count` will be `0`.
* `into`: `String` filename (relative to the directory containing the `index.arrow` file) naming the Arrow stream file containing the events recorded for this simulation segment, in long format (see below).

!!! tip
    From Julia, this file can for example be loaded using *Arrow.jl* and used with *DataFrames.jl*:
    ```julia
    using DataFrames
    import Arrow

    d = DataFrame(Arrow.Table("index.arrow"))
    ```

    From Python, this file can for example be loaded using *pyarrow* and used with *pandas*:
    ```python
    import pyarrow
    import pandas

    with open("index.arrow", "rb") as file:
      with pyarrow.ipc.open_file(file) as reader:
            d = reader.read_pandas()
    ```

### Event streams

An `events.stream.arrow` stream (with one row per event) contains the following columns:

* `i`: `Int64` foreign key into `index.arrow` referencing the simulation segment that produced this event.
* `t`: `Float64` holding the time point of this event (in simulation time, ignoring any branching).
* `name`: `Symbol` naming the molecular species changing through this event (potentially using `.` as a subsystem separator).
* `value`: The new `Int64` value after this event.

The trajectories are thus recorded in long format.
For each simulation segment, if events were recorded in the first place, the stream will contain an initializing event for each molecular species and then one event for each jump in the trajectory.

Typically, multiple events belong to the same simulation segment, which can be joined on `i` with `index.arrow` to access e.g. the segment interval.
However, it is sometimes possible to avoid such joins by specifying the schedule to break out events into multiple events streams (that is, setting `"into"` appropriately); see [`Schedule`](@ref).
For example it may make sense to use one stream per simulation branch.

!!! tip
    From Julia, the `events.stream.arrow` files can be used like described above. From Python, `pyarrow.ipc.open_stream` must be used instead of `pyarrow.ipc.open_file` because it does not automatically detect the streaming results format.
