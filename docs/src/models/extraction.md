# Extraction

The package provides various [molecular counts resampling models](@ref "Counts resampling primitives") that can be assembled into observation models to mimic experimental protocols for extracting and measuring materials from cells.
A small set of such extraction schemes [is predefined](@ref "Predefined extraction schemes").
Although it is not a requirement, they are `Instant` and directly [change the
system state](@ref "Non-destructive extraction").
For an example, see `examples/specification/extraction.schedule.json`.

More realistic extraction schemes need to be defined by hand.

## Predefined extraction schemes

```@docs
Models.Extraction
Models.Extraction.simple_transcriptome
Models.Extraction.amplified_transcriptome
Models.Extraction.simple_proteome
```

Since these models are implemented as [`Schedule`](@ref)s, which normally call back `trace` to produce output for each simulation segment, but here we would like to treat the extraction as a unitary step and are not interested in the intermediate steps, we wrap the specification to suppress that output:

```@docs
Models.Extraction.with_intermediate_output_suppressed
```

## Non-destructive extraction

Invoking an extraction scheme will directly modify the current state of the simulated system.
This corresponds to the assumption that extraction from a real system (like in scRNAseq) would physically destroy it.
Extraction therefore typically ends a sampled trajectory.

To model non-destructive observation that does not affect the trajectory, a simulation schedule may be instructed to branch before the application of the observation. For example,
```
[
    <pre>,
    {"branch": true, "step": [
        <extraction>
    ]},
    <post>
]
```
will run the `<pre>` step(s), then branch off to apply an instant `<extraction>` step, and then return to the stem and proceed with the `<post>` step(s). The branched model may itself be a [`Schedule`](@ref) and thus for example continue regulation for a while before the extraction, and it is also possible to simulate multiple extractions from the same state, or to regularly branch-and-extract until the simulation time budget is exhausted.

## Combining multiple modalities

Similarly, it is possible to invoke multiple distinct extraction schemes on separate branches and to then merge the results back together. This is written as
```
[
    {"branch": true, "step": [
        <extraction1>,
        <extraction2>
    ]},
    {"{merge}": "+"}
]
```
and it can for example be used to simulate a naive multi-omics protocol; see also [`Models.Scheduling.Merge`](@ref).

## Counts resampling primitives

For reference, these are the resampling primitives that are used to construct extraction schemes:

```@autodocs
Modules = [Models.Resampling]
```
