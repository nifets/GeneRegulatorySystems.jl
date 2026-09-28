module Script

using ArgParse

DESCRIPTION = """
Convert the results produced by the experiment tool to wide format (counts
matrices).
"""

EPILOG = """
This will produce a counts matrix with trajectory components (dimensions) as
columns, and observations as rows. One observation is produced for each distinct
label among the included simulation segments, with values determined at the end
of the final corresponding segment. By default, this means that intermediate
observations are only possible if they are distinctly labeled. If this was not
specified in the simulation schedule, the `--relabel-with-segment-id` option may
be used to uniquely relabel all simulation segments.



The segments and trajectory components to include can (and typically need to) be
filtered, for example by selecting only the segments corresponding to extraction
models. This is easiest to achieve by already preparing the simulation schedule
appropriately. Besides labeling the segments, the schedule might redirect any
targeted segments into a separate events stream by defining `"into"`. Segments
can then simply be chosen by the `--model`, `--events` and `--labeled` options.
Likewise, the trajectory dimensions to include in observations can be selected
through the `--kinds` option.



If the identity of an observation's final segment is ambiguous (because of
unresolved trajectory branching), an error will be raised. This can be
side-stepped through the `--relabel-with-segment-id` option. Dimensions that
have no events in the observation's final segment but that have events in any
other included segment will be zeroed.
"""

settings() = @add_arg_table! ArgParseSettings(
    prog = "export",
    description = DESCRIPTION,
    epilog = EPILOG,
    autofix_names = true,
    exit_after_help = false,
) begin
    "source"
        required = true
        help = "Set the path to the results location to operate on."

    "sink"
        help = """
            Define where to place the exported results. \
            If not provided, write to standard output.\
        """

    "--format", "-x"
        arg_type = Symbol
        help = """
            Choose the output format: one of `csv`, `tsv`, `arrow` or `h5ad`. \
            If missing, it will be determined from the SINK filename \
            extension, falling back to `csv` if SINK is also missing.\
        """

    "--path", "-p"
        arg_type = String
        action = :append_arg
        help = """
            Filter for simulation segments that have a path starting with this \
            prefix.\
        """

    "--model", "-m"
        arg_type = String
        action = :append_arg
        help = """
            Filter for simulation segments that have a model path starting \
            with this prefix. This option can be repeated to include \
            additional segments.\
        """

    "--labeled", "-l"
        arg_type = Regex
        action = :append_arg
        help = """
            Filter for simulation segments that have a label matching this \
            regular expression anywhere. This option can be repeated to \
            include additional segments.\
        """

    "--relabel-with-segment-id", "-s"
        action = :store_true
        help = """
            Automatically append each segment's id to its label. \
            (This prevents non-final segments from being filtered out.)\
        """

    "--relabel-with-timestamp", "-t"
        action = :store_true
        help = """
            Automatically append each segment's end timestamp to its label \
            *after* non-final segments have been filtered out.\
        """

    "--events", "-e"
        arg_type = Regex
        action = :append_arg
        help = """
            Filter for simulation segments that have a an events stream name \
            matching this regular expression anywhere. This option can be \
            repeated to include additional event streams.\
        """

    "--kinds", "-k"
        arg_type = Regex
        action = :append_arg
        help = """
            Filter for trajectory dimensions that have names matching this \
            regular expression anywhere. This option can be repeated to \
            include additional dimensions. If such an expression contains a \
            single capture group, matching dimensions will be renamed to the \
            captured match. If such a capture group is named, the matching \
            dimensions will be collected into a thus-named layer; matches are \
            attempted in order, and the first match determines which layer the \
            dimension will be attached to. If the output format supports it, \
            layers are directly represented in the output. (Currently, this only \
            applies only to the `h5ad` format.) The set of (renamed) \
            dimensions must be homogenous between all layers, or an error will \
            be raised.\
        """

    "--dry"
        action = :store_true
        help = "Don't export, only print matching simulation segments."
end

function run(arguments = ARGS)
    parsed = @something(
        parse_args(arguments, settings(), as_symbols = true),
        return 1
    )

    @eval import ExportTool
    Tool = @invokelatest (@__MODULE__).ExportTool
    try @invokelatest Tool.main(; parsed...)
    catch e
        if e == :help
            ArgParse.show_help(settings(); exit_when_done = false)
            return 1
        else
            rethrow()
        end
    end
end

end
