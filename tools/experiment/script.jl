module Script

using ArgParse

DESCRIPTION = """
Run a simulation experiment according to the provided specifications and collect
the results into files at LOCATION.
"""

EPILOG = """
The specifications need to be defined according the JSON-based schedule language
as described in the documentation. They will be copied to LOCATION and loaded
from a new top-level schedule named `experiment.schedule.json` that will also
include the package and Julia versions used to create the results.
Besides the specifications, results will include an index file `index.arrow`
recording the contained simulation segments as well as one or more Arrow events
streams containing the actual trajectory data in long format.



For more details, refer to the documentation.
"""

settings() = @add_arg_table! ArgParseSettings(
    prog = "experiment",
    description = DESCRIPTION,
    epilog = EPILOG,
    exit_after_help = false
) begin
    "--location", "-l"
        default = "results/{TIMESTAMP}/"
        help = """
            Set the path to the results location to operate on. \
            It specifies a directory if it ends on a `/` and a common filename \
            prefix otherwise. \
            The placeholder `{TIMESTAMP}` will be replaced by a current ISO \
            timestamp. \
            See more info below.\
        """

    "--prepare"
        action = :store_true
        help = "Only prepare the LOCATION, but don't run simulations."

    "--seed"
        default = "seed"
        help = "Set the random seed to define in the root specification."

    "--simulate"
        action = :store_true
        help = "Run simulations for a previously prepared LOCATION."

    "--progress"
        default = :bars
        arg_type = Symbol
        help = """
            Specify how to report progress: \
            with `--progress=bars` show progress bars; \
            with `--progress=simple` just print messages; \
            with `--progress=timing` write JSON with timing information; \
            with `--progress=none` don't show progress.\
        """

    "--dry"
        action = :store_true
        help = "Don't simulate, only print planned simulation segments."

    "specifications"
        nargs = '*'
        arg_type = String
        help = """
            Define the (JSON) schedule specification file(s) to be run. \
            See more info below.\
        """
end

function run(arguments = ARGS)
    parsed = @something(
        parse_args(arguments, settings(), as_symbols = true),
        return 1
    )

    if (
        !parsed[:prepare] &&
        !parsed[:simulate] &&
        isempty(parsed[:specifications])
    )
        ArgParse.show_help(settings(); exit_when_done = false)
        return 0
    end

    @eval import ExperimentTool
    Tool = @invokelatest (@__MODULE__).ExperimentTool
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
