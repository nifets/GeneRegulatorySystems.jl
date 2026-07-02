module Script

using ArgParse

DESCRIPTION = """
An experimental GUI application to explore the results produced by the
experiment tool.
"""

settings() = @add_arg_table! ArgParseSettings(
    prog = "inspect",
    description = DESCRIPTION,
    autofix_names = true,
    exit_after_help = false,
) begin
    "location"
        required = true
        help = "Set the path to the results location to operate on."

    "--channel", "-c"
        help = "Preset the text field to filter segments by events stream."

    "--items-prefix", "-i"
        help = "Preset the text field to filter segments by path prefix."

    "--label-pattern", "-l"
        help = "Preset the text field to filter segments by label."

    "--displays", "-d"
        default = "selector,trajectory,model,legend,info"
        help = "Specify which plots to show."

    "--kinds", "-k"
        default = "activity,mrnas,proteins"
        help = """
            Specify which trajectory components to show. These may include the
            categories `activity`, `elongations`, `premrnas`, `mrnas` and
            `proteins` as well as individual species names defined as part of
            additional mass-action reactions in the simulation schedule.\
        """

    "--size", "-s"
        default = "1280x720"
        help = "The initial size (in px) of the window."

    "--ppi", "-r"
        default = 96.0
        arg_type = Float64
        help = "The resolution (per 96px) to use when saving the figure."

    "--out", "-o"
        help = "A file to save to, instead of displaying the window."

    "--no-wait-for-close"
        action = :store_false
        dest_name = "wait_for_close"
        help = """
            Close the window immediately after loading. \
            This option is used during precompilation.\
        """
end

function run(arguments = ARGS)
    parsed = @something(
        parse_args(arguments, settings(), as_symbols = true),
        return 1
    )

    @eval import InspectTool
    Tool = @invokelatest (@__MODULE__).InspectTool
    @invokelatest Tool.main(; parsed...)
end

end
