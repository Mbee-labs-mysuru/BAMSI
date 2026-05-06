add_test(NAME bamsi_cli_smoke
    COMMAND $<TARGET_FILE:bamsi>
)
set_tests_properties(bamsi_cli_smoke PROPERTIES
    PASS_REGULAR_EXPRESSION "usage: bamsi <subcommand>"
    LABELS "cli;smoke"
)

add_test(NAME bamsi_version_smoke
    COMMAND $<TARGET_FILE:bamsi> version
)
set_tests_properties(bamsi_version_smoke PROPERTIES
    PASS_REGULAR_EXPRESSION "0\\.1\\.0"
    LABELS "cli;version;smoke"
)

add_test(NAME bamsi_invalid_subcommand
    COMMAND $<TARGET_FILE:bamsi> definitely-invalid
)
set_tests_properties(bamsi_invalid_subcommand PROPERTIES
    WILL_FAIL TRUE
    LABELS "cli;error-path"
)

add_test(NAME bamsi_help_smoke
    COMMAND $<TARGET_FILE:bamsi> --help
)
set_tests_properties(bamsi_help_smoke PROPERTIES
    PASS_REGULAR_EXPRESSION "subcommands: version, build, count, exists, locate, region-count, region-exists, reconstruct, info, verify"
    LABELS "cli;help;smoke"
)

add_test(NAME bamsi_short_help_smoke
    COMMAND $<TARGET_FILE:bamsi> -h
)
set_tests_properties(bamsi_short_help_smoke PROPERTIES
    PASS_REGULAR_EXPRESSION "subcommands: version, build, count, exists, locate, region-count, region-exists, reconstruct, info, verify"
    LABELS "cli;help;smoke"
)

add_test(NAME bamsi_help_subcommand_smoke
    COMMAND $<TARGET_FILE:bamsi> help
)
set_tests_properties(bamsi_help_subcommand_smoke PROPERTIES
    PASS_REGULAR_EXPRESSION "subcommands: version, build, count, exists, locate, region-count, region-exists, reconstruct, info, verify"
    LABELS "cli;help;smoke"
)

add_test(NAME bamsi_multi_arg_error
    COMMAND $<TARGET_FILE:bamsi> foo bar
)
set_tests_properties(bamsi_multi_arg_error PROPERTIES
    WILL_FAIL TRUE
    LABELS "cli;error-path"
)
