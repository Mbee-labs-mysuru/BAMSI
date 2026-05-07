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

# ---- BAMSI unit tests ----

add_executable(bamsi_unit_tests
    unit/version_test.cpp
)

target_include_directories(bamsi_unit_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src
)

target_link_libraries(bamsi_unit_tests PRIVATE
    bamsi-core
    gtest_main
)

add_test(NAME bamsi_unit_tests
    COMMAND $<TARGET_FILE:bamsi_unit_tests>
)

set_tests_properties(bamsi_unit_tests PROPERTIES
    LABELS "unit"
)

add_executable(bamsi_integration_tests
    # integration/clismoketest.cpp
    integration/version_cli_test.cpp
    # other integration tests...
)

target_link_libraries(bamsi_integration_tests PRIVATE
    bamsi-core
    gtest_main
)

add_test(
    NAME bamsi_integration_tests
    COMMAND $<TARGET_FILE:bamsi_integration_tests>
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
)

# ---- TIER 2 Integration Tests (Contract §S4) ----

add_executable(bamsi_tier2_tests
    test_tier2_integration.cpp
)

target_include_directories(bamsi_tier2_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(bamsi_tier2_tests PRIVATE
    bamsi-core
    gtest_main
)

add_test(
    NAME bamsi_tier2_tests
    COMMAND bamsi_tier2_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(bamsi_tier2_tests PROPERTIES
    LABELS "tier2;integration;codec;error-handling"
    TIMEOUT 120
)

# ---- Error Handling Sweep Tests (Contract §0.8) ----

add_executable(bamsi_error_sweep_tests
    test_error_sweep.cpp
)

target_include_directories(bamsi_error_sweep_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(bamsi_error_sweep_tests PRIVATE
    bamsi-core
    gtest_main
)

add_test(
    NAME bamsi_error_sweep_tests
    COMMAND bamsi_error_sweep_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(bamsi_error_sweep_tests PROPERTIES
    LABELS "error-handling;stage4"
    TIMEOUT 120
)

# ---- Findings & Fixes Regression Tests (Execution Plan §4.3.3) ----

add_executable(bamsi_ff_regression_tests
    test_findings_fixes_regression.cpp
)

target_include_directories(bamsi_ff_regression_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(bamsi_ff_regression_tests PRIVATE
    bamsi-core
    gtest_main
)

add_test(
    NAME bamsi_ff_regression_tests
    COMMAND bamsi_ff_regression_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(bamsi_ff_regression_tests PROPERTIES
    LABELS "stage4;regression;findings-fixes"
    TIMEOUT 180
)

# ---- V5 Feature Tests (ISA, SARange, Bidirectional) ----

add_executable(bamsi_v5_tests
    test_v5_features.cpp
)

target_include_directories(bamsi_v5_tests PRIVATE
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(bamsi_v5_tests PRIVATE
    bamsi-core
    gtest_main
)

add_test(
    NAME bamsi_v5_tests
    COMMAND bamsi_v5_tests
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)

set_tests_properties(bamsi_v5_tests PROPERTIES
    LABELS "v5;isa;sarange;enhanced"
    TIMEOUT 120
)
