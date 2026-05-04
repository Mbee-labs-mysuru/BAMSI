add_library(bamsi-core STATIC
    cabi/version.cpp
    cli/app.cpp
    cli/dispatch.cpp
    cli/help.cpp
    cli/version.cpp
    format/placeholder.cpp
)

target_include_directories(bamsi-core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/generated/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(bamsi-core
    PRIVATE
        sdsl
        Threads::Threads
        z
)

target_compile_features(bamsi-core PUBLIC cxx_std_20)
