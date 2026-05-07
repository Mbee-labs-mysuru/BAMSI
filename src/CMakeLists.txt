add_library(bamsi-core STATIC
    # CLI / App
    cabi/version.cpp
    cli/app.cpp
    cli/build.cpp
    cli/dispatch.cpp
    cli/help.cpp
    cli/version.cpp

    # Build Pipeline Stages
    ingest/ingest.cpp
    ordering/ordering.cpp
    seqbuilder/seqbuilder.cpp
    sais/sais.cpp
    fmindex/fmindex.cpp
    streamencode/streamencode.cpp
    streamencode/mtf_rle_arith.cpp
    windows/windows.cpp
    bitvectors/bitvectors.cpp
    mapping/mapping.cpp
    query/query.cpp
    seal/seal.cpp
    sarange/sarange.cpp

    # Format (.bsi reader)
    format/format.cpp

    # C ABI
    cabi/cabi.cpp
)

target_include_directories(bamsi-core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/generated/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/external/htslib
        ${CMAKE_SOURCE_DIR}/external/deps/include
)

target_link_libraries(bamsi-core
    PUBLIC
        sdsl
    PRIVATE
        libsais
        libzstd_static
        xxhash
        htslib_imported
        OpenSSL::Crypto
)

target_compile_features(bamsi-core PUBLIC cxx_std_20)
