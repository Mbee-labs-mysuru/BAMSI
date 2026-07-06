add_library(bamsix-core STATIC
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
    reconstruct/reconstruct.cpp

    # Format (.bsi reader)
    format/format.cpp

    # C ABI
    cabi/cabi.cpp
)

target_include_directories(bamsix-core
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/generated/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/external/htslib
        ${CMAKE_SOURCE_DIR}/external/deps/include
)

target_link_libraries(bamsix-core
    PUBLIC
        sdsl
        OpenSSL::Crypto
    PRIVATE
        libsais
        libzstd_static
        xxhash
        htslib_imported
)

target_compile_features(bamsix-core PUBLIC cxx_std_20)
