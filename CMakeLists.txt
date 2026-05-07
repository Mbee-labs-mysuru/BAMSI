cmake_minimum_required(VERSION 3.20)

project(bamsi VERSION 0.1.0 LANGUAGES C CXX)

set(BAMSI_VERSION "0.1.0")
set(BAMSI_FORMAT_VERSION 6)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/bamsi/config.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/generated/include/bamsi/config.hpp
    @ONLY
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_C_STANDARD 11)

include(GNUInstallDirs)
find_package(Threads REQUIRED)

# ─── OpenSSL for SHA-256 ──────────────────────────────────────────────────────
find_package(OpenSSL REQUIRED)

# ─── sdsl-lite ────────────────────────────────────────────────────────────────
add_subdirectory(external/sdsl-lite)
set_property(TARGET sdsl PROPERTY CXX_STANDARD 17)
set_property(TARGET sdsl PROPERTY CXX_STANDARD_REQUIRED ON)
set_property(TARGET sdsl PROPERTY CXX_EXTENSIONS OFF)

# ─── libsais ─────────────────────────────────────────────────────────────────
add_subdirectory(external/libsais)

# ─── zstd (only the library, no programs) ─────────────────────────────────────
set(ZSTD_BUILD_PROGRAMS OFF CACHE BOOL "Don't build zstd programs" FORCE)
set(ZSTD_BUILD_TESTS OFF CACHE BOOL "Don't build zstd tests" FORCE)
set(ZSTD_BUILD_SHARED OFF CACHE BOOL "Don't build shared lib" FORCE)
set(ZSTD_BUILD_STATIC ON CACHE BOOL "Build static lib" FORCE)
add_subdirectory(external/zstd/build/cmake)

# ─── xxHash (header-only with XXH_INLINE_ALL) ────────────────────────────────
add_library(xxhash INTERFACE)
target_include_directories(xxhash INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/external/xxHash)

# ─── Local deps (bzip2 built within project) ─────────────────────────────────
# bzip2 headers and static library for htslib CRAM support.
# Kept in external/deps/ to avoid system dependency.
set(BAMSI_DEPS_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/deps)

add_library(bzip2_local STATIC IMPORTED GLOBAL)
set_target_properties(bzip2_local PROPERTIES
    IMPORTED_LOCATION ${BAMSI_DEPS_DIR}/lib/libbz2.a
)
target_include_directories(bzip2_local INTERFACE ${BAMSI_DEPS_DIR}/include)

# ─── htslib (pre-built static library, curl disabled) ────────────────────────
# htslib rebuilt without libcurl for project isolation.
# Rebuild: cd external/htslib && make clean && make lib-static \
#   NONCONFIGURE_OBJS="" \
#   CFLAGS="-g -O2 -fPIC -I$(pwd)/../deps/include" -j4
add_library(htslib_imported STATIC IMPORTED GLOBAL)
set_target_properties(htslib_imported PROPERTIES
    IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/external/htslib/libhts.a
)
target_include_directories(htslib_imported INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/external/htslib
)
# htslib transitive dependencies (all from system — these are C runtime libs)
target_link_libraries(htslib_imported INTERFACE
    bzip2_local
    z
    lzma
    Threads::Threads
)

# ─── Main library and executable ─────────────────────────────────────────────
add_subdirectory(src)

add_executable(bamsi
    src/cli/main.cpp
)

target_link_libraries(bamsi
    PRIVATE
        bamsi-core
)

include(CTest)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()

install(TARGETS bamsi-core bamsi
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/generated/include/bamsi/config.hpp
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/bamsi
)
