cmake_minimum_required(VERSION 3.20)

project(bamsi VERSION 0.1.0 LANGUAGES CXX)

set(BAMSI_VERSION "0.1.0")
set(BAMSI_FORMAT_VERSION 1)

configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/include/bamsi/config.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/generated/include/bamsi/config.hpp
    @ONLY
)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(GNUInstallDirs)
find_package(Threads REQUIRED)

add_subdirectory(external/sdsl-lite)
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
