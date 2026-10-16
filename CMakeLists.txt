# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>

cmake_minimum_required(VERSION 3.16)
project(spclust VERSION 0.3.3 LANGUAGES CXX C
        DESCRIPTION "Agglomerative hierarchical clustering over a sparse similarity matrix")

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------
# Note: the linkage is a **run-time** argument (spclust::Linkage at construction), so the
# library has no compile-time options - one binary covers all three. The CLI uses --linkage.
# Defaults are ON only for a standalone build of this repository, OFF under add_subdirectory.
# option() is used deliberately: it respects the cached value, so -DSPCLUST_BUILD_*=OFF works.
# (This used to be set(... ON), which silently overwrote an explicit OFF.)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(_spclust_default ON)
else()
    set(_spclust_default OFF)
endif()
option(SPCLUST_BUILD_TOOLS "build the spclust_cli command-line tool" ${_spclust_default})
option(SPCLUST_BUILD_TESTS "build the tests and register them with ctest" ${_spclust_default})

# ---------------------------------------------------------------------------
# The library (no dependency on tools/)
# ---------------------------------------------------------------------------
set(SPCLUST_LIB_SOURCES src/hac.cpp src/scipy.cpp src/spclust_c.cpp)
file(GLOB SPCLUST_PUBLIC_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/spclust/*.h)
add_library(spclust ${SPCLUST_LIB_SOURCES} ${SPCLUST_PUBLIC_HEADERS})
add_library(spclust::spclust ALIAS spclust)

target_compile_features(spclust PUBLIC cxx_std_17)
target_include_directories(spclust
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
if(MSVC)
    target_compile_options(spclust PRIVATE /W4)
else()
    target_compile_options(spclust PRIVATE -Wall -Wextra)
endif()

# ---------------------------------------------------------------------------
# The command-line tool (the library's first user, owning all IO and presentation)
# ---------------------------------------------------------------------------
if(SPCLUST_BUILD_TOOLS)
    add_executable(spclust_cli tools/spclust_cli.cpp tools/readdb.cpp)
    target_link_libraries(spclust_cli PRIVATE spclust::spclust)
    # The CLI's own compile-time options (not part of the library): default cut, debug level, input bound
    target_include_directories(spclust_cli PRIVATE tools)
    target_compile_definitions(spclust_cli PRIVATE MIN_SIM_CLUST=50 DEBUG=1)
    if(NOT MSVC)
        target_compile_options(spclust_cli PRIVATE -Wall -Wextra)
    endif()
endif()

# ---------------------------------------------------------------------------
# Install / find_package
# ---------------------------------------------------------------------------
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(TARGETS spclust EXPORT spclustTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
if(SPCLUST_BUILD_TOOLS)
    install(TARGETS spclust_cli RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()

install(EXPORT spclustTargets
        FILE spclustTargets.cmake
        NAMESPACE spclust::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/spclust)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/spclustConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/spclustConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/spclustConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/spclust)
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/spclustConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/spclustConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/spclust)

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
if(SPCLUST_BUILD_TESTS)
    enable_testing()
    find_package(Python3 COMPONENTS Interpreter REQUIRED)

    # Library-level property test: the three cut entry points (expectations from leaf sets, bottom-up)
    if(TARGET spclust_cli)
        set(cli_bin ${CMAKE_CURRENT_BINARY_DIR}/spclust_cli)
    else()
        # With no CLI target there is nothing for the CLI-level tests to run
        set(cli_bin "")
    endif()

    add_executable(cut_property tests/test_cut.cpp)
    target_link_libraries(cut_property PRIVATE spclust::spclust)
    add_test(NAME cut_property COMMAND cut_property)

    # C API: this test is compiled by a **C compiler**, so it also proves the header is plain C.
    # CMake uses the C compiler for .c and links with C++ (the library is C++).
    add_executable(c_api tests/test_c_api.c)
    target_link_libraries(c_api PRIVATE spclust::spclust)
    set_property(TARGET c_api PROPERTY C_STANDARD 11)
    add_test(NAME c_api COMMAND c_api)

    # FFI: build the same library again as a shared object and call it from Python via ctypes
    add_library(spclust_shared SHARED ${SPCLUST_LIB_SOURCES} ${SPCLUST_PUBLIC_HEADERS})
    target_include_directories(spclust_shared PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
    target_compile_features(spclust_shared PUBLIC cxx_std_17)
    # The Python binding on top of the same shared library.
    add_test(NAME python_binding
             COMMAND ${CMAKE_COMMAND} -E env PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python
                     ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_python.py
                     --lib $<TARGET_FILE:spclust_shared>
                     --cli $<TARGET_FILE:spclust_cli>)
    set_tests_properties(python_binding PROPERTIES
                         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    add_test(NAME c_api_ffi
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_c_api.py
                     --lib $<TARGET_FILE:spclust_shared>)
    set_tests_properties(c_api_ffi PROPERTIES
                         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    # Differential test: an independent oracle, merge by merge, plus the CLI contract tests.
    # The linkage is a run-time argument now, so one binary with --linkage 0/1/2 covers all three.
    if(cli_bin)
        add_test(NAME differential
                 COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py
                         --bin ${cli_bin} --min-sim 50)
        add_test(NAME cli
                 COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cli.py
                         --bin ${cli_bin})
        # scipy interop: hand the --scipy export to the real scipy. Skips when numpy/scipy are absent.
        add_test(NAME scipy_interop
                 COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_scipy.py
                         --bin ${cli_bin})
        set_tests_properties(differential cli scipy_interop PROPERTIES
                             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    endif()
endif()
