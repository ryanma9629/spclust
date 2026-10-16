# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>

cmake_minimum_required(VERSION 3.16)
project(spclust VERSION 0.2.0 LANGUAGES CXX
        DESCRIPTION "基于稀疏相似度矩阵的凝聚型层次聚类")

# ---------------------------------------------------------------------------
# 选项
# ---------------------------------------------------------------------------
# 链接方式目前仍是编译期常量（换它要重新编译本库）。0=average 1=max(单链接) 2=min(nonzero)
set(SPCLUST_LINKAGE "0" CACHE STRING "linkage: 0=average 1=max 2=min(nonzero)")
set_property(CACHE SPCLUST_LINKAGE PROPERTY STRINGS 0 1 2)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(SPCLUST_BUILD_TOOLS ON)
    set(SPCLUST_BUILD_TESTS ON)
endif()
option(SPCLUST_BUILD_TOOLS "构建命令行工具 spclust_cli" ${SPCLUST_BUILD_TOOLS})
option(SPCLUST_BUILD_TESTS "构建测试用的三种 linkage 二进制并注册 ctest" ${SPCLUST_BUILD_TESTS})

# ---------------------------------------------------------------------------
# 库（只有 src/hac.cpp；不依赖 tools/）
# ---------------------------------------------------------------------------
add_library(spclust src/hac.cpp)
add_library(spclust::spclust ALIAS spclust)

target_compile_features(spclust PUBLIC cxx_std_17)
target_compile_definitions(spclust PRIVATE LINKAGE_MTD=${SPCLUST_LINKAGE})
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
# 命令行工具（库的第一个使用者，负责所有 IO 与展示）
# ---------------------------------------------------------------------------
if(SPCLUST_BUILD_TOOLS)
    add_executable(spclust_cli tools/spclust_cli.cpp tools/readdb.cpp)
    target_link_libraries(spclust_cli PRIVATE spclust::spclust)
    # CLI 自己的编译期选项（不属于库）：默认切树阈值、调试级别、输入上界
    target_include_directories(spclust_cli PRIVATE tools)
    target_compile_definitions(spclust_cli PRIVATE MIN_SIM_CLUST=50 DEBUG=1)
    if(NOT MSVC)
        target_compile_options(spclust_cli PRIVATE -Wall -Wextra)
    endif()
endif()

# ---------------------------------------------------------------------------
# 安装 / find_package
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
# 测试：三种 linkage 各一个二进制 + ctest
# ---------------------------------------------------------------------------
if(SPCLUST_BUILD_TESTS)
    enable_testing()
    find_package(Python3 COMPONENTS Interpreter REQUIRED)

    set(test_bins "")
    foreach(m 0 1 2)
        set(bin "spclust_cli_linkage_${m}")
        add_executable(${bin} tools/spclust_cli.cpp tools/readdb.cpp src/hac.cpp)
        target_compile_definitions(${bin} PRIVATE
            LINKAGE_MTD=${m} MIN_SIM_CLUST=50 DEBUG=1)
        target_include_directories(${bin} PRIVATE include tools)
        list(APPEND test_bins ${bin})
    endforeach()

    # 差分测试：独立 oracle 逐次合并比对 + CLI 契约测试
    add_test(NAME differential
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py
                     --bin-prefix ${CMAKE_CURRENT_BINARY_DIR}/spclust_cli_linkage_ --min-sim 50)
    add_test(NAME cli
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cli.py
                     --bin ${CMAKE_CURRENT_BINARY_DIR}/spclust_cli)
    set_tests_properties(differential cli PROPERTIES
                         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()
