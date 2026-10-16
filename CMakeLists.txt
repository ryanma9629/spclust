# SPDX-License-Identifier: MIT
# Copyright (c) 2022-2026 Ryan Ma <pningma@sina.com>

cmake_minimum_required(VERSION 3.16)
project(spclust VERSION 0.3.1 LANGUAGES CXX
        DESCRIPTION "基于稀疏相似度矩阵的凝聚型层次聚类")

# ---------------------------------------------------------------------------
# 选项
# ---------------------------------------------------------------------------
# 注：链接方式是**运行时**的（spclust::Linkage，构造 HAC 时给），所以库本身没有
# 任何编译期选项 —— 一个二进制就能跑三种链接。命令行工具用 --linkage 选。
# 默认值只在"本仓库独立构建"时开；被 add_subdirectory 进来时默认关。
# 注意用 option() —— 它尊重用户已经在缓存的取值，所以 -DSPCLUST_BUILD_*=OFF 一定生效。
# （以前这里是 set(... ON)，会把用户显式传的 OFF 覆盖掉。）
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(_spclust_default ON)
else()
    set(_spclust_default OFF)
endif()
option(SPCLUST_BUILD_TOOLS "构建命令行工具 spclust_cli" ${_spclust_default})
option(SPCLUST_BUILD_TESTS "构建测试并注册 ctest" ${_spclust_default})

# ---------------------------------------------------------------------------
# 库（只有 src/hac.cpp；不依赖 tools/）
# ---------------------------------------------------------------------------
add_library(spclust src/hac.cpp src/scipy.cpp)
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
# 测试
# ---------------------------------------------------------------------------
if(SPCLUST_BUILD_TESTS)
    enable_testing()
    find_package(Python3 COMPONENTS Interpreter REQUIRED)

    # 库级属性测试：三个切树接口（期望值用"从合并树自底向上收集叶子"独立算出）
    if(TARGET spclust_cli)
        set(cli_bin ${CMAKE_CURRENT_BINARY_DIR}/spclust_cli)
    else()
        # 只编库时，测试需要一个带 CLI 的二进制
        set(cli_bin "")
    endif()

    add_executable(cut_property tests/test_cut.cpp)
    target_link_libraries(cut_property PRIVATE spclust::spclust)
    add_test(NAME cut_property COMMAND cut_property)

    # 差分测试：独立 oracle 逐次合并比对 + CLI 契约测试。
    # 链接方式现在是运行时的，所以一个二进制 --linkage 0/1/2 就覆盖三种。
    if(cli_bin)
        add_test(NAME differential
                 COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_tests.py
                         --bin ${cli_bin} --min-sim 50)
        add_test(NAME cli
                 COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cli.py
                         --bin ${cli_bin})
        # scipy 互操作：把 --scipy 的导出交给真的 scipy 读。没装 numpy/scipy 时自动跳过。
        add_test(NAME scipy_interop
                 COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_scipy.py
                         --bin ${cli_bin})
        set_tests_properties(differential cli scipy_interop PROPERTIES
                             WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    endif()
endif()
