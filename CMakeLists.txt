# set appropriate cmake version to enable CUDA
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

# set path to CUDA compiler
set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")
project(keccaksolver CXX CUDA)

# set a default build type if none was specified
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message("Setting build type to 'Release' as none was specified.")
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif ()

# ========================================================================== #
# language standard requirements
# ========================================================================== #

if (NOT DEFINED CMAKE_C_STANDARD)
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED true)
endif ()

if (NOT DEFINED CMAKE_CUDA_STANDARD)
    set(CMAKE_CUDA_STANDARD 11)
    set(CMAKE_CUDA_STANDARD_REQUIRED true)
endif ()

# ========================================================================== #
# compilation flags
# ========================================================================== #

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(CheckCXXCompilerFlag)

#find_package(OpenMP)
#if (OPENMP_FOUND)
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -pthread")
#endif ()

check_cxx_compiler_flag("-Wall" COMPILER_CXX_ALL_WARN)
if (COMPILER_CXX_ALL_WARN)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall"
            CACHE STRING "Show all warning messages" FORCE)
endif ()

check_cxx_compiler_flag("-Wextra" COMPILER_CXX_EXTRA_WARN)
if (COMPILER_CXX_EXTRA_WARN)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra"
            CACHE STRING "Show extra warning messages" FORCE)
endif ()

check_cxx_compiler_flag("-pipe" COMPILER_CXX_PIPE)
if (COMPILER_CXX_PIPE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pipe"
            CACHE STRING "Speed up compilation by piping" FORCE)
endif ()

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler '${CMAKE_C_FLAGS}'"
        CACHE STRING "Pass the same C flags to back-end compiler of nvcc" FORCE)

check_cxx_compiler_flag("-mavx" COMPILER_CXX_AVX)
if (COMPILER_CXX_AVX)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mavx"
            CACHE STRING "Optimize with AVX instructions" FORCE)
endif ()

check_cxx_compiler_flag("-mavx2" COMPILER_CXX_AVX2)
if (COMPILER_CXX_AVX2)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mavx2"
            CACHE STRING "Optimize with AVX2 instructions" FORCE)
endif ()

check_cxx_compiler_flag("-march=native" COMPILER_CXX_ARCH)
if (COMPILER_CXX_ARCH)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native"
            CACHE STRING "Optimize based on the architecture" FORCE)
endif ()

check_cxx_compiler_flag("-mtune=native" COMPILER_CXX_TUNE)
if (COMPILER_CXX_TUNE)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -mtune=native"
            CACHE STRING "Fine-tune the program based on the CPU" FORCE)
endif ()

check_cxx_compiler_flag("-fomit-frame-pointer" COMPILER_CXX_NO_FP)
if (COMPILER_CXX_NO_FP)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -fomit-frame-pointer"
            CACHE STRING "Omit frame pointers" FORCE)
endif ()

check_cxx_compiler_flag("-O3" COMPILER_CXX_OPTIMIZATION)
if (COMPILER_CXX_OPTIMIZATION)
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3"
            CACHE STRING "Trade floating point precision for speed" FORCE)
endif ()

set(CMAKE_CUDA_FLAGS_RELEASE
        "${CMAKE_CUDA_FLAGS_RELEASE} -Xcompiler '${CMAKE_CXX_FLAGS_RELEASE}'"
        CACHE STRING "Pass the same CXX release flags to back-end compiler of nvcc"
        FORCE)


# ========================================================================== #
# dat file
# ========================================================================== #
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/data/append_analysis.dat ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/data/mq_analysis.dat ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/data/lin_analysis.dat ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/data/variate_analysis.dat ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

execute_process(
        COMMAND bash -c "mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/testBench/ \
        && cp ${CMAKE_CURRENT_SOURCE_DIR}/test/testBench/testMessage.txt ${CMAKE_CURRENT_BINARY_DIR}/testBench/ \
        && cp ${CMAKE_CURRENT_SOURCE_DIR}/test/testBench/test_append_analysis.dat ${CMAKE_CURRENT_BINARY_DIR}/testBench/ \
        && cp ${CMAKE_CURRENT_SOURCE_DIR}/test/testBench/test_mq_analysis.dat ${CMAKE_CURRENT_BINARY_DIR}/testBench/ \
        && cp ${CMAKE_CURRENT_SOURCE_DIR}/test/testBench/test_variate_analysis.dat ${CMAKE_CURRENT_BINARY_DIR}/testBench/ \
        && cp ${CMAKE_CURRENT_SOURCE_DIR}/test/testBench/test_lin_analysis.dat ${CMAKE_CURRENT_BINARY_DIR}/testBench/"
        OUTPUT_VARIABLE outVar
)

# ========================================================================== #
# source code
# ========================================================================== #

set(CXX_SRCS
        src/options.cpp
        src/utils.cpp
        src/threadpool.cpp
        src/kmath.cpp
        src/keccak.cpp)

set(CUDA_SRCS
        src/keccak.cu
        src/solver.cu
        src/cuda_utils.cu)


# ========================================================================== #
# main program
# ========================================================================== #

add_executable(${PROJECT_NAME} ${CXX_SRCS} ${CUDA_SRCS} main.cpp)

target_link_libraries(${PROJECT_NAME}
        PRIVATE m
        PRIVATE pthread)

target_include_directories(${PROJECT_NAME}
        PRIVATE ${PROJECT_SOURCE_DIR}/include
        PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

set_target_properties(${PROJECT_NAME} PROPERTIES
        LINKER_LANGUAGE CXX
        CUDA_SEPARABLE_COMPILATION ON)

add_executable(testRandMessage ${CXX_SRCS} ${CUDA_SRCS} test/testRandMessage.cpp)
target_include_directories(testRandMessage
        PRIVATE ${PROJECT_SOURCE_DIR}/include
        PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})

add_executable(testPreimageAttack ${CXX_SRCS} ${CUDA_SRCS} test/testPreimageAttack.cpp)
target_include_directories(testPreimageAttack
        PRIVATE ${PROJECT_SOURCE_DIR}/include
        PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
target_compile_definitions(testPreimageAttack
        PRIVATE TEST_PRE=1)

add_executable(testConstraintCorrectness test/testConstraintCorrectness.cpp)

#add_executable(testMQIO debug/testMQIO.cpp)
#add_executable(testAp debug/testAp.cpp)
#add_executable(testMQ debug/testMQ.cu src/keccak.cpp)
#target_include_directories(testMQ
#        PRIVATE ${PROJECT_SOURCE_DIR}/include
#        PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
