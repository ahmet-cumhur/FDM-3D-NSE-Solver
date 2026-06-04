cmake_minimum_required(VERSION 3.20)

project(3D_NSE_SOLVER LANGUAGES Fortran)

# =========================
# Options
# =========================
option(USE_IBM "Enable immersed boundary method (IBM)" OFF)
option(USE_IBM_G "Enable second-order immersed boundary method" OFF)
option(USE_OPENMP "Enable OpenMP parallelization" ON)
# =========================
# Sources
# =========================
set(SOURCES
    main.f90
    init.f90
    boundary.f90
    io.f90
    step.f90
    3d_fftw.f90
    thomas.f90
    rk3_step.f90
    rk3_rhs.f90
    richardson.f90
)

# Add IBM sources only if enabled
if(USE_IBM OR USE_IBM_G)
    list(APPEND SOURCES
        ibm.f90
    )
endif()

add_executable(main ${SOURCES})
# =========================
# OPENMP
# =========================
if(USE_OPENMP)
    find_package(OpenMP REQUIRED)

    target_link_libraries(main PRIVATE OpenMP::OpenMP_Fortran)
    target_compile_definitions(main PRIVATE USE_OPENMP)
endif()
# =========================
# Module directory
# =========================
set_target_properties(main PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod
)

# =========================
# Compiler flags
# =========================
target_compile_options(main PRIVATE -O3 -march=native)

# Enable preprocessing depending on compiler
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(main PRIVATE -cpp -ffree-line-length-none)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    target_compile_options(main PRIVATE -fpp)
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    target_compile_options(main PRIVATE -Mpreprocess)
endif()

# =========================
# Define macros
# =========================
if(USE_IBM)
    target_compile_definitions(main PRIVATE USE_IBM)
endif()
if(USE_IBM_G)
    target_compile_definitions(main PRIVATE USE_IBM_G)
endif()

# =========================
# FFTW
# =========================

# Optional user-provided FFTW root directory, e.g.
# cmake -B build -DFFTW3_ROOT=/path/to/fftw
set(FFTW3_ROOT "" CACHE PATH "Root directory of FFTW3 installation")

find_path(FFTW3_INCLUDE_DIR
    NAMES fftw3.f03
    HINTS
        ${FFTW3_ROOT}/include
        ${FFTW3_ROOT}
    PATH_SUFFIXES include
)

find_library(FFTW3_LIB
    NAMES fftw3
    HINTS
        ${FFTW3_ROOT}/lib
        ${FFTW3_ROOT}/lib64
        ${FFTW3_ROOT}
    PATH_SUFFIXES lib lib64
)

if(NOT FFTW3_INCLUDE_DIR)
    message(FATAL_ERROR "Could not find fftw3.f03. Set -DFFTW3_ROOT=/path/to/fftw or -DFFTW3_INCLUDE_DIR=/path/to/include")
endif()

if(NOT FFTW3_LIB)
    message(FATAL_ERROR "Could not find FFTW3 library. Set -DFFTW3_ROOT=/path/to/fftw or -DFFTW3_LIB=/path/to/libfftw3.so")
endif()

message(STATUS "FFTW3 include dir: ${FFTW3_INCLUDE_DIR}")
message(STATUS "FFTW3 library: ${FFTW3_LIB}")

target_include_directories(main PRIVATE ${FFTW3_INCLUDE_DIR})
target_link_libraries(main PRIVATE ${FFTW3_LIB})


# =========================
# Optional: include dirs (only if really needed)
# =========================
# target_include_directories(main PRIVATE /usr/include)
