cmake_minimum_required(VERSION 3.20)

project(3D_NSE_SOLVER LANGUAGES Fortran)

add_executable(main
    main.f90
    initi.f90
    step.f90
    boun.f90
    th_mod.f90
    3d_fftw.f90
    data_save.f90
    ibm_1st.f90
    step_ibm.f90
)

set_target_properties(main PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod
)

target_include_directories(main PRIVATE C:/msys64/ucrt64/include)
target_link_directories(main PRIVATE C:/msys64/ucrt64/lib)
target_link_libraries(main PRIVATE fftw3)
target_compile_options(main PRIVATE -O2)