cmake_minimum_required(VERSION 3.20)
project(gpec_cartesian_query_tests LANGUAGES C Fortran)
set(CMAKE_BUILD_TYPE Release CACHE STRING "Build configuration")
get_filename_component(GPEC_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/.." ABSOLUTE)
if(NOT CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  message(FATAL_ERROR "This standalone native test adapter supports GNU Fortran")
endif()
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
find_path(NETCDF_INCLUDE netcdf.mod REQUIRED)
find_library(NETCDFF_LIBRARY netcdff REQUIRED)
find_library(NETCDF_LIBRARY netcdf REQUIRED)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/version.inc"
  "      CHARACTER(len=*), PARAMETER :: version='cartesian-query-tests'\n")

function(gpec_library name)
  file(GLOB src "${GPEC_SOURCE}/${name}/*.f" "${GPEC_SOURCE}/${name}/*.F"
    "${GPEC_SOURCE}/${name}/*.f90" "${GPEC_SOURCE}/${name}/*.F90"
    "${GPEC_SOURCE}/${name}/*.c")
  list(FILTER src EXCLUDE REGEX "/(pentrc\\.F90|slayer\\.f)$")
  add_library(g_${name} STATIC ${src})
  set_target_properties(g_${name} PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/${name}")
  target_include_directories(g_${name} PUBLIC
    "${CMAKE_CURRENT_BINARY_DIR}/modules/${name}")
  target_include_directories(g_${name} PRIVATE "${NETCDF_INCLUDE}"
    "${CMAKE_CURRENT_BINARY_DIR}")
  target_compile_options(g_${name} PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:-ffixed-line-length-none;-fallow-argument-mismatch;-frecursive;-fopenmp>)
  target_link_libraries(g_${name} PUBLIC ${ARGN} "${NETCDFF_LIBRARY}"
    "${NETCDF_LIBRARY}" ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  target_link_options(g_${name} PUBLIC -fopenmp)
endfunction()
gpec_library(equil)
gpec_library(lsode g_equil)
gpec_library(vacuum g_equil g_lsode)
gpec_library(pentrc g_equil g_lsode)
gpec_library(coil g_equil)
gpec_library(slayer g_equil g_lsode)
foreach(code dcon gpec)
  file(GLOB src "${GPEC_SOURCE}/${code}/*.f" "${GPEC_SOURCE}/${code}/*.F")
  add_executable(${code} ${src})
  set_target_properties(${code} PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules/${code}")
  target_include_directories(${code} PRIVATE "${NETCDF_INCLUDE}"
    "${CMAKE_CURRENT_BINARY_DIR}" "${CMAKE_CURRENT_BINARY_DIR}/modules/${code}")
  target_compile_options(${code} PRIVATE -ffixed-line-length-none
    -fallow-argument-mismatch -frecursive -fopenmp)
  target_link_libraries(${code} PRIVATE g_pentrc g_vacuum g_slayer g_coil)
endforeach()

set(QUERY_CASE "${CMAKE_CURRENT_BINARY_DIR}/query_case")
file(MAKE_DIRECTORY "${QUERY_CASE}")
configure_file(fixtures/vac.in "${QUERY_CASE}/vac.in" COPYONLY)
execute_process(COMMAND "${Python3_EXECUTABLE}"
  "${CMAKE_CURRENT_SOURCE_DIR}/generate_query_fixture.py" "${QUERY_CASE}"
  RESULT_VARIABLE FIXTURE_STATUS)
if(NOT FIXTURE_STATUS EQUAL 0)
  message(FATAL_ERROR "Synthetic native query fixture generation failed")
endif()
enable_testing()
add_executable(test_vacuum_query test_vacuum_query.f90)
target_link_libraries(test_vacuum_query PRIVATE g_vacuum)
target_compile_options(test_vacuum_query PRIVATE -ffree-line-length-none -fcheck=all)
add_test(NAME vacuum_query_exact COMMAND test_vacuum_query)
set_tests_properties(vacuum_query_exact PROPERTIES
  WORKING_DIRECTORY "${QUERY_CASE}" TIMEOUT 90)
add_executable(test_cartesian_postprocessors test_cartesian_postprocessors.f90)
target_link_libraries(test_cartesian_postprocessors PRIVATE g_vacuum)
add_test(NAME cartesian_postprocessors_finite_and_raw COMMAND test_cartesian_postprocessors)
foreach(mode chebyshev divzero)
  add_test(NAME cartesian_postprocessors_reject_${mode}
    COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/expect_rejection.py"
      "$<TARGET_FILE:test_cartesian_postprocessors>" "${mode}")
endforeach()
set_tests_properties(cartesian_postprocessors_finite_and_raw
  cartesian_postprocessors_reject_chebyshev cartesian_postprocessors_reject_divzero
  PROPERTIES TIMEOUT 15)
