cmake_minimum_required(VERSION 3.16)

project(RobustCause VERSION 0.1.0 LANGUAGES CXX)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

option(ROBUSTCAUSE_BUILD_EXAMPLES "Build example programs" ON)
option(ROBUSTCAUSE_BUILD_TESTS "Build smoke tests" ON)
option(ROBUSTCAUSE_BUILD_CAPI "Build the C ABI layer for language bindings" ON)

find_package(Armadillo REQUIRED)

set(ROBUSTCAUSE_SOURCES
  src/robust_adstock.cpp
  src/robust_mmm.cpp
  src/robust_rlm_hc.cpp
  src/robust_sc.cpp
  src/robust_s_estimator.cpp
)

if(ROBUSTCAUSE_BUILD_CAPI)
  list(APPEND ROBUSTCAUSE_SOURCES src/capi.cpp)
endif()

add_library(robustcause ${ROBUSTCAUSE_SOURCES})
add_library(RobustCause::robustcause ALIAS robustcause)

target_compile_features(robustcause PUBLIC cxx_std_17)
target_include_directories(robustcause
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(robustcause PUBLIC armadillo)
set_target_properties(robustcause PROPERTIES
  OUTPUT_NAME robustcause
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  EXPORT_NAME robustcause
)

if(ROBUSTCAUSE_BUILD_EXAMPLES)
  add_executable(example_s_formula examples/example_s_formula.cpp)
  target_link_libraries(example_s_formula PRIVATE RobustCause::robustcause)

  add_executable(example_rlm_hc examples/example_rlm_hc.cpp)
  target_link_libraries(example_rlm_hc PRIVATE RobustCause::robustcause)

  add_executable(example_adstock examples/example_adstock.cpp)
  target_link_libraries(example_adstock PRIVATE RobustCause::robustcause)

  add_executable(example_mmm examples/example_mmm.cpp)
  target_link_libraries(example_mmm PRIVATE RobustCause::robustcause)

  add_executable(example_sc examples/example_sc.cpp)
  target_link_libraries(example_sc PRIVATE RobustCause::robustcause)

  if(ROBUSTCAUSE_BUILD_CAPI)
    add_executable(example_capi examples/example_capi.cpp)
    target_link_libraries(example_capi PRIVATE RobustCause::robustcause)
  endif()
endif()

if(ROBUSTCAUSE_BUILD_TESTS)
  enable_testing()
  add_executable(robustcause_smoke tests/smoke.cpp)
  target_link_libraries(robustcause_smoke PRIVATE RobustCause::robustcause)
  add_test(NAME robustcause_smoke COMMAND robustcause_smoke)
endif()

install(
  TARGETS robustcause
  EXPORT RobustCauseTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

install(
  EXPORT RobustCauseTargets
  FILE RobustCauseTargets.cmake
  NAMESPACE RobustCause::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RobustCause
)

configure_package_config_file(
  cmake/RobustCauseConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/RobustCauseConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RobustCause
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/RobustCauseConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/RobustCauseConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/RobustCauseConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/RobustCause
)
