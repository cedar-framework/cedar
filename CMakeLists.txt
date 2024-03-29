cmake_minimum_required(VERSION 3.12)
project(cedar VERSION 0.1.0 LANGUAGES CXX C Fortran)

set(CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

set(INCLUDE_INSTALL_DIR include/ CACHE STRING "")
set(CEDAR_TARGETS_EXPORT_NAME       "${PROJECT_NAME}Targets")
set(CEDAR_CONFIG_INSTALL_DIR        "lib/cmake/${PROJECT_NAME}")
set(CEDAR_CMAKE_CONFIG_TEMPLATE     "cmake/config.cmake.in")
set(CEDAR_CMAKE_CONFIG_DIR          "${CMAKE_CURRENT_BINARY_DIR}/cmake_config")
set(CEDAR_CMAKE_VERSION_CONFIG_FILE "${CEDAR_CMAKE_CONFIG_DIR}/${PROJECT_NAME}ConfigVersion.cmake")
set(CEDAR_CMAKE_PROJECT_CONFIG_FILE "${CEDAR_CMAKE_CONFIG_DIR}/${PROJECT_NAME}Config.cmake")

include(options)
include(testing)
include(ExternalProject)

include_directories(include)
install(DIRECTORY include/ DESTINATION include)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/include)
include_directories(${CMAKE_BINARY_DIR}/include)

ExternalProject_Add(
  tausch
  PREFIX ${CMAKE_BINARY_DIR}/tausch
  GIT_REPOSITORY https://github.com/luspi/tausch.git
  GIT_TAG 425f062eb5661ea16696b2ffc0813a7e6d4eb4e4
  UPDATE_DISCONNECTED 1
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  LOG_DOWNLOAD ON
)
ExternalProject_Get_Property(tausch source_dir)
set(tausch_location ${source_dir})

ExternalProject_Add(
  json
  PREFIX ${CMAKE_BINARY_DIR}/json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG v3.1.2
  UPDATE_DISCONNECTED 1
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  LOG_DOWNLOAD ON
)
ExternalProject_Get_Property(json source_dir)
set(json_location ${source_dir}/include)

FIND_PACKAGE(Boost 1.44 REQUIRED COMPONENTS system filesystem)

FIND_PACKAGE(BLAS REQUIRED)
FIND_PACKAGE(LAPACK REQUIRED)

find_package(MPI REQUIRED)

if (ENABLE_PLANE_AGG)
  find_package(Argobots REQUIRED)
endif()

add_subdirectory(src)
add_library(cedar ${cedar_src})
add_dependencies(cedar json tausch)
target_compile_definitions(cedar PUBLIC
  MSG_DOUBLE
  NO_SECOND_UNDERSCORE
  $<$<CONFIG:DEBUG>:DEBUG;BOUNDS_CHECK>)
target_compile_features(cedar PUBLIC cxx_std_11)
target_compile_options(cedar PUBLIC
  $<$<COMPILE_LANGUAGE:CXX>:-std=c++11>
  $<$<COMPILE_LANGUAGE:Fortran>:-fallow-argument-mismatch>
  $<$<CONFIG:DEBUG>:-O0;-Wall;-g>
  $<$<CONFIG:RELEASE>:-O3>)
set_target_properties(cedar PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_link_libraries(cedar PUBLIC
  ${Boost_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
  ${MPI_C_LIBRARIES}
  ${MPI_CXX_LIBRARIES}
  ${MPI_Fortran_LIBRARIES})
include_directories(${json_location})
include_directories(${tausch_location})
target_compile_definitions(cedar PUBLIC WITH_TAUSCH)
if (ENABLE_PLANE_AGG)
  target_link_libraries(cedar PUBLIC ${Argobots_LIBRARIES})
  target_include_directories(cedar PUBLIC ${Argobots_INCLUDE_DIRS})
  target_compile_definitions(cedar PUBLIC PLANE_AGG)
endif()
if (ENABLE_3D)
  target_compile_definitions(cedar PUBLIC ENABLE_3D)
endif()

target_include_directories(cedar PUBLIC
  $<INSTALL_INTERFACE:include/>
  ${MPI_C_INCLUDE_PATH}
  ${MPI_Fortran_INCLUDE_PATH}
  ${Boost_INCLUDE_DIR})

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${CEDAR_CMAKE_VERSION_CONFIG_FILE}
  COMPATIBILITY SameMajorVersion
)
configure_package_config_file(
  ${CEDAR_CMAKE_CONFIG_TEMPLATE}
  ${CEDAR_CMAKE_PROJECT_CONFIG_FILE}
  INSTALL_DESTINATION ${CEDAR_CONFIG_INSTALL_DIR}
  PATH_VARS INCLUDE_INSTALL_DIR
)

install(
  FILES ${CEDAR_CMAKE_PROJECT_CONFIG_FILE} ${CEDAR_CMAKE_VERSION_CONFIG_FILE}
  DESTINATION ${CEDAR_CONFIG_INSTALL_DIR}
  )

install(DIRECTORY ${json_location}/nlohmann DESTINATION "include")

install(TARGETS cedar DESTINATION "lib"
  EXPORT ${CEDAR_TARGETS_EXPORT_NAME}
  )

install(EXPORT ${CEDAR_TARGETS_EXPORT_NAME}
  DESTINATION ${CEDAR_CONFIG_INSTALL_DIR}
)


if(ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

if(ENABLE_UNIT_TESTS)
  find_package(PythonLibs REQUIRED)
  find_package(NumPy REQUIRED)
  message(STATUS ${PYTHON_INCLUDE_DIR})
  list(APPEND cedar-deps ${PYTHON_LIBRARY})
  include_directories(${PYTHON_INCLUDE_DIR} ${PYTHON_NUMPY_INCLUDE_DIR})

  add_subdirectory(test)
endif()

# add_subdirectory(app)
