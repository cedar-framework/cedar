find_path(Argobots_INCLUDE_DIR NAMES abt.h DOC "Argobots include directory")

find_library(Argobots_LIBRARY NAMES abt DOC "Argobots library")

find_package_handle_standard_args(Argobots
  REQUIRED_VARS Argobots_INCLUDE_DIR Argobots_LIBRARY)

if (Argobots_FOUND)
  set(Argobots_LIBRARIES ${Argobots_LIBRARY})
  set(Argobots_INCLUDE_DIRS ${Argobots_INCLUDE_DIR})
endif()

mark_as_advanced(Argobots_INCLUDE_DIR Argobots_LIBRARY)
