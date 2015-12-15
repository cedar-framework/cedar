project(boxmg)

enable_language(C CXX Fortran)

cinch_add_application_directory(app)
cinch_add_library_target(boxmg src)

set(CINCH_HEADER_SUFFIXES "\\.h")

add_definitions(-DMSG_DOUBLE -DNO_SECOND_UNDERSCORE -fPIC -DBOUNDS_CHECK)
include_directories(include)
install(DIRECTORY include/ DESTINATION include)
