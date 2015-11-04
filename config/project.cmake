project(boxmg)

enable_language(C CXX Fortran)

cinch_add_application_directory(app)
cinch_add_library_target(boxmg-common src/common)
cinch_add_library_target(boxmg-2d src/2d)

set(CINCH_HEADER_SUFFIXES "\\.h")

add_definitions(-DMSG_DOUBLE -DNO_SECOND_UNDERSCORE -fPIC)
