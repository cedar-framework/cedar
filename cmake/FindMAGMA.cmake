if( NOT MAGMA_ROOT AND NOT $ENV{MAGMA_ROOT} STREQUAL "" )
    set( MAGMA_ROOT $ENV{MAGMA_ROOT} )
    set(MAGMA_LIBRARY_DIRS ${MAGMA_ROOT}/lib)
    set(MAGMA_INCLUDE_DIRS ${MAGMA_ROOT}/include)
	find_library(
      MAGMA_LIBRARIES
      NAMES magma
      PATHS ${MAGMA_ROOT}
      PATH_SUFFIXES lib
      NO_DEFAULT_PATH
	  )
    set(MAGMA_FOUND TRUE)
else()
    set(MAGMA_FOUND FALSE)
endif()
