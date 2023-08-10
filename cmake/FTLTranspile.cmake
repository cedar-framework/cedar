function (FTLTranspile output_sources input_sources)
  set(output_sources_local)
  foreach(input_path ${${input_sources}})
    file(REAL_PATH ${input_path} full_path BASE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
    file(RELATIVE_PATH rel_path ${CMAKE_SOURCE_DIR} ${full_path})
    get_filename_component(file_name ${full_path} NAME_WE)
    set(output_cpp ${CMAKE_CURRENT_BINARY_DIR}/${input_path}.cpp)
    set(output_hpp ${CMAKE_CURRENT_BINARY_DIR}/${input_path}.hpp)

    add_custom_command(OUTPUT ${output_cpp}
      COMMAND ${ftl_exe} ${full_path} -i ftl/Cedar.hpp -o ${output_cpp} --no-cpu
      DEPENDS ${full_path}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Transpiling Fortran source ${full_path}"
      VERBATIM
    )

    add_custom_command(OUTPUT ${output_hpp}
      COMMAND ${ftl_exe} ${full_path} -d -o ${output_hpp}
      DEPENDS ${full_path}
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Generating header from function definitions in ${rel_path}"
      VERBATIM
    )

    list(APPEND output_sources_local ${output_cpp})
    list(APPEND output_sources_local ${output_hpp})
  endforeach()
  set(output_combined)
  list(APPEND output_combined ${${output_sources}})
  list(APPEND output_combined ${output_sources_local})

  set(${output_sources} ${output_combined} PARENT_SCOPE)
endfunction()
