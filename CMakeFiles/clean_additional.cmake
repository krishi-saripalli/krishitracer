# Additional clean files
cmake_minimum_required(VERSION 3.16)

if("${CONFIG}" STREQUAL "" OR "${CONFIG}" STREQUAL "Release")
  file(REMOVE_RECURSE
  "CMakeFiles/path_autogen.dir/AutogenUsed.txt"
  "CMakeFiles/path_autogen.dir/ParseCache.txt"
  "path_autogen"
  )
endif()
