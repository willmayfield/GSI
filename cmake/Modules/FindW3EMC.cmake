# This module looks for environment variables detailing where W3EMC lib is
# If variables are not set, W3EMC will be built from external source 
include(ExternalProject)
if(NOT BUILD_W3EMC )
  if(DEFINED ENV{W3EMC_LIB4})
    set(W3EMC_LIBRARY $ENV{W3EMC_LIB4} )
    set(W3EMC_LIBRARY $ENV{W3EMC_LIBd} )
    set(W3EMCINC $ENV{W3EMC_INCd} )
    set(W3EMC4INC $ENV{W3EMC_INC4} )
    if( CORE_LIBRARIES )
      list( APPEND CORE_LIBRARIES ${W3EMC_LIBRARY} )
      list( APPEND CORE_INCS ${W3EMCINC} )
    else()
      set( CORE_LIBRARIES ${W3EMC_LIBRARY} )
      set( CORE_INCS ${W3EMCINC} )
    endif()
  endif()
else()
  set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})
  ExternalProject_Add(NCEPLIBS-w3emc 
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/libsrc/w3emc 
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
  execute_process(COMMAND grep "set(VERSION" CMakeLists.txt WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/libsrc/w3emc OUTPUT_VARIABLE LIBVERSION)
  string(REPLACE "set(VERSION " "" LIBVERSION ${LIBVERSION})
  string(REPLACE ")" "" LIBVERSION ${LIBVERSION})
  string(REPLACE "\n" "" LIBVERSION ${LIBVERSION})
  message("w3emc version is ${LIBVERSION}")
  set( W3EMC_LIBRARY ${PROJECT_BINARY_DIR}/lib/libw3emc_${LIBVERSION}_d.a )
  set( W3EMC_4_LIBRARY ${PROJECT_BINARY_DIR}/lib/libw3emc_${LIBVERSION}_4.a )
  if( CORE_BUILT )
      list( APPEND CORE_BUILT ${W3EMC_LIBRARY} )
      list( APPEND EXT_BUILT NCEPLIBS-w3emc )
  else()
      set( CORE_BUILT ${W3EMC_LIBRARY} )
      set( EXT_BUILT NCEPLIBS-w3emc )
  endif()
endif( )

set( W3EMC_LIBRARY_PATH ${W3EMC_LIBRARY} CACHE STRING "W3EMC Library Location" )
set( W3EMC_4_LIBRARY_PATH ${W3EMC_4_LIBRARY} CACHE STRING "W3EMC_4 Library Location" )
set( W3EMC_INCLUDE_PATH ${W3EMCINC} CACHE STRING "W3EMC Include Location" )
set( W3EMC_INCLUDE_4_PATH ${W3EMC4INC} CACHE STRING "W3EMC_4 Include Location" )

