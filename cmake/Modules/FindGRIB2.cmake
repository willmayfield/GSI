# - Find the GRIB2 modules
set( NO_DEFAULT_PATH )

if(DEFINED ENV{G2_LIB4})
  set(G2_LIB4 $ENV{G2_LIB4} )
else()
  find_library( G2_LIB4  
      NAMES libg2_4.a
      HINTS
        $ENV{GRIB2PATH}/g2/lib
      ${NO_DEFAULT_PATH})
endif()

if(DEFINED ENV{G2TMPL_LIB})
  set(G2TMPL_LIB $ENV{G2TMPL_LIB} )
else()
  find_library( G2TMPL_LIB
      NAMES libg2tmpl.a
      HINTS
        $ENV{GRIB2PATH}/g2tmpl/lib
      ${NO_DEFAULT_PATH})
endif()

if(DEFINED ENV{GRIB2INC})
  set(GRIB2INC $ENV{GRIB2INC} )
else()
find_path( GRIB2INC
    NAMES grib_mod.mod
    HINTS
       $ENV{GRIB2PATH}/g2/include
    ${NO_DEFAULT_PATH})
endif()

set( GRIB2_LIBRARY ${G2_LIB4} ${G2TMPL_LIB} CACHE STRING "GRIB2 Library Location" )
set( GRIB2_INCLUDE_PATH ${GRIB2INC} CACHE STRING "GRIB2 Include Location" )

# - Find the JASPER modules
set( NO_DEFAULT_PATH )
find_library( JASPER_LIB
    NAMES libjasper.so
    HINTS
       /usr/lib64
    ${NO_DEFAULT_PATH})

# - Find the PNG modules
set( NO_DEFAULT_PATH )
find_library( PNG_LIB
    NAMES libpng.so
    HINTS
       /usr/lib64
    ${NO_DEFAULT_PATH})

# - Find the Z modules
set( NO_DEFAULT_PATH )
find_library( Z_LIB
    NAMES libz.so
    HINTS
       /usr/lib64
    ${NO_DEFAULT_PATH})
