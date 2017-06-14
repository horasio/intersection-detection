
cmake_minimum_required(VERSION 2.6)

# set default cmake build type to Release (None Debug Release RelWithDebInfo MinSizeRel)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(CPP11_FLAG -std=c++11)

if( APPLE )
	set(LIBS ${LIBS} "-L/opt/local/lib")
	include_directories(/opt/local/include)
else()
	#set(LIBS ${LIBS} "-L/usr/local/lib")
endif()

set(LIBS ${LIBS} "-lsqlite3" "-lmpfr" "-lgmp")

find_package(CGAL)

if (CGAL_FOUND)
  include_directories(${CGAL_INCLUDE_DIRS})

  set(TET_SOURCES
    convexes.cpp
    gjk.cpp
    kmeans.cpp
    tettest.cpp
    vec.cpp
    CGALStuff.cpp
    )
  set_source_files_properties(${TET_SOURCES} PROPERTIES COMPILE_FLAGS "${CPP11_FLAG}")
  add_executable(tet ${TET_SOURCES})
  message("${CGAL_LIBRARY}")
  target_link_libraries (tet ${CGAL_LIBRARY} ${LIBS})

endif (CGAL_FOUND)
