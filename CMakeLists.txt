cmake_minimum_required(VERSION 2.8)
project(girdap)


if (CMAKE_VERSION VERSION_LESS "3.1")
  if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
  endif ()
else ()
  set (CMAKE_CXX_STANDARD 11)
endif ()


include_directories(include)

set(SOURCES  ./include/field/Var.cpp 
             ./include/grid2/Cell.cpp 
             ./include/grid2/Quad.cpp  
             ./include/grid2/Grid.cpp  
             ./include/grid2/operator.cpp  
             ./include/grid2/Vertex.cpp
             ./include/grid2/Line.cpp
             )

add_executable(girdap src/main.cpp ${SOURCES})

install(TARGETS girdap RUNTIME DESTINATION bin)
