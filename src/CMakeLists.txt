set(SOURCES  field/Var.cpp 
             grid/Cell.cpp 
             grid/Quad.cpp  
             grid/Grid.cpp  
             grid/operator.cpp  
             grid/operator2.cpp  
             grid/Vertex.cpp
             grid/Line.cpp
	     grid/Block1.cpp
             )

add_library(girdapcore STATIC ${SOURCES})
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)
