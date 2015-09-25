#! /bin/sh
SOURCES='src/main_triint.cpp 
           include/field/Var.cpp
           include/grid2/Cell.cpp
           include/grid2/Quad.cpp
           include/grid2/Grid.cpp
           include/grid2/operator.cpp
           include/grid2/Vertex.cpp
           include/grid2/Line.cpp'
#           include/grid1/Hexa.cpp'
#          include/grid2/Vertex.cpp 
#          include/grid2/Cell.cpp 
#          include/grid2/Grid.cpp'
#         include/grid/type/Line.cpp include/grid/type/Tri.cpp
#         include/grid/type/Hexa.cpp'
#c++ -g -I include/ -std=c++0x src/main1.cpp include/grid/Vertex.cpp include/grid/Cell.cpp include/grid/quad/Quad.cpp -o cart3
#c++ -I include/ -std=c++11 $SOURCES -o cart3
clang++ -g -std=c++11 -stdlib=libc++ -I include $SOURCES -o cart3
