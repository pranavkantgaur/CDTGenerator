g++ -std=c++0x -g rply/rply.cpp cdtGen.cpp -I /usr/include/vtk-5.8 -I /usr/local/include/pcl-1.7/ -lgmp -lmpfr -lCGAL -lpcl_common -lpcl_visualization -lboost_system -lpcl_io  -o output
