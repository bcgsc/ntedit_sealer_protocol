btllib_path=$(shell dirname $(realpath $(MAKEFILE_LIST)))/../btllib

CXX=g++
CXXFLAGS=-I $(btllib_path)/install/include -fopenmp -pthread -O3 -std=c++17
LDFLAGS=-L $(btllib_path)/install/lib -lbtllib

%: %.cpp utils.cpp
	$(CXX) $^ -o $@ $(CXXFLAGS) $(LDFLAGS)

bins: build_index build_targeted_bfs

.PHONY: clean

clean:
	$(RM) build_index build_targeted_bfs
