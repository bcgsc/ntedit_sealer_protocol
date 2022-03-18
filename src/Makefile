btllib_path=$(shell dirname $(realpath $(MAKEFILE_LIST)))/../btllib

CXX=g++
CXXFLAGS=-I $(btllib_path)/include -fopenmp -pthread -O3
LDFLAGS=-L $(btllib_path)/lib -lbtllib

%: %.cpp
	$(CXX) $< -o $@ $(CXXFLAGS) $(LDFLAGS)

bins: build_index get_contig_reads build_targeted_bfs

.PHONY: clean

clean:
	$(RM) build_index get_contig_reads build_targeted_bfs