DEBUG?=0
D?=2
R?=32
CXX=mpic++
FLAGS=-std=c++20 -DPTDIM=$(D) -DFPSIZE=$(R)
INCS=-I./include

ifeq ($(DEBUG),1)
FLAGS+=-DDEBUG -O0 -g -fsanitize=address -fno-omit-frame-pointer
else
FLAGS+=-O2
endif

all: genpoints build_rgraph graph_diff point_info

genpoints: genpoints.cpp include
	$(CXX) -o $@ $(FLAGS) $(INCS) $<

build_rgraph: build_rgraph.cpp include
	$(CXX) -o $@ $(FLAGS) $(INCS) $<

graph_diff: graph_diff.cpp include
	$(CXX) -o $@ $(FLAGS) $(INCS) $<

point_info: point_info.cpp include
	$(CXX) -o $@ $(FLAGS) $(INCS) $<

clean:
	rm -rf *.dSYM genpoints build_rgraph graph_diff point_info
