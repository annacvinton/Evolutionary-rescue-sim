CXX=g++

CPPFLAGS=-I. -std=c++11 -O2
DEPS = Population.hpp Individual.hpp SpatialContainer.hpp
OBJS = main.o Population.o Individual.o SpatialContainer.o

%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CPPFLAGS)

iiasa: $(OBJS)
	$(CXX) -o $@ $^ $(CPPFLAGS)
