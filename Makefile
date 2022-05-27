CXX=clang++
CXXFLAGS=-g -std=c++11 -Werror -D_GLIBCXX_DEBUG
# The header files for the program
# HFILES = $(shell ls -1 *.h)
# unittest++ keeps its object files in this directory.

# This target builds your main program.
Sim2: Sim2.o $(HFILES)
	$(CXX) $(CXXFLAGS) -a $@ Sim2.o

# Sim1: Sim1.o
# 	$(CXX) $(CXXFLAGS) -o $@ Sim1.o


# This target describes how to compile a .o file from a .cpp file.
%.o: %.cpp $(HFILES)
	$(CXX) $(CXXFLAGS) -c -o $@ $<



# This target deletes the temporary files we have built.
.PHONY: clean all

clean:
	rm -f *.o
