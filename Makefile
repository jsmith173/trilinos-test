include $(HOME)/Trilinos-install/include/Makefile.export.Trilinos

CC = $(Trilinos_C_COMPILER)
CFLAGS = $(Trilinos_C_COMPILER_FLAGS)

CXX = $(Trilinos_CXX_COMPILER)
CPPFLAGS = $(Trilinos_INCLUDE_DIRS)
CXXFLAGS = $(Trilinos_CXX_COMPILER_FLAGS) 
LDFLAGS = $(Trilinos_LIBRARY_DIRS)
LDLIBS = $(Trilinos_LIBRARIES)
	

cxx_main: cxx_main.o 
	$(CXX) $(LDFLAGS) $(LDLIBS) -o cxx_main cxx_main.o 

cxx_main.o: cxx_main.cpp 
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c cxx_main.cpp

