# Compiler
CC := g++

# Building this file is the main objective
TARGET := bin/simulate_pellet_combustion

# Directories
SRCDIR := src
INCDIR := include
BINDIR := bin
TESTDIR := test
BUILDDIR := build

# Finding all source and object files
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) $(LIBDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

# Flags required for compiler
CFLAGS := -fopenmp
LIB := -lm
INC := -I $(INCDIR)

$(TARGET) : $(OBJECTS)
	@echo "\nLinking all...";
	$(CC) $(OBJECTS) $(CFLAGS) $(LIB) -o $(TARGET)

$(BUILDDIR)/main.o : $(SRCDIR)/main.cpp $(INCDIR)/Kinetics.hpp $(INCDIR)/TDMA_solver.hpp $(INCDIR)/Thermodynamic_Properties.hpp $(INCDIR)/Temperature_Step_Iterator.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling main...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/main.cpp -o $(BUILDDIR)/main.o

$(BUILDDIR)/TDMA_solver.o : $(SRCDIR)/TDMA_solver.cpp $(INCDIR)/TDMA_solver.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling TDMA_solver...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/TDMA_solver.cpp -o $(BUILDDIR)/TDMA_solver.o

$(BUILDDIR)/Thermodynamic_Properties.o : $(SRCDIR)/Thermodynamic_Properties.cpp $(INCDIR)/Thermodynamic_Properties.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling Thermodynamic_Properties...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Thermodynamic_Properties.cpp -o $(BUILDDIR)/Thermodynamic_Properties.o

$(BUILDDIR)/Kinetics.o : $(SRCDIR)/Kinetics.cpp $(INCDIR)/Kinetics.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling Kinetics...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Kinetics.cpp -o $(BUILDDIR)/Kinetics.o

$(BUILDDIR)/Temperature_Step_Iterator.o : $(SRCDIR)/Temperature_Step_Iterator.cpp $(INCDIR)/Temperature_Step_Iterator.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling Temperature_Step_Iterator...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Temperature_Step_Iterator.cpp -o $(BUILDDIR)/Temperature_Step_Iterator.o

# Builds the example for TDMA solver
test_TDMA_solver : $(BUILDDIR)/TDMA_example.o $(BUILDDIR)/TDMA_solver.o
	@echo "\nBuilding TDMA_example...";
	$(CC) $(CFLAGS) $(LIB) $(BUILDDIR)/TDMA_example.o $(BUILDDIR)/TDMA_solver.o -o $(BINDIR)/TDMA_example

$(BUILDDIR)/TDMA_example.o : $(TESTDIR)/TDMA_example.cpp $(INCDIR)/TDMA_solver.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling TDMA_example...";
	$(CC) $(CFLAGS) $(INC) -c $(TESTDIR)/TDMA_example.cpp -o $(BUILDDIR)/TDMA_example.o

clean:
	@echo "\nCleaning..."; 
	@echo "$(RM) -r $(BUILDDIR)/*.o $(TARGET)"; $(RM) -r $(BUILDDIR)/* $(TARGET)