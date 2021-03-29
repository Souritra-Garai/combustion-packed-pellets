# Compiler
CC := g++

# Directories
SRCDIR := src
INCDIR := include
BINDIR := bin
BUILDDIR := build

# Building this file is the main objective
TARGET := $(BINDIR)/simulate_pellet_diffusion

# Finding all source and object files
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) $(LIBDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

# Flags required for compiler
CFLAGS := -fopenmp
LIB := -lm
INC := -I $(INCDIR)

$(TARGET) : $(OBJECTS)
	@mkdir -p $(BINDIR);
	@echo "\nLinking all...";
	$(CC) $(OBJECTS) $(CFLAGS) $(LIB) -o $(TARGET)

$(BUILDDIR)/main.o : $(SRCDIR)/main.cpp $(INCDIR)/Temperature_Iterator.hpp $(INCDIR)/Thermodynamic_Properties.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling main...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/main.cpp -o $(BUILDDIR)/main.o

$(BUILDDIR)/Temperature_Iterator.o : $(SRCDIR)/Temperature_Iterator.cpp $(INCDIR)/Temperature_Iterator.hpp $(INCDIR)/TDMA_solver.hpp $(INCDIR)/Kinetics.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling Temperature_Iterator...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Temperature_Iterator.cpp -o $(BUILDDIR)/Temperature_Iterator.o

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

clean:
	@echo "\nCleaning..."; 
	@echo "$(RM) -r $(BUILDDIR)/*.o $(TARGET)"; $(RM) -r $(BUILDDIR) $(BINDIR)