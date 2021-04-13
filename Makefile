# Compiler
CC := g++

# Directories
SRCDIR := src
INCDIR := include
BINDIR := bin
BUILDDIR := build

# Building this file is the main objective
TARGET := $(BINDIR)/simulate_pellet_combustion

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

$(BUILDDIR)/main.o : $(SRCDIR)/main.cpp $(INCDIR)/Combustion_Problem.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling main...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/main.cpp -o $(BUILDDIR)/main.o

$(BUILDDIR)/Combustion_Problem.o : $(SRCDIR)/Combustion_Problem.cpp $(INCDIR)/Combustion_Problem.hpp $(INCDIR)/TDMA_solver.hpp $(INCDIR)/Kinetics/Reaction.hpp $(INCDIR)/Thermo_Physical_Properties/Combustion_Pellet.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling Combustion_Problem...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Combustion_Problem.cpp -o $(BUILDDIR)/Combustion_Problem.o

$(BUILDDIR)/TDMA_solver.o : $(SRCDIR)/TDMA_solver.cpp $(INCDIR)/TDMA_solver.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling TDMA_solver...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/TDMA_solver.cpp -o $(BUILDDIR)/TDMA_solver.o

$(BUILDDIR)/Combustion_Pellet.o : $(SRCDIR)/Combustion_Pellet.cpp $(INCDIR)/Thermo_Physical_Properties/Combustion_Pellet.hpp $(INCDIR)/Thermo_Physical_Properties/Pellet_Properties.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling Combustion_Pellet...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Combustion_Pellet.cpp -o $(BUILDDIR)/Combustion_Pellet.o

$(BUILDDIR)/Pellet_Properties.o : $(SRCDIR)/Pellet_Properties.cpp $(INCDIR)/Thermo_Physical_Properties/Pellet_Properties.hpp $(INCDIR)/Thermo_Physical_Properties/Coated_Particle.hpp $(INCDIR)/Thermo_Physical_Properties/Substance.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling Pellet_Properties...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Pellet_Properties.cpp -o $(BUILDDIR)/Pellet_Properties.o

$(BUILDDIR)/Coated_Particle.o : $(SRCDIR)/Coated_Particle.cpp $(INCDIR)/Thermo_Physical_Properties/Coated_Particle.hpp $(INCDIR)/Thermo_Physical_Properties/Substance.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling Coated_Particle...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Coated_Particle.cpp -o $(BUILDDIR)/Coated_Particle.o

$(BUILDDIR)/Substance.o : $(SRCDIR)/Substance.cpp $(INCDIR)/Thermo_Physical_Properties/Substance.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling Substance...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Substance.cpp -o $(BUILDDIR)/Substance.o

$(BUILDDIR)/Reaction.o : $(SRCDIR)/Reaction.cpp $(INCDIR)/Kinetics/Reaction.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling Reaction...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Reaction.cpp -o $(BUILDDIR)/Reaction.o

$(BUILDDIR)/Kinetics.o : $(SRCDIR)/Kinetics.cpp $(INCDIR)/Kinetics/Kinetics.hpp
	@mkdir -p $(BUILDDIR);
	@echo "\nCompiling Kinetics...";
	$(CC) $(CFLAGS) $(INC) -c $(SRCDIR)/Kinetics.cpp -o $(BUILDDIR)/Kinetics.o

clean:
	@echo "\nCleaning..."; 
	@echo "$(RM) -r $(BUILDDIR)/*.o $(TARGET)"; $(RM) -r $(BUILDDIR) $(BINDIR)