# Compiler
FC = mpif90

# Paths (Modify if necessary)
CUDECOMP_DIR = /home/milton/MHIT36_cuDecomp/cuDecomp-main/build/

# Compiler and Linker Flags
FFLAGS = -cuda -I$(CUDECOMP_DIR)/include
LDFLAGS = -L$(CUDECOMP_DIR)/lib -L/usr/local/cuda/lib64
LIBS = -lcudecomp_fort -lcudecomp -cuda

# Source and Object Files
SRCS = basic_usage.f90
OBJS = $(SRCS:.f90=.o)

# Output Executable
TARGET = my_program

# Default Rule
all: $(TARGET)

# Linking
$(TARGET): $(OBJS)
	$(FC) $(LDFLAGS) -o $(TARGET) $(OBJS) $(LIBS)

# Compilation Rule
%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

# Clean Rule
clean:
	rm -f $(TARGET) $(OBJS) *.mod
