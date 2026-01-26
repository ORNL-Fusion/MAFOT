# ---- macOS SDK (needed so libc headers like <time.h> are found) ----
SYSROOT := $(shell xcrun --sdk macosx --show-sdk-path)
# Align deployment target across C++, Fortran, and link steps
SDKVER := $(shell xcrun --sdk macosx --show-sdk-version)
export MACOSX_DEPLOYMENT_TARGET := $(SDKVER)

CFLAGS   += -mmacosx-version-min=$(SDKVER)
F90FLAGS += -mmacosx-version-min=$(SDKVER)
LDFLAGS  += -mmacosx-version-min=$(SDKVER)

# ---- Folders ----
MAFOT_DIR  = $(PWD)
BIN_DIR    = $(MAFOT_DIR)/build/bin
LIB_DIR    = $(MAFOT_DIR)/build/lib

# ---- Support codes ----
M3DC1 = False
VMEC  = False
SIESTA= False

# ---- Compilers ----
# Use Open-MPI's wrapper; we’ll point mpicxx at Homebrew LLVM via OMPI_CXX.
CXX      = /opt/homebrew/bin/mpicxx
F90      = /opt/homebrew/bin/gfortran

# ---- Fortran runtime (portable discovery) ----
# Find where Homebrew put libgfortran and libgcc_s
GFORTRAN      := $(F90)
LIBGFORTRAN   := $(shell $(GFORTRAN) -print-file-name=libgfortran.dylib)
LIBGCC_S      := $(shell $(GFORTRAN) -print-file-name=libgcc_s.1.dylib)
FLIBDIR       := $(dir $(LIBGFORTRAN))
LIBGCCDIR     := $(dir $(LIBGCC_S))

# Link flags & rpaths so the executable can find them at runtime
FLIBS         := -L$(FLIBDIR) -lgfortran -Wl,-rpath,$(FLIBDIR)
# Some GCC builds also need libgcc_s at runtime; include if it exists
ifneq ($(LIBGCC_S),libgcc_s.1.dylib)
  FLIBS += -L$(LIBGCCDIR) -lgcc_s.1 -Wl,-rpath,$(LIBGCCDIR)
endif


# ---- Compilation flags (the recipe uses CFLAGS, so put everything here) ----
CFLAGS   = -O3 -fPIC -std=c++11 \
           -isysroot $(SYSROOT) \
           -Xpreprocessor -fopenmp -w \
           -Ulinux
F90FLAGS = -O3 -fPIC -isysroot $(SYSROOT)

# ---- Linker ----
ARCH     = ar crus
LDD      = $(CXX) -dynamiclib
LDFLAGS  = -isysroot $(SYSROOT)

# ---- Additional Defines ----
DEFS     =

# Blitz from Homebrew
BLITZINCLUDE = -I$(MAFOT_DIR)/third_party/blitz++-0.9.orig/
BLITZLIBS    = -L$(MAFOT_DIR)/third_party/blitz++-0.9.orig/lib -lblitz

# OpenMP runtime for clang
OMPLIBS      = -L/opt/homebrew/opt/libomp/lib -lomp

# HDF5 / NetCDF (C)
HDF5_CFLAGS   := $(shell pkg-config --cflags hdf5)
HDF5_LIBS     := $(shell pkg-config --libs   hdf5)
NETCDF_CFLAGS := $(shell pkg-config --cflags netcdf)
NETCDF_LIBS   := $(shell pkg-config --libs   netcdf)

# NetCDF-Fortran
NETCDFF_FLAGS := $(shell nf-config --fflags)
NETCDFF_LIBS  := $(shell nf-config --flibs)

# Variables MAFOT makefiles typically consume
NETCDFINCLUDE = $(NETCDF_CFLAGS)
HDF5LIBS      = $(HDF5_LIBS)
NETCDFLIBS    = $(NETCDF_LIBS) $(NETCDFF_LIBS)

# MPI include/libs come from mpicxx automatically
OMPINCLUDE    = -I/opt/homebrew/opt/libomp/include
M3DC1INCLUDE  =
M3DC1LIBS     =

