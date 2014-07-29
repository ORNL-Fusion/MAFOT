# ---- include platform-specific options ----
include ./make.inc


# ---- set directories ----
OBJDIR := $(MAFOT_DIR)/release
DIRS = $(OBJDIR)/d3d $(OBJDIR)/iter $(OBJDIR)/nstx $(OBJDIR)/mast $(LIB_DIR) $(BIN_DIR)


# ---- M3DC1 setup ----
ifeq ($(M3DC1),True)
   LIBS = -L$(LIB_DIR) -ltrip3d -lla_string $(BLITZLIBS) $(FLIBS) $(M3DC1LIBS) $(HDF5LIBS)
   INCLUDE = -I$(MAFOT_DIR)/include $(BLITZINCLUDE) $(M3DC1INCLUDE)
   DEFINES = -Dlinux -Dm3dc1
else
   LIBS = -L$(LIB_DIR) -ltrip3d -lla_string $(BLITZLIBS) $(FLIBS)
   INCLUDE = -I$(MAFOT_DIR)/include $(BLITZINCLUDE)
   DEFINES = -Dlinux
endif


# ---- Sources ----
SRCS = la_string.cxx
OBJS = $(SRCS:.cxx=.o)
OBJS := $(addprefix $(OBJDIR)/, $(OBJS))

MPISRCS = laminar_mpi.cxx foot_mpi.cxx plot_mpi.cxx
MPIOBJS = $(MPISRCS:.cxx=.o)
MPIOBJS_D3D = $(addprefix $(OBJDIR)/d3d/, $(MPIOBJS))
MPIOBJS_ITER = $(addprefix $(OBJDIR)/iter/, $(MPIOBJS))
MPIOBJS_NSTX = $(addprefix $(OBJDIR)/nstx/, $(MPIOBJS))
MPIOBJS_MAST = $(addprefix $(OBJDIR)/mast/, $(MPIOBJS))

SERSRCS = fix.cxx man.cxx plot.cxx
SEROBJS = $(SERSRCS:.cxx=.o)
SEROBJS_D3D = $(addprefix $(OBJDIR)/d3d/, $(SEROBJS))
SEROBJS_ITER = $(addprefix $(OBJDIR)/iter/, $(SEROBJS))
SEROBJS_NSTX = $(addprefix $(OBJDIR)/nstx/, $(SEROBJS))
SEROBJS_MAST = $(addprefix $(OBJDIR)/mast/, $(SEROBJS))

FSRCS = biotloop.f circleb.f d3icoils.f nstxecgeom.f polygonb.f \
        d3ccoils.f d3pferr.f ellipints.f itericoilsgeom.f masteccoilsgeom.f masticoilsgeom.f
FOBJS = $(FSRCS:.f=.o)
FOBJS := $(addprefix $(OBJDIR)/, $(FOBJS))


# ---- Common Targets ----
all : $(DIRS) d3d iter nstx mast

.PHONY : d3d
d3d : dtplot dtfix dtman dtlaminar_mpi dtfoot_mpi dtplot_mpi $(DIRS)

.PHONY : iter
iter : iterplot iterfix iterman iterlaminar_mpi iterfoot_mpi iterplot_mpi $(DIRS)

.PHONY : nstx 
nstx : nstxplot nstxfix nstxman nstxlaminar_mpi nstxfoot_mpi nstxplot_mpi $(DIRS)

.PHONY : mast 
mast : mastplot mastfix mastman mastlaminar_mpi mastfoot_mpi mastplot_mpi $(DIRS)

libtrip3d.a : $(FOBJS) $(DIRS)
	$(ARCH) $(LIB_DIR)/$@ $(FOBJS)

libla_string.a : $(OBJS) $(DIRS)
	$(ARCH) $(LIB_DIR)/$@ $(OBJS)

clean :
	rm -rf $(OBJDIR)
	
$(DIRS) : 
	mkdir -p $@


# ---- D3D Targets ----
dtplot : $(OBJDIR)/d3d/plot.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/d3d/plot.o -o $(BIN_DIR)/$@ $(LIBS)

dtfix : $(OBJDIR)/d3d/fix.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/d3d/fix.o -o $(BIN_DIR)/$@ $(LIBS)

dtman : $(OBJDIR)/d3d/man.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/d3d/man.o -o $(BIN_DIR)/$@ $(LIBS)

dtlaminar_mpi : $(OBJDIR)/d3d/laminar_mpi.o libla_string.a libtrip3d.a
	$(CXX) -fopenmp $(LDFLAGS) $(OBJDIR)/d3d/laminar_mpi.o -o $(BIN_DIR)/$@ $(OMPLIBS) $(LIBS)

dtfoot_mpi : $(OBJDIR)/d3d/foot_mpi.o libla_string.a libtrip3d.a
	$(CXX) -fopenmp $(LDFLAGS) $(OBJDIR)/d3d/foot_mpi.o -o $(BIN_DIR)/$@ $(OMPLIBS) $(LIBS)

dtplot_mpi : $(OBJDIR)/d3d/plot_mpi.o libla_string.a libtrip3d.a
	$(CXX) -fopenmp $(LDFLAGS) $(OBJDIR)/d3d/plot_mpi.o -o $(BIN_DIR)/$@ $(OMPLIBS) $(LIBS)


# ---- ITER Targets ----
iterplot : $(OBJDIR)/iter/plot.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/iter/plot.o -o $(BIN_DIR)/$@ $(LIBS)

iterfix : $(OBJDIR)/iter/fix.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/iter/fix.o -o $(BIN_DIR)/$@ $(LIBS)

iterman : $(OBJDIR)/iter/man.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/iter/man.o -o $(BIN_DIR)/$@ $(LIBS)

iterlaminar_mpi : $(OBJDIR)/iter/laminar_mpi.o libla_string.a libtrip3d.a
	$(CXX) -fopenmp $(LDFLAGS) $(OBJDIR)/iter/laminar_mpi.o -o $(BIN_DIR)/$@ $(OMPLIBS) $(LIBS)

iterfoot_mpi : $(OBJDIR)/iter/foot_mpi.o libla_string.a libtrip3d.a
	$(CXX) -fopenmp $(LDFLAGS) $(OBJDIR)/iter/foot_mpi.o -o $(BIN_DIR)/$@ $(OMPLIBS) $(LIBS)

iterplot_mpi : $(OBJDIR)/iter/plot_mpi.o libla_string.a libtrip3d.a
	$(CXX) -fopenmp $(LDFLAGS) $(OBJDIR)/iter/plot_mpi.o -o $(BIN_DIR)/$@ $(OMPLIBS) $(LIBS)


# ---- NSTX Targets ----
nstxplot : $(OBJDIR)/nstx/plot.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/nstx/plot.o -o $(BIN_DIR)/$@ $(LIBS)

nstxfix : $(OBJDIR)/nstx/fix.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/nstx/fix.o -o $(BIN_DIR)/$@ $(LIBS)

nstxman : $(OBJDIR)/nstx/man.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/nstx/man.o -o $(BIN_DIR)/$@ $(LIBS)

nstxlaminar_mpi : $(OBJDIR)/nstx/laminar_mpi.o libla_string.a libtrip3d.a
	$(CXX) -fopenmp $(LDFLAGS) $(OBJDIR)/nstx/laminar_mpi.o -o $(BIN_DIR)/$@ $(OMPLIBS) $(LIBS)

nstxfoot_mpi : $(OBJDIR)/nstx/foot_mpi.o libla_string.a libtrip3d.a
	$(CXX) -fopenmp $(LDFLAGS) $(OBJDIR)/nstx/foot_mpi.o -o $(BIN_DIR)/$@ $(OMPLIBS) $(LIBS)

nstxplot_mpi : $(OBJDIR)/nstx/plot_mpi.o libla_string.a libtrip3d.a
	$(CXX) -fopenmp $(LDFLAGS) $(OBJDIR)/nstx/plot_mpi.o -o $(BIN_DIR)/$@ $(OMPLIBS) $(LIBS)


# ---- MAST Targets ----
mastplot : $(OBJDIR)/mast/plot.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/mast/plot.o -o $(BIN_DIR)/$@ $(LIBS)

mastfix : $(OBJDIR)/mast/fix.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/mast/fix.o -o $(BIN_DIR)/$@ $(LIBS)

mastman : $(OBJDIR)/mast/man.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/mast/man.o -o $(BIN_DIR)/$@ $(LIBS)

mastlaminar_mpi : $(OBJDIR)/mast/laminar_mpi.o libla_string.a libtrip3d.a
	$(CXX) -fopenmp $(LDFLAGS) $(OBJDIR)/mast/laminar_mpi.o -o $(BIN_DIR)/$@ $(OMPLIBS) $(LIBS)

mastfoot_mpi : $(OBJDIR)/mast/foot_mpi.o libla_string.a libtrip3d.a
	$(CXX) -fopenmp $(LDFLAGS) $(OBJDIR)/mast/foot_mpi.o -o $(BIN_DIR)/$@ $(OMPLIBS) $(LIBS)

mastplot_mpi : $(OBJDIR)/mast/plot_mpi.o libla_string.a libtrip3d.a
	$(CXX) -fopenmp $(LDFLAGS) $(OBJDIR)/mast/plot_mpi.o -o $(BIN_DIR)/$@ $(OMPLIBS) $(LIBS)


# ---- Compile ----
$(OBJS) : $(OBJDIR)/%.o : $(MAFOT_DIR)/src/libla_string/%.cxx
	$(CXX) -c $(CFLAGS) $(INCLUDE) $(DEFINES) $< -o $@

$(FOBJS) : $(OBJDIR)/%.o : $(MAFOT_DIR)/src/libtrip3d/%.f
	$(F90) -c $(F90FLAGS) $< -o $@

$(MPIOBJS_D3D) : $(OBJDIR)/d3d/%.o : $(MAFOT_DIR)/src/%.cxx
	$(CXX) -c $(CFLAGS) $(OMPFLAGS) $(INCLUDE) $(OMPINCLUDE) $(DEFINES) -DD3D $< -o $@

$(MPIOBJS_ITER) : $(OBJDIR)/iter/%.o : $(MAFOT_DIR)/src/%.cxx
	$(CXX) -c $(CFLAGS) $(OMPFLAGS) $(INCLUDE) $(OMPINCLUDE) $(DEFINES) -DITER $< -o $@

$(MPIOBJS_MAST) : $(OBJDIR)/mast/%.o : $(MAFOT_DIR)/src/%.cxx
	$(CXX) -c $(CFLAGS) $(OMPFLAGS) $(INCLUDE) $(OMPINCLUDE) $(DEFINES) -DMAST $< -o $@

$(MPIOBJS_NSTX) : $(OBJDIR)/nstx/%.o : $(MAFOT_DIR)/src/%.cxx
	$(CXX) -c $(CFLAGS) $(OMPFLAGS) $(INCLUDE) $(OMPINCLUDE) $(DEFINES) -DNSTX $< -o $@

$(SEROBJS_D3D) : $(OBJDIR)/d3d/%.o : $(MAFOT_DIR)/src/%.cxx
	$(CXX) -c $(CFLAGS) $(INCLUDE) $(DEFINES) -DD3D $< -o $@

$(SEROBJS_ITER) : $(OBJDIR)/iter/%.o : $(MAFOT_DIR)/src/%.cxx
	$(CXX) -c $(CFLAGS) $(INCLUDE) $(DEFINES) -DITER $< -o $@

$(SEROBJS_MAST) : $(OBJDIR)/mast/%.o : $(MAFOT_DIR)/src/%.cxx
	$(CXX) -c $(CFLAGS) $(INCLUDE) $(DEFINES) -DMAST $< -o $@

$(SEROBJS_NSTX) : $(OBJDIR)/nstx/%.o : $(MAFOT_DIR)/src/%.cxx
	$(CXX) -c $(CFLAGS) $(INCLUDE) $(DEFINES) -DNSTX $< -o $@

