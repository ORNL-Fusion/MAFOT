# ---- include platform-specific options ----
include ./make.inc

#CFLAGS += -std=c++11


# ---- set directories ----
OBJDIR := $(MAFOT_DIR)/release
DIRS = $(OBJDIR)/d3d $(OBJDIR)/iter $(OBJDIR)/nstx $(OBJDIR)/mast $(OBJDIR)/cmod $(LIB_DIR) $(BIN_DIR)


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


# ---- Other Defines ----
ifdef DEFS
   DEFINES += $(DEFS)
endif


# ---- SIESTA support ----
ifdef SIESTA
ifeq ($(SIESTA),True)
   DEFINES += -DUSE_SIESTA
   VMEC = True
endif
endif


# ---- VMEC support ----
ifdef VMEC
ifeq ($(VMEC),True)
   LIBS += $(NETCDFLIBS)
   INCLUDE += $(NETCDFINCLUDE)
   DEFINES += -DUSE_XFIELD
endif
endif


# ---- Sources ----
SRCS = la_string.cxx
OBJS = $(SRCS:.cxx=.o)
OBJS := $(addprefix $(OBJDIR)/, $(OBJS))
DEPS = $(OBJS:.o=.d)

MPISRCS = laminar_mpi.cxx foot_mpi.cxx plot_mpi.cxx
MPIOBJS = $(MPISRCS:.cxx=.o)
MPIOBJS_D3D = $(addprefix $(OBJDIR)/d3d/, $(MPIOBJS))
MPIOBJS_ITER = $(addprefix $(OBJDIR)/iter/, $(MPIOBJS))
MPIOBJS_NSTX = $(addprefix $(OBJDIR)/nstx/, $(MPIOBJS))
MPIOBJS_MAST = $(addprefix $(OBJDIR)/mast/, $(MPIOBJS))
MPIOBJS_CMOD = $(addprefix $(OBJDIR)/cmod/, $(MPIOBJS))
MPIDEPS_D3D = $(MPIOBJS_D3D:.o=.d)
MPIDEPS_ITER = $(MPIOBJS_ITER:.o=.d)
MPIDEPS_NSTX = $(MPIOBJS_NSTX:.o=.d)
MPIDEPS_MAST = $(MPIOBJS_MAST:.o=.d)
MPIDEPS_CMOD = $(MPIOBJS_CMOD:.o=.d)

SERSRCS = fix.cxx man.cxx plot.cxx structure.cxx
SEROBJS = $(SERSRCS:.cxx=.o)
SEROBJS_D3D = $(addprefix $(OBJDIR)/d3d/, $(SEROBJS))
SEROBJS_ITER = $(addprefix $(OBJDIR)/iter/, $(SEROBJS))
SEROBJS_NSTX = $(addprefix $(OBJDIR)/nstx/, $(SEROBJS))
SEROBJS_MAST = $(addprefix $(OBJDIR)/mast/, $(SEROBJS))
SEROBJS_CMOD = $(addprefix $(OBJDIR)/cmod/, $(SEROBJS))
SERDEPS_D3D = $(SEROBJS_D3D:.o=.d)
SERDEPS_ITER = $(SEROBJS_ITER:.o=.d)
SERDEPS_NSTX = $(SEROBJS_NSTX:.o=.d)
SERDEPS_MAST = $(SEROBJS_MAST:.o=.d)
SERDEPS_CMOD = $(SEROBJS_CMOD:.o=.d)

FSRCS = biotloop.f circleb.f d3icoils.f nstxecgeom.f polygonb.f d3busgeom.f\
        d3ccoils.f d3pferr.f ellipints.f itericoilsgeom.f masteccoilsgeom.f masticoilsgeom.f
FOBJS = $(FSRCS:.f=.o)
FOBJS := $(addprefix $(OBJDIR)/, $(FOBJS))


# ---- Common Targets ----
all : $(DIRS) d3d iter nstx mast cmod gui xpand d3dplot

.PHONY : d3d
d3d : $(DIRS) dtplot dtfix dtman dtlaminar_mpi dtfoot_mpi dtplot_mpi dtstructure

.PHONY : iter
iter : $(DIRS) iterplot iterfix iterman iterlaminar_mpi iterfoot_mpi iterplot_mpi iterstructure

.PHONY : nstx 
nstx : $(DIRS) nstxplot nstxfix nstxman nstxlaminar_mpi nstxfoot_mpi nstxplot_mpi nstxstructure

.PHONY : mast 
mast : $(DIRS) mastplot mastfix mastman mastlaminar_mpi mastfoot_mpi mastplot_mpi maststructure

.PHONY : cmod 
cmod : $(DIRS) cmodplot cmodfix cmodman cmodlaminar_mpi cmodfoot_mpi cmodplot_mpi cmodstructure

.PHONY : plot
plot : $(DIRS) dtplot_mpi iterplot_mpi nstxplot_mpi mastplot_mpi cmodplot_mpi

.PHONY : fix
fix : $(DIRS) dtfix iterfix nstxfix mastfix cmodfix

.PHONY : man
man : $(DIRS) dtman iterman nstxman mastman cmodman

.PHONY : laminar
laminar : $(DIRS) dtlaminar_mpi iterlaminar_mpi nstxlaminar_mpi mastlaminar_mpi cmodlaminar_mpi

.PHONY : foot
foot : $(DIRS) dtfoot_mpi iterfoot_mpi nstxfoot_mpi mastfoot_mpi cmodfoot_mpi

.PHONY : structure
structure : $(DIRS) dtstructure iterstructure nstxstructure maststructure cmodstructure

.PHONY : gui
gui : $(MAFOT_DIR)/python/mafot_gui.py
	rm -f $(BIN_DIR)/mafot_gui.py
	ln $(MAFOT_DIR)/python/mafot_gui.py $(BIN_DIR)/mafot_gui.py
	chmod +x $(BIN_DIR)/mafot_gui.py

.PHONY : xpand
xpand : $(DIRS) xpand_mpi $(MAFOT_DIR)/python/use_xpand.py
	rm -f $(BIN_DIR)/use_xpand.py
	ln $(MAFOT_DIR)/python/use_xpand.py $(BIN_DIR)/use_xpand.py
	chmod +x $(BIN_DIR)/use_xpand.py

.PHONY : d3dplot
d3dplot : $(MAFOT_DIR)/python/d3dplot.py
	rm -f $(BIN_DIR)/d3dplot.py
	ln $(MAFOT_DIR)/python/d3dplot.py $(BIN_DIR)/d3dplot.py
	chmod +x $(BIN_DIR)/d3dplot.py

libtrip3d.a : $(DIRS) $(FOBJS) 
	$(ARCH) $(LIB_DIR)/$@ $(FOBJS)

libla_string.a : $(DIRS) $(OBJS) 
	$(ARCH) $(LIB_DIR)/$@ $(OBJS)

clean :
	rm -rf $(OBJDIR)
	
$(DIRS) : 
	mkdir -p $@


# ---- Targets ----
xpand_mpi : $(MAFOT_DIR)/src/xpand_mpi.cxx
ifdef VMEC
	$(CXX) -c $(CFLAGS) $(OMPFLAGS) $(INCLUDE) $(OMPINCLUDE) $(DEFINES) $< -o $(OBJDIR)/xpand_mpi.o
	$(CXX) -fopenmp $(LDFLAGS) $(OBJDIR)/xpand_mpi.o -o $(BIN_DIR)/$@ $(OMPLIBS) $(LIBS)
else
	@echo "-----------------------------------"
	@echo "VMEC not supported"
	@echo "-----------------------------------"
endif


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

dtstructure : $(OBJDIR)/d3d/structure.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/d3d/structure.o -o $(BIN_DIR)/$@ $(LIBS)
	
fi_prepare : $(OBJDIR)/d3d/fi_prepare.o libla_string.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/d3d/fi_prepare.o -o $(BIN_DIR)/$@ $(LIBS)

$(OBJDIR)/d3d/fi_prepare.o : $(OBJDIR)/d3d/%.o : $(MAFOT_DIR)/src/%.cxx
	$(CXX) -c $(CFLAGS) $(INCLUDE) $(DEFINES) -DD3D $< -o $@

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

iterstructure : $(OBJDIR)/iter/structure.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/iter/structure.o -o $(BIN_DIR)/$@ $(LIBS)


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

nstxstructure : $(OBJDIR)/nstx/structure.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/nstx/structure.o -o $(BIN_DIR)/$@ $(LIBS)


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

maststructure : $(OBJDIR)/mast/structure.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/mast/structure.o -o $(BIN_DIR)/$@ $(LIBS)


# ---- CMOD Targets ----
cmodplot : $(OBJDIR)/cmod/plot.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/cmod/plot.o -o $(BIN_DIR)/$@ $(LIBS)

cmodfix : $(OBJDIR)/cmod/fix.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/cmod/fix.o -o $(BIN_DIR)/$@ $(LIBS)

cmodman : $(OBJDIR)/cmod/man.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/cmod/man.o -o $(BIN_DIR)/$@ $(LIBS)

cmodlaminar_mpi : $(OBJDIR)/cmod/laminar_mpi.o libla_string.a libtrip3d.a
	$(CXX) -fopenmp $(LDFLAGS) $(OBJDIR)/cmod/laminar_mpi.o -o $(BIN_DIR)/$@ $(OMPLIBS) $(LIBS)

cmodfoot_mpi : $(OBJDIR)/cmod/foot_mpi.o libla_string.a libtrip3d.a
	$(CXX) -fopenmp $(LDFLAGS) $(OBJDIR)/cmod/foot_mpi.o -o $(BIN_DIR)/$@ $(OMPLIBS) $(LIBS)

cmodplot_mpi : $(OBJDIR)/cmod/plot_mpi.o libla_string.a libtrip3d.a
	$(CXX) -fopenmp $(LDFLAGS) $(OBJDIR)/cmod/plot_mpi.o -o $(BIN_DIR)/$@ $(OMPLIBS) $(LIBS)

cmodstructure : $(OBJDIR)/cmod/structure.o libla_string.a libtrip3d.a
	$(CXX) $(LDFLAGS) $(OBJDIR)/cmod/structure.o -o $(BIN_DIR)/$@ $(LIBS)


# ---- Include Dependencies ----
-include $(DEPS)
-include $(SERDEPS_D3D)
-include $(SERDEPS_ITER)
-include $(SERDEPS_NSTX)
-include $(SERDEPS_MAST)
-include $(SERDEPS_CMOD)
-include $(MPIDEPS_D3D)
-include $(MPIDEPS_ITER)
-include $(MPIDEPS_NSTX)
-include $(MPIDEPS_MAST)
-include $(MPIDEPS_CMOD)


# ---- Compile ----
$(OBJS) : $(OBJDIR)/%.o : $(MAFOT_DIR)/src/libla_string/%.cxx
	$(CXX) -c $(CFLAGS) $(INCLUDE) $(DEFINES) $< -o $@

$(FOBJS) : $(OBJDIR)/%.o : $(MAFOT_DIR)/src/libtrip3d/%.f
	$(F90) -c $(F90FLAGS) $< -o $@

$(MPIOBJS_D3D) : $(OBJDIR)/d3d/%.o : $(MAFOT_DIR)/src/%.cxx
	$(CXX) -c $(CFLAGS) -MMD $(OMPFLAGS) $(INCLUDE) $(OMPINCLUDE) $(DEFINES) -DD3D $< -o $@

$(MPIOBJS_ITER) : $(OBJDIR)/iter/%.o : $(MAFOT_DIR)/src/%.cxx
	$(CXX) -c $(CFLAGS) -MMD $(OMPFLAGS) $(INCLUDE) $(OMPINCLUDE) $(DEFINES) -DITER $< -o $@

$(MPIOBJS_MAST) : $(OBJDIR)/mast/%.o : $(MAFOT_DIR)/src/%.cxx
	$(CXX) -c $(CFLAGS) -MMD $(OMPFLAGS) $(INCLUDE) $(OMPINCLUDE) $(DEFINES) -DMAST $< -o $@

$(MPIOBJS_NSTX) : $(OBJDIR)/nstx/%.o : $(MAFOT_DIR)/src/%.cxx
	$(CXX) -c $(CFLAGS) -MMD $(OMPFLAGS) $(INCLUDE) $(OMPINCLUDE) $(DEFINES) -DNSTX $< -o $@
	
$(MPIOBJS_CMOD) : $(OBJDIR)/cmod/%.o : $(MAFOT_DIR)/src/%.cxx
	$(CXX) -c $(CFLAGS) -MMD $(OMPFLAGS) $(INCLUDE) $(OMPINCLUDE) $(DEFINES) -DCMOD $< -o $@

$(SEROBJS_D3D) : $(OBJDIR)/d3d/%.o : $(MAFOT_DIR)/src/%.cxx
	$(CXX) -c $(CFLAGS) -MMD $(INCLUDE) $(DEFINES) -DD3D $< -o $@

$(SEROBJS_ITER) : $(OBJDIR)/iter/%.o : $(MAFOT_DIR)/src/%.cxx
	$(CXX) -c $(CFLAGS) -MMD $(INCLUDE) $(DEFINES) -DITER $< -o $@

$(SEROBJS_MAST) : $(OBJDIR)/mast/%.o : $(MAFOT_DIR)/src/%.cxx
	$(CXX) -c $(CFLAGS) -MMD $(INCLUDE) $(DEFINES) -DMAST $< -o $@

$(SEROBJS_NSTX) : $(OBJDIR)/nstx/%.o : $(MAFOT_DIR)/src/%.cxx
	$(CXX) -c $(CFLAGS) -MMD $(INCLUDE) $(DEFINES) -DNSTX $< -o $@

$(SEROBJS_CMOD) : $(OBJDIR)/cmod/%.o : $(MAFOT_DIR)/src/%.cxx
	$(CXX) -c $(CFLAGS) -MMD $(INCLUDE) $(DEFINES) -DCMOD $< -o $@

