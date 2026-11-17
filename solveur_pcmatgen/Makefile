#===============================================================================
# solveur_pcmatgen -- solveur BLR + MPI (LU ou LDLt), A et B generees par
# PCMATGEN. Dossier autonome : `make` suffit, sans rien compiler ailleurs.
#
# Prerequis : ifx + MKL + Intel MPI (source /opt/intel/oneapi/setvars.sh).
#
# Options de compilation (combinables) :
#   make OPENMP=0   sans OpenMP, MKL sequentiel (recommande en MPI pur :
#                   memoire mesuree plus faible et plus fiable)
#   make CHECK=0    sans verification des bornes de tableaux (plus rapide,
#                   pour les calculs de production)
#   make FAST=1     -O3 + jeu d'instructions de la machine (FAST_ARCH,
#                   defaut -xHost ; sur un cluster, compiler sur un noeud de
#                   calcul ou fixer FAST_ARCH=-xCORE-AVX512 par exemple)
# Chaque combinaison a son propre dossier d'objets et son propre executable
# (suffixes _noomp, _nocheck, _fast), make ne detectant pas un changement
# d'options.
#===============================================================================
MPIFC ?= mpiifx
MKLROOT ?= /opt/intel/oneapi/mkl/latest

OPENMP ?= 1
ifeq ($(OPENMP),0)
  OMPF   =
  MKLTHR = -lmkl_sequential
  SUFFIX = _noomp
else
  OMPF   = -qopenmp
  MKLTHR = -lmkl_intel_thread
  SUFFIX =
endif

CHECK ?= 1
ifeq ($(CHECK),0)
  CHECKF =
  SUFFIX := $(SUFFIX)_nocheck
else
  CHECKF = -check bounds
endif

FAST ?= 0
FAST_ARCH ?= -xHost
ifeq ($(FAST),1)
  OPTLEVEL = -O3 $(FAST_ARCH)
  SUFFIX := $(SUFFIX)_fast
else
  OPTLEVEL = -O2
endif

OBJDIR = build$(SUFFIX)
EXE    = solveur_blr$(SUFFIX)
FFLAGS = $(OPTLEVEL) $(OMPF) $(CHECKF) -module $(OBJDIR)
# ScaLAPACK/BLACS de MKL : fournissent NUMROC, ICEIL et PXERBLA appeles par
# PCMATGEN (aucun contexte BLACS n'est cree).
LIBS   = -L$(MKLROOT)/lib -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 \
         -lmkl_intel_lp64 $(MKLTHR) -lmkl_core -lpthread -lm -ldl

# Ordre de compilation = ordre des dependances entre modules.
MODS = lib/blr_kinds_mod lib/progress_bar_mod lib/bench_utils_mod lib/blr_compress_mod \
       lib/blr_matrix_mod lib/mpi_grid_mod lib/mpi_lu_solve_2d_mod lib/mpi_blr_lu_pivot_bcast_mod \
       lib/mpi_blr_matrix_mod lib/mpi_blr_ldlt_mod \
       app/matrix_gen_mod app/options_mod app/checkpoint_io_mod app/profiling_mod app/residual_mod
OBJS = $(addprefix $(OBJDIR)/,$(notdir $(addsuffix .o,$(MODS)))) \
       $(OBJDIR)/pcmatgen.o $(OBJDIR)/scalapack_randgen.o

all: $(EXE)

$(EXE): $(OBJS) $(OBJDIR)/solveur_blr.o
	$(MPIFC) $(FFLAGS) -o $@ $^ $(LIBS)

# Chaque module depend du precedent dans MODS : compilation sequentielle
# dans le bon ordre, sans liste de dependances a maintenir a la main.
define MODRULE
$(OBJDIR)/$(notdir $(1)).o: $(1).f90 $(2) | $(OBJDIR)
	$$(MPIFC) $$(FFLAGS) -c $$< -o $$@
endef
PREV :=
$(foreach m,$(MODS),$(eval $(call MODRULE,$(m),$(PREV)))$(eval PREV := $(OBJDIR)/$(notdir $(m)).o))

$(OBJDIR)/pcmatgen.o: scalapack/pcmatgen.f | $(OBJDIR)
	$(MPIFC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/scalapack_randgen.o: scalapack/scalapack_randgen.f90 | $(OBJDIR)
	$(MPIFC) $(FFLAGS) -c $< -o $@

$(OBJDIR)/solveur_blr.o: app/solveur_blr.f90 $(OBJS) | $(OBJDIR)
	$(MPIFC) $(FFLAGS) -c $< -o $@

$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -rf build build_* solveur_blr solveur_blr_*

.PHONY: all clean
