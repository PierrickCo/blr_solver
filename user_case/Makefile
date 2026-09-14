# Build de la version utilisateur (example_blr_mpi_pivot_user.f90 +
# ses modules cli_options_mod/checkpoint_io_mod/profiling_mod/
# residual_mod/user_matrix_mod) -- REUTILISE les .o/.mod DEJA CONSTRUITS
# par le Makefile principal de blr_solver/ (../lib/*.o, ../lib/*.mod).
# Executer d'abord `make` (ou `make examples`) a la racine de blr_solver
# AVANT ce Makefile-ci.
#
# Prerequis : source /opt/intel/oneapi/setvars.sh (voir ../INSTALL.md).

ROOT = ..
LIBDIR = $(ROOT)/lib

MKLROOT = /opt/intel/oneapi/mkl/2026.1
IMPI_LIBDIR = /opt/intel/oneapi/mpi/2021.18/lib

MPIFC = mpiifx
# -module . : les .mod PROPRES a ce dossier sont ecrits ICI, jamais dans
# ../lib -- -I$(LIBDIR) ajoute seulement ../lib a la recherche, pour lire
# les .mod deja construits par le Makefile racine (mpi_blr_matrix_mod.mod,
# etc.) sans les dupliquer ni les modifier.
OPTF = -O2 -qopenmp -warn all -check bounds -module . -I$(LIBDIR)

LIBMKL = -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
MPILIBS = -L$(IMPI_LIBDIR) $(LIBMKL) -lpthread -lm -ldl

# Meme jeu de .o que example_blr_mpi_pivot dans le Makefile racine (la
# version utilisateur reutilise exactement la meme bibliotheque BLR+MPI).
LIBOBJS = $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o \
          $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o \
          $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o \
          $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o \
          $(LIBDIR)/bench_utils_mod.o

# Modules propres a ce dossier, DANS L'ORDRE DE DEPENDANCE (chacun ne
# depend que de LIBOBJS ci-dessus, sauf residual_mod qui n'a besoin
# d'aucun autre module de ce dossier -- l'ordre choisi ici est simplement
# celui qui simplifie la lecture, pas une contrainte de compilation).
LOCALOBJS = cli_options_mod.o checkpoint_io_mod.o profiling_mod.o residual_mod.o user_matrix_mod.o

all: example_blr_mpi_pivot_user

example_blr_mpi_pivot_user: $(LOCALOBJS) example_blr_mpi_pivot_user.o
	$(MPIFC) -o $@ $(OPTF) $(LIBOBJS) $(LOCALOBJS) example_blr_mpi_pivot_user.o $(MPILIBS)

cli_options_mod.o: cli_options_mod.f90
	$(MPIFC) $(OPTF) -c cli_options_mod.f90 -o cli_options_mod.o

checkpoint_io_mod.o: checkpoint_io_mod.f90
	$(MPIFC) $(OPTF) -c checkpoint_io_mod.f90 -o checkpoint_io_mod.o

profiling_mod.o: profiling_mod.f90
	$(MPIFC) $(OPTF) -c profiling_mod.f90 -o profiling_mod.o

residual_mod.o: residual_mod.f90
	$(MPIFC) $(OPTF) -c residual_mod.f90 -o residual_mod.o

user_matrix_mod.o: user_matrix_mod.f90
	$(MPIFC) $(OPTF) -c user_matrix_mod.f90 -o user_matrix_mod.o

example_blr_mpi_pivot_user.o: example_blr_mpi_pivot_user.f90 $(LOCALOBJS)
	$(MPIFC) $(OPTF) -c example_blr_mpi_pivot_user.f90 -o example_blr_mpi_pivot_user.o

clean:
	$(RM) *.o *.mod example_blr_mpi_pivot_user
