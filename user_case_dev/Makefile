# Build ENTIEREMENT AUTONOME -- ne depend plus du Makefile racine de
# blr_solver/ ni de ../lib : blr_lib/ contient sa PROPRE copie des
# modules partages du solveur (vendorisee depuis ../src/, voir
# blr_lib/README si present) ; scalapack/ contient sa propre copie de
# PCMATGEN (depuis ../../scalapack_ooc_ref/) et les 4 routines qui lui
# manquent, reimplementees. Seule dependance externe restante : MKL +
# Intel MPI (bibliotheques systeme, hors du projet).
#
# `make clean && make` fonctionne ici SEUL, sans etape prealable ailleurs.
#
# Prerequis : source /opt/intel/oneapi/setvars.sh (voir ../INSTALL.md).

MKLROOT = /opt/intel/oneapi/mkl/2026.1
IMPI_LIBDIR = /opt/intel/oneapi/mpi/2021.18/lib

MPIFC = mpiifx
# -module . : tous les .mod sont ecrits ICI (racine de user_case_dev),
# quel que soit le sous-dossier du .f90 source -- pas besoin de -I
# supplementaire, tout le monde cherche/ecrit au meme endroit.
OPTF = -O2 -qopenmp -warn all -check bounds -module .

LIBMKL = -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
SCALAPACKLIBS = -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64
MPILIBS = -L$(IMPI_LIBDIR) -L$(MKLROOT)/lib $(SCALAPACKLIBS) $(LIBMKL) -lpthread -lm -ldl

# "all" DOIT rester la toute PREMIERE regle du fichier : make sans
# argument construit la premiere regle rencontree, pas forcement "all"
# (piege reellement rencontre en ecrivant ce Makefile -- `make` avait
# construit uniquement blr_kinds_mod.o).
all: example_blr_mpi_pivot_user

# --- blr_lib/ : copie locale des modules partages du solveur, DANS
# L'ORDRE DE DEPENDANCE (chaque .o depend des .o des modules qu'il USE).
# mkl_interfaces_mod.o (interfaces explicites CGEMM/CTRSM/CGETRF/...,
# voir blr_lib/mkl_interfaces_mod.f90) est une dependance transverse,
# necessaire a quiconque appelle du BLAS/LAPACK/MKL sans passer par un
# module de plus haut niveau.
LIBOBJS = mkl_interfaces_mod.o blr_kinds_mod.o progress_bar_mod.o thread_control_mod.o \
          blr_arena_mod.o blr_compress_mod.o blr_ooc_mod.o \
          mpi_grid_mod.o mpi_lu_pivot_mod.o mpi_blr_lu_pivot_bcast_mod.o \
          bench_utils_mod.o mpi_lu_solve_2d_mod.o blr_matrix_mod.o \
          mpi_blr_matrix_mod.o

mkl_interfaces_mod.o: blr_lib/mkl_interfaces_mod.f90
	$(MPIFC) $(OPTF) -c blr_lib/mkl_interfaces_mod.f90 -o mkl_interfaces_mod.o

blr_kinds_mod.o: blr_lib/blr_kinds_mod.f90
	$(MPIFC) $(OPTF) -c blr_lib/blr_kinds_mod.f90 -o blr_kinds_mod.o

progress_bar_mod.o: blr_lib/progress_bar_mod.f90
	$(MPIFC) $(OPTF) -c blr_lib/progress_bar_mod.f90 -o progress_bar_mod.o

thread_control_mod.o: blr_lib/thread_control_mkl.f90 mkl_interfaces_mod.o
	$(MPIFC) $(OPTF) -c blr_lib/thread_control_mkl.f90 -o thread_control_mod.o

blr_arena_mod.o: blr_lib/blr_arena_mod.f90 blr_kinds_mod.o
	$(MPIFC) $(OPTF) -c blr_lib/blr_arena_mod.f90 -o blr_arena_mod.o

blr_compress_mod.o: blr_lib/blr_compress_mod.f90 blr_kinds_mod.o blr_arena_mod.o mkl_interfaces_mod.o
	$(MPIFC) $(OPTF) -c blr_lib/blr_compress_mod.f90 -o blr_compress_mod.o

blr_ooc_mod.o: blr_lib/blr_ooc_mod.f90 blr_kinds_mod.o blr_compress_mod.o
	$(MPIFC) $(OPTF) -c blr_lib/blr_ooc_mod.f90 -o blr_ooc_mod.o

mpi_grid_mod.o: blr_lib/mpi_grid_mod.f90 blr_kinds_mod.o
	$(MPIFC) $(OPTF) -c blr_lib/mpi_grid_mod.f90 -o mpi_grid_mod.o

mpi_lu_pivot_mod.o: blr_lib/mpi_lu_pivot_mod.f90 mpi_grid_mod.o progress_bar_mod.o mkl_interfaces_mod.o
	$(MPIFC) $(OPTF) -c blr_lib/mpi_lu_pivot_mod.f90 -o mpi_lu_pivot_mod.o

mpi_blr_lu_pivot_bcast_mod.o: blr_lib/mpi_blr_lu_pivot_bcast_mod.f90 mpi_grid_mod.o blr_compress_mod.o blr_arena_mod.o mpi_lu_pivot_mod.o mkl_interfaces_mod.o
	$(MPIFC) $(OPTF) -c blr_lib/mpi_blr_lu_pivot_bcast_mod.f90 -o mpi_blr_lu_pivot_bcast_mod.o

bench_utils_mod.o: blr_lib/bench_utils_mod.f90
	$(MPIFC) $(OPTF) -c blr_lib/bench_utils_mod.f90 -o bench_utils_mod.o

mpi_lu_solve_2d_mod.o: blr_lib/mpi_lu_solve_2d_mod.f90 mpi_grid_mod.o mpi_lu_pivot_mod.o progress_bar_mod.o mkl_interfaces_mod.o
	$(MPIFC) $(OPTF) -c blr_lib/mpi_lu_solve_2d_mod.f90 -o mpi_lu_solve_2d_mod.o

blr_matrix_mod.o: blr_lib/blr_matrix_mod.f90 blr_compress_mod.o blr_arena_mod.o progress_bar_mod.o thread_control_mod.o mkl_interfaces_mod.o
	$(MPIFC) $(OPTF) -c blr_lib/blr_matrix_mod.f90 -o blr_matrix_mod.o

mpi_blr_matrix_mod.o: blr_lib/mpi_blr_matrix_mod.f90 mpi_grid_mod.o blr_compress_mod.o blr_arena_mod.o blr_ooc_mod.o blr_matrix_mod.o mpi_blr_lu_pivot_bcast_mod.o progress_bar_mod.o bench_utils_mod.o mkl_interfaces_mod.o
	$(MPIFC) $(OPTF) -c blr_lib/mpi_blr_matrix_mod.f90 -o mpi_blr_matrix_mod.o

# --- scalapack/ : PCMATGEN (copie fidele, jamais modifiee) + les 4
# routines de generateur "a saut" qui lui manquent (reimplementees).
RANDGENOBJS = pcmatgen.o scalapack_randgen.o

pcmatgen.o: scalapack/pcmatgen.f
	$(MPIFC) $(OPTF) -c scalapack/pcmatgen.f -o pcmatgen.o

scalapack_randgen.o: scalapack/scalapack_randgen.f90
	$(MPIFC) $(OPTF) -c scalapack/scalapack_randgen.f90 -o scalapack_randgen.o

# --- Modules propres a cette version utilisateur.
LOCALOBJS = cli_options_mod.o checkpoint_io_mod.o profiling_mod.o residual_mod.o user_matrix_mod.o

cli_options_mod.o: cli/cli_options_mod.f90
	$(MPIFC) $(OPTF) -c cli/cli_options_mod.f90 -o cli_options_mod.o

checkpoint_io_mod.o: checkpoint/checkpoint_io_mod.f90 mpi_blr_matrix_mod.o
	$(MPIFC) $(OPTF) -c checkpoint/checkpoint_io_mod.f90 -o checkpoint_io_mod.o

profiling_mod.o: profiling/profiling_mod.f90 mpi_blr_matrix_mod.o bench_utils_mod.o
	$(MPIFC) $(OPTF) -c profiling/profiling_mod.f90 -o profiling_mod.o

residual_mod.o: residual/residual_mod.f90 mpi_blr_matrix_mod.o mpi_lu_solve_2d_mod.o blr_matrix_mod.o
	$(MPIFC) $(OPTF) -c residual/residual_mod.f90 -o residual_mod.o

user_matrix_mod.o: matrix/user_matrix_mod.f90 mpi_blr_matrix_mod.o $(RANDGENOBJS)
	$(MPIFC) $(OPTF) -c matrix/user_matrix_mod.f90 -o user_matrix_mod.o

# --- Programme principal + executable.
example_blr_mpi_pivot_user.o: example_blr_mpi_pivot_user.f90 $(LOCALOBJS)
	$(MPIFC) $(OPTF) -c example_blr_mpi_pivot_user.f90 -o example_blr_mpi_pivot_user.o

example_blr_mpi_pivot_user: $(LIBOBJS) $(RANDGENOBJS) $(LOCALOBJS) example_blr_mpi_pivot_user.o
	$(MPIFC) -o $@ $(OPTF) $(LIBOBJS) $(RANDGENOBJS) $(LOCALOBJS) example_blr_mpi_pivot_user.o $(MPILIBS)

clean:
	$(RM) *.o *.mod *__genmod.f90 example_blr_mpi_pivot_user
