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
# Tous les .o et .mod vont dans OBJDIR (pas la racine) -- seul
# l'executable final reste ici (cf. consigne : la fonction doit rester
# dans le dossier apres make). -module $(OBJDIR) sert les DEUX roles a la
# fois chez ifx : ecrire les .mod produits ET les chercher pour les USE,
# quel que soit le sous-dossier du .f90 source -- pas besoin de -I
# supplementaire.
OBJDIR = build
OPTF = -O2 -qopenmp -warn all -check bounds -module $(OBJDIR)

LIBMKL = -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
SCALAPACKLIBS = -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64
MPILIBS = -L$(IMPI_LIBDIR) -L$(MKLROOT)/lib $(SCALAPACKLIBS) $(LIBMKL) -lpthread -lm -ldl

# "all" DOIT rester la toute PREMIERE regle du fichier : make sans
# argument construit la premiere regle rencontree, pas forcement "all"
# (piege reellement rencontre en ecrivant ce Makefile -- `make` avait
# construit uniquement blr_kinds_mod.o).
all: example_blr_mpi_pivot_user

# Cree OBJDIR avant toute compilation -- prerequis "order-only" (apres le
# |, dans chaque regle .o ci-dessous) : garantit que le dossier existe
# sans jamais forcer de recompilation quand son horodatage change (ce
# qu'un prerequis normal ferait a chaque nouveau fichier cree dedans).
$(OBJDIR):
	mkdir -p $(OBJDIR)

# --- blr_lib/ : copie locale des modules partages du solveur, DANS
# L'ORDRE DE DEPENDANCE (chaque .o depend des .o des modules qu'il USE).
# mkl_interfaces_mod.o (interfaces explicites CGEMM/CTRSM/CGETRF/...,
# voir blr_lib/mkl_interfaces_mod.f90) est une dependance transverse,
# necessaire a quiconque appelle du BLAS/LAPACK/MKL sans passer par un
# module de plus haut niveau.
LIBOBJS = $(OBJDIR)/mkl_interfaces_mod.o $(OBJDIR)/blr_kinds_mod.o $(OBJDIR)/progress_bar_mod.o $(OBJDIR)/thread_control_mod.o \
          $(OBJDIR)/blr_arena_mod.o $(OBJDIR)/blr_compress_mod.o $(OBJDIR)/blr_ooc_mod.o \
          $(OBJDIR)/mpi_grid_mod.o $(OBJDIR)/mpi_lu_pivot_mod.o $(OBJDIR)/mpi_blr_lu_pivot_bcast_mod.o \
          $(OBJDIR)/bench_utils_mod.o $(OBJDIR)/mpi_lu_solve_2d_mod.o $(OBJDIR)/blr_matrix_mod.o \
          $(OBJDIR)/mpi_blr_matrix_mod.o

$(OBJDIR)/mkl_interfaces_mod.o: blr_lib/mkl_interfaces_mod.f90 | $(OBJDIR)
	$(MPIFC) $(OPTF) -c blr_lib/mkl_interfaces_mod.f90 -o $@

$(OBJDIR)/blr_kinds_mod.o: blr_lib/blr_kinds_mod.f90 | $(OBJDIR)
	$(MPIFC) $(OPTF) -c blr_lib/blr_kinds_mod.f90 -o $@

$(OBJDIR)/progress_bar_mod.o: blr_lib/progress_bar_mod.f90 | $(OBJDIR)
	$(MPIFC) $(OPTF) -c blr_lib/progress_bar_mod.f90 -o $@

$(OBJDIR)/thread_control_mod.o: blr_lib/thread_control_mkl.f90 $(OBJDIR)/mkl_interfaces_mod.o | $(OBJDIR)
	$(MPIFC) $(OPTF) -c blr_lib/thread_control_mkl.f90 -o $@

$(OBJDIR)/blr_arena_mod.o: blr_lib/blr_arena_mod.f90 $(OBJDIR)/blr_kinds_mod.o | $(OBJDIR)
	$(MPIFC) $(OPTF) -c blr_lib/blr_arena_mod.f90 -o $@

$(OBJDIR)/blr_compress_mod.o: blr_lib/blr_compress_mod.f90 $(OBJDIR)/blr_kinds_mod.o $(OBJDIR)/blr_arena_mod.o $(OBJDIR)/mkl_interfaces_mod.o | $(OBJDIR)
	$(MPIFC) $(OPTF) -c blr_lib/blr_compress_mod.f90 -o $@

$(OBJDIR)/blr_ooc_mod.o: blr_lib/blr_ooc_mod.f90 $(OBJDIR)/blr_kinds_mod.o $(OBJDIR)/blr_compress_mod.o | $(OBJDIR)
	$(MPIFC) $(OPTF) -c blr_lib/blr_ooc_mod.f90 -o $@

$(OBJDIR)/mpi_grid_mod.o: blr_lib/mpi_grid_mod.f90 $(OBJDIR)/blr_kinds_mod.o | $(OBJDIR)
	$(MPIFC) $(OPTF) -c blr_lib/mpi_grid_mod.f90 -o $@

$(OBJDIR)/mpi_lu_pivot_mod.o: blr_lib/mpi_lu_pivot_mod.f90 $(OBJDIR)/mpi_grid_mod.o $(OBJDIR)/progress_bar_mod.o $(OBJDIR)/mkl_interfaces_mod.o | $(OBJDIR)
	$(MPIFC) $(OPTF) -c blr_lib/mpi_lu_pivot_mod.f90 -o $@

$(OBJDIR)/mpi_blr_lu_pivot_bcast_mod.o: blr_lib/mpi_blr_lu_pivot_bcast_mod.f90 $(OBJDIR)/mpi_grid_mod.o $(OBJDIR)/blr_compress_mod.o $(OBJDIR)/blr_arena_mod.o $(OBJDIR)/mpi_lu_pivot_mod.o $(OBJDIR)/mkl_interfaces_mod.o | $(OBJDIR)
	$(MPIFC) $(OPTF) -c blr_lib/mpi_blr_lu_pivot_bcast_mod.f90 -o $@

$(OBJDIR)/bench_utils_mod.o: blr_lib/bench_utils_mod.f90 | $(OBJDIR)
	$(MPIFC) $(OPTF) -c blr_lib/bench_utils_mod.f90 -o $@

$(OBJDIR)/mpi_lu_solve_2d_mod.o: blr_lib/mpi_lu_solve_2d_mod.f90 $(OBJDIR)/mpi_grid_mod.o $(OBJDIR)/mpi_lu_pivot_mod.o $(OBJDIR)/progress_bar_mod.o $(OBJDIR)/mkl_interfaces_mod.o | $(OBJDIR)
	$(MPIFC) $(OPTF) -c blr_lib/mpi_lu_solve_2d_mod.f90 -o $@

$(OBJDIR)/blr_matrix_mod.o: blr_lib/blr_matrix_mod.f90 $(OBJDIR)/blr_compress_mod.o $(OBJDIR)/blr_arena_mod.o $(OBJDIR)/progress_bar_mod.o $(OBJDIR)/thread_control_mod.o $(OBJDIR)/mkl_interfaces_mod.o | $(OBJDIR)
	$(MPIFC) $(OPTF) -c blr_lib/blr_matrix_mod.f90 -o $@

$(OBJDIR)/mpi_blr_matrix_mod.o: blr_lib/mpi_blr_matrix_mod.f90 $(OBJDIR)/mpi_grid_mod.o $(OBJDIR)/blr_compress_mod.o $(OBJDIR)/blr_arena_mod.o $(OBJDIR)/blr_ooc_mod.o $(OBJDIR)/blr_matrix_mod.o $(OBJDIR)/mpi_blr_lu_pivot_bcast_mod.o $(OBJDIR)/progress_bar_mod.o $(OBJDIR)/bench_utils_mod.o $(OBJDIR)/mkl_interfaces_mod.o | $(OBJDIR)
	$(MPIFC) $(OPTF) -c blr_lib/mpi_blr_matrix_mod.f90 -o $@

# --- scalapack/ : PCMATGEN (copie fidele, jamais modifiee) + les 4
# routines de generateur "a saut" qui lui manquent (reimplementees).
RANDGENOBJS = $(OBJDIR)/pcmatgen.o $(OBJDIR)/scalapack_randgen.o

$(OBJDIR)/pcmatgen.o: scalapack/pcmatgen.f | $(OBJDIR)
	$(MPIFC) $(OPTF) -c scalapack/pcmatgen.f -o $@

$(OBJDIR)/scalapack_randgen.o: scalapack/scalapack_randgen.f90 | $(OBJDIR)
	$(MPIFC) $(OPTF) -c scalapack/scalapack_randgen.f90 -o $@

# --- Modules propres a cette version utilisateur.
LOCALOBJS = $(OBJDIR)/cli_options_mod.o $(OBJDIR)/checkpoint_io_mod.o $(OBJDIR)/profiling_mod.o $(OBJDIR)/residual_mod.o $(OBJDIR)/user_matrix_mod.o

$(OBJDIR)/cli_options_mod.o: cli/cli_options_mod.f90 | $(OBJDIR)
	$(MPIFC) $(OPTF) -c cli/cli_options_mod.f90 -o $@

$(OBJDIR)/checkpoint_io_mod.o: checkpoint/checkpoint_io_mod.f90 $(OBJDIR)/mpi_blr_matrix_mod.o | $(OBJDIR)
	$(MPIFC) $(OPTF) -c checkpoint/checkpoint_io_mod.f90 -o $@

$(OBJDIR)/profiling_mod.o: profiling/profiling_mod.f90 $(OBJDIR)/mpi_blr_matrix_mod.o $(OBJDIR)/bench_utils_mod.o | $(OBJDIR)
	$(MPIFC) $(OPTF) -c profiling/profiling_mod.f90 -o $@

$(OBJDIR)/residual_mod.o: residual/residual_mod.f90 $(OBJDIR)/mpi_blr_matrix_mod.o $(OBJDIR)/mpi_lu_solve_2d_mod.o $(OBJDIR)/blr_matrix_mod.o | $(OBJDIR)
	$(MPIFC) $(OPTF) -c residual/residual_mod.f90 -o $@

$(OBJDIR)/user_matrix_mod.o: matrix/user_matrix_mod.f90 $(OBJDIR)/mpi_blr_matrix_mod.o $(RANDGENOBJS) | $(OBJDIR)
	$(MPIFC) $(OPTF) -c matrix/user_matrix_mod.f90 -o $@

# --- Programme principal + executable -- l'executable reste A LA RACINE
# de user_case_dev/ (seul son .o va dans OBJDIR avec tout le reste).
$(OBJDIR)/example_blr_mpi_pivot_user.o: example_blr_mpi_pivot_user.f90 $(LOCALOBJS) | $(OBJDIR)
	$(MPIFC) $(OPTF) -c example_blr_mpi_pivot_user.f90 -o $@

example_blr_mpi_pivot_user: $(LIBOBJS) $(RANDGENOBJS) $(LOCALOBJS) $(OBJDIR)/example_blr_mpi_pivot_user.o
	$(MPIFC) -o $@ $(OPTF) $(LIBOBJS) $(RANDGENOBJS) $(LOCALOBJS) $(OBJDIR)/example_blr_mpi_pivot_user.o $(MPILIBS)

clean:
	$(RM) -r $(OBJDIR)
	$(RM) *__genmod.f90 example_blr_mpi_pivot_user
