# Solveur BLR "from scratch" -- ligne VIVANTE (5 solveurs de reference :
# A/B/MPI/C/BLR+MPI, voir README.md). Le code mis de cote (LDLt, ancienne
# lignee MPI a stockage dense) vit dans legacy/, avec son propre Makefile
# (voir legacy/README.md) -- construire ce Makefile-ci EN PREMIER, legacy/
# reutilise directement ses .o/.mod (desormais dans lib/, voir MAIN dans
# legacy/Makefile).
#
# Prerequis : source /opt/intel/oneapi/setvars.sh (voir INSTALL.md pour le
# chemin alternatif gcc+OpenBLAS+OpenMPI).
#
# RANGEMENT (demande explicite) : tous les .o/.mod generes vont dans lib/,
# tous les executables (tests/benchmarks a la racine ET exemples de
# examples/) vont dans bin/ -- plus rien ne pollue la racine du projet.
# Les deux repertoires sont crees automatiquement (voir plus bas) des que
# make est invoque, sans etape manuelle prealable.

MKLROOT = /opt/intel/oneapi/mkl/2026.1
IMPI_LIBDIR = /opt/intel/oneapi/mpi/2021.18/lib

FC   = ifx
# -module $(LIBDIR) : sortie ET recherche des .mod dans lib/ (equivalent
# du -J de gfortran) -- inclus dans OPTF (reutilise pour compiler ET
# lier) pour s'appliquer a TOUS les appels sans dupliquer le flag partout.
# OPENMP=0 (`make OPENMP=0`) : compile SANS OpenMP et lie MKL SEQUENTIEL
# (pas de libiomp5). Recommande pour les runs MPI purs (1 thread par
# processus) : avec ifx, des que libiomp5 est lie, tous les ALLOCATE passent
# par l'allocateur du runtime OpenMP d'Intel, qui ne rend jamais la memoire
# liberee -- mesure : RSS de factorisation -45 % sur une matrice peu
# compressible, meme temps (voir CLUSTER.md). Sorties dans des repertoires
# separes (suffixe _noomp) pour ne jamais melanger les deux variantes.
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

# CHECK=0 (`make CHECK=0`) : compile SANS verification des bornes de
# tableaux. Utile au developpement, couteuse a l'execution -- mesure
# (N=40000, NB=500) : assemblage x2.4 (LDLt) a x3.5 (LU), factorisation LU
# x1.7, resolution x1.1 a x1.4. Recommande pour les calculs de production ;
# combinable avec OPENMP=0. Sorties dans des repertoires separes (suffixe
# _nocheck) : make ne detecte pas un changement d'options.
CHECK ?= 1
ifeq ($(CHECK),0)
  CHECKF =
  SUFFIX := $(SUFFIX)_nocheck
else
  CHECKF = -check bounds
endif

OPTF = -O2 $(OMPF) -warn all $(CHECKF) -module $(LIBDIR)

# -L$(IMPI_LIBDIR) EN PREMIER dans MPILIBS : sans ca, le linker peut
# resoudre "-lmpi" vers un OpenMPI systeme non lie (/usr/lib/x86_64-
# linux-gnu/libmpi.so.40, installe par un paquet apt sans rapport) au lieu
# de celui d'Intel MPI (libmpi.so.12) -- plantage immediat au premier
# appel MPI. Piege rencontre et corrige lors de l'exploration STRUMPACK
# (voir la memoire du projet), applique ici de maniere preventive.
MPIFC = mpiifx

SRC  = src
TEST = test
EX   = examples
USERDIR = user_case
BINDIR = bin$(SUFFIX)
LIBDIR = lib$(SUFFIX)

# Cree les deux repertoires de sortie AU MOMENT DU PARSING (avant que la
# moindre regle ne s'execute) -- evite d'ajouter un prerequis d'ordre
# (| $(BINDIR)) a chacune des ~65 regles de ce fichier.
$(shell mkdir -p $(BINDIR) $(LIBDIR))

INCS = -I$(SRC)

LIBMKL = -L$(MKLROOT)/lib -lmkl_intel_lp64 $(MKLTHR) -lmkl_core
LIBS   = $(LIBMKL) -lpthread -lm -ldl
MPILIBS = -L$(IMPI_LIBDIR) $(LIBS)

# test_openmp_dense teste OpenMP lui-meme : exclu en OPENMP=0.
ifeq ($(OPENMP),0)
  OMP_TESTS =
else
  OMP_TESTS = $(BINDIR)/test_openmp_dense
endif

all: $(BINDIR)/test_block_lu $(BINDIR)/test_blr_compress $(BINDIR)/test_mpi_grid $(BINDIR)/test_mpi_lu_pivot $(OMP_TESTS) $(BINDIR)/test_mpi_lu_solve_2d $(BINDIR)/test_blr_matrix_progress $(BINDIR)/test_blr_workspace $(BINDIR)/test_mpi_blr_matrix $(BINDIR)/test_mpi_blr_matrix_pivot $(BINDIR)/test_mpi_blr_matrix_pivot_arena $(BINDIR)/test_mpi_blr_matrix_pivot_ooc $(BINDIR)/bench_seq $(BINDIR)/bench_mpi $(BINDIR)/bench_blr_pipeline $(BINDIR)/bench_blr_mpi_pivot_compressibility $(BINDIR)/bench_blr_mpi_pivot_ratio $(BINDIR)/bench_blr_mpi_pivot_arena $(BINDIR)/bench_blr_mpi_pivot_ooc_compare $(BINDIR)/bench_blr_mpi_pivot_ooc_massive_rhs $(BINDIR)/bench_blr_mpi_pivot_ooc_geometric

$(LIBDIR)/blr_kinds_mod.o: $(SRC)/blr_kinds_mod.f90
	$(FC) $(OPTF) -c $(SRC)/blr_kinds_mod.f90 -o $(LIBDIR)/blr_kinds_mod.o

# --- Reference C : LU par blocs dense, sans compression (NB=N degenere en
# un seul CGETRF monolithique, voir bench_seq.f90) -----------------------
$(BINDIR)/test_block_lu: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/block_lu_mod.o $(LIBDIR)/test_block_lu.o
	$(FC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/block_lu_mod.o $(LIBDIR)/test_block_lu.o $(LIBS)

$(LIBDIR)/block_lu_mod.o: $(SRC)/block_lu_mod.f90 $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/progress_bar_mod.o
	$(FC) $(OPTF) -c $(SRC)/block_lu_mod.f90 -o $(LIBDIR)/block_lu_mod.o

$(LIBDIR)/test_block_lu.o: $(TEST)/test_block_lu.f90 $(LIBDIR)/block_lu_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/test_block_lu.f90 -o $(LIBDIR)/test_block_lu.o

# --- Bibliotheque d'algebre par blocs full/low-rank (base de A/B et de
# BLR+MPI) ----------------------------------------------------------------
$(BINDIR)/test_blr_compress: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/test_blr_compress.o
	$(FC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/test_blr_compress.o $(LIBS)

$(BINDIR)/bench_rsvd_compress: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/bench_rsvd_compress.o
	$(FC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/bench_rsvd_compress.o $(LIBS)

$(LIBDIR)/bench_rsvd_compress.o: $(TEST)/bench_rsvd_compress.f90 $(LIBDIR)/blr_compress_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/bench_rsvd_compress.f90 -o $(LIBDIR)/bench_rsvd_compress.o

$(LIBDIR)/blr_arena_mod.o: $(SRC)/blr_arena_mod.f90 $(LIBDIR)/blr_kinds_mod.o
	$(FC) $(OPTF) -c $(SRC)/blr_arena_mod.f90 -o $(LIBDIR)/blr_arena_mod.o

$(LIBDIR)/blr_compress_mod.o: $(SRC)/blr_compress_mod.f90 $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/blr_arena_mod.o
	$(FC) $(OPTF) -c $(SRC)/blr_compress_mod.f90 -o $(LIBDIR)/blr_compress_mod.o

$(LIBDIR)/blr_ooc_mod.o: $(SRC)/blr_ooc_mod.f90 $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/blr_compress_mod.o
	$(FC) $(OPTF) -c $(SRC)/blr_ooc_mod.f90 -o $(LIBDIR)/blr_ooc_mod.o

$(LIBDIR)/test_blr_compress.o: $(TEST)/test_blr_compress.f90 $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/test_blr_compress.f90 -o $(LIBDIR)/test_blr_compress.o

# --- Grille de processus + distribution bloc-cyclique 2D (base MPI) -----
$(BINDIR)/test_mpi_grid: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/test_mpi_grid.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/test_mpi_grid.o $(MPILIBS)

$(LIBDIR)/mpi_grid_mod.o: $(SRC)/mpi_grid_mod.f90 $(LIBDIR)/blr_kinds_mod.o
	$(MPIFC) $(OPTF) -c $(SRC)/mpi_grid_mod.f90 -o $(LIBDIR)/mpi_grid_mod.o

$(LIBDIR)/test_mpi_grid.o: $(TEST)/test_mpi_grid.f90 $(LIBDIR)/mpi_grid_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/test_mpi_grid.f90 -o $(LIBDIR)/test_mpi_grid.o

# --- Reference MPI (base) : factorisation LU dense distribuee, AVEC pivot
# Reutilise block_lu_mod.o (compile via FC=ifx) tel quel a l'edition de
# liens : ifx et mpiifx sont le meme compilateur sous-jacent (mpiifx
# ajoute seulement les includes/libs MPI), .o/.mod directement compatibles.
$(BINDIR)/test_mpi_lu_pivot: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/block_lu_mod.o $(LIBDIR)/test_mpi_lu_pivot.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/block_lu_mod.o $(LIBDIR)/test_mpi_lu_pivot.o $(MPILIBS)

$(LIBDIR)/mpi_lu_pivot_mod.o: $(SRC)/mpi_lu_pivot_mod.f90 $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/progress_bar_mod.o
	$(MPIFC) $(OPTF) -c $(SRC)/mpi_lu_pivot_mod.f90 -o $(LIBDIR)/mpi_lu_pivot_mod.o

$(LIBDIR)/test_mpi_lu_pivot.o: $(TEST)/test_mpi_lu_pivot.f90 $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/block_lu_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/test_mpi_lu_pivot.f90 -o $(LIBDIR)/test_mpi_lu_pivot.o

# --- mpi_blr_lu_pivot_bcast_mod.o : garde UNIQUEMENT pour BCAST_BLR_BLOCK,
# reutilisee par mpi_blr_matrix_mod.f90 (diffusion de panneaux DEJA
# compresses) -- sa propre routine de factorisation top-level (proposition
# 3 historique) n'est plus testee ici, voir legacy/ pour son test dedie.
$(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o: $(SRC)/mpi_blr_lu_pivot_bcast_mod.f90 $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o
	$(MPIFC) $(OPTF) -c $(SRC)/mpi_blr_lu_pivot_bcast_mod.f90 -o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o

# --- Case "dense + OpenMP" : pas de directive OMP a nous (un seul gros
# CGEMM par panneau, deja thread par MKL) -- ce test MESURE que
# l'acceleration existe quand meme plutot que de l'affirmer.
$(BINDIR)/test_openmp_dense: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/block_lu_mod.o $(LIBDIR)/test_openmp_dense.o
	$(FC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/block_lu_mod.o $(LIBDIR)/test_openmp_dense.o $(LIBS)

$(LIBDIR)/test_openmp_dense.o: $(TEST)/test_openmp_dense.f90 $(LIBDIR)/block_lu_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/test_openmp_dense.f90 -o $(LIBDIR)/test_openmp_dense.o

# --- Benchmarks de reference (A/B/MPI/C, N=10000/NRHS=100 par defaut) ---
$(LIBDIR)/bench_utils_mod.o: $(SRC)/bench_utils_mod.f90
	$(FC) $(OPTF) -c $(SRC)/bench_utils_mod.f90 -o $(LIBDIR)/bench_utils_mod.o

# reference C
$(BINDIR)/bench_seq: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/block_lu_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/bench_seq.o
	$(FC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/block_lu_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/bench_seq.o $(LIBS)

$(LIBDIR)/bench_seq.o: $(TEST)/bench_seq.f90 $(LIBDIR)/block_lu_mod.o $(LIBDIR)/bench_utils_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/bench_seq.f90 -o $(LIBDIR)/bench_seq.o

# reference MPI (dense+pivot + resolution 2D)
$(BINDIR)/bench_mpi: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/block_lu_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/bench_mpi.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/block_lu_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/bench_mpi.o $(MPILIBS)

$(LIBDIR)/bench_mpi.o: $(TEST)/bench_mpi.f90 $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/block_lu_mod.o $(LIBDIR)/bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/bench_mpi.f90 -o $(LIBDIR)/bench_mpi.o

# --- Resolution MPI distribuee, vraie distribution 2D du second membre --
$(BINDIR)/test_mpi_lu_solve_2d: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/test_mpi_lu_solve_2d.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/test_mpi_lu_solve_2d.o $(MPILIBS)

$(LIBDIR)/mpi_lu_solve_2d_mod.o: $(SRC)/mpi_lu_solve_2d_mod.f90 $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/progress_bar_mod.o
	$(MPIFC) $(OPTF) -c $(SRC)/mpi_lu_solve_2d_mod.f90 -o $(LIBDIR)/mpi_lu_solve_2d_mod.o

$(LIBDIR)/test_mpi_lu_solve_2d.o: $(TEST)/test_mpi_lu_solve_2d.f90 $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/test_mpi_lu_solve_2d.f90 -o $(LIBDIR)/test_mpi_lu_solve_2d.o

# --- Solveur "BLR+MPI" (frontiere du projet) : stockage local en
# BLR_BLOCK_T (jamais dense), panneaux L/U diffuses deja compresses,
# resolution 2D complete -- variantes SANS et AVEC pivot confine au bloc
# diagonal --------------------------------------------------------------
$(BINDIR)/test_mpi_blr_matrix: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/test_mpi_blr_matrix.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/test_mpi_blr_matrix.o $(MPILIBS)

$(LIBDIR)/mpi_blr_matrix_mod.o: $(SRC)/mpi_blr_matrix_mod.f90 $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/bench_utils_mod.o
	$(MPIFC) $(OPTF) -c $(SRC)/mpi_blr_matrix_mod.f90 -o $(LIBDIR)/mpi_blr_matrix_mod.o

$(LIBDIR)/test_mpi_blr_matrix.o: $(TEST)/test_mpi_blr_matrix.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/test_mpi_blr_matrix.f90 -o $(LIBDIR)/test_mpi_blr_matrix.o

$(BINDIR)/test_mpi_blr_matrix_pivot: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/test_mpi_blr_matrix_pivot.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/test_mpi_blr_matrix_pivot.o $(MPILIBS)

$(LIBDIR)/test_mpi_blr_matrix_pivot.o: $(TEST)/test_mpi_blr_matrix_pivot.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/test_mpi_blr_matrix_pivot.f90 -o $(LIBDIR)/test_mpi_blr_matrix_pivot.o

$(BINDIR)/test_mpi_blr_matrix_pivot_arena: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/test_mpi_blr_matrix_pivot_arena.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/test_mpi_blr_matrix_pivot_arena.o $(MPILIBS)

$(LIBDIR)/test_mpi_blr_matrix_pivot_arena.o: $(TEST)/test_mpi_blr_matrix_pivot_arena.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/test_mpi_blr_matrix_pivot_arena.f90 -o $(LIBDIR)/test_mpi_blr_matrix_pivot_arena.o

$(BINDIR)/test_mpi_blr_matrix_pivot_ooc: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/test_mpi_blr_matrix_pivot_ooc.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/test_mpi_blr_matrix_pivot_ooc.o $(MPILIBS)

$(LIBDIR)/test_mpi_blr_matrix_pivot_ooc.o: $(TEST)/test_mpi_blr_matrix_pivot_ooc.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/test_mpi_blr_matrix_pivot_ooc.f90 -o $(LIBDIR)/test_mpi_blr_matrix_pivot_ooc.o

$(BINDIR)/bench_blr_mpi_pivot_ooc_compare: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/bench_blr_mpi_pivot_ooc_compare.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/bench_blr_mpi_pivot_ooc_compare.o $(MPILIBS)

$(LIBDIR)/bench_blr_mpi_pivot_ooc_compare.o: $(TEST)/bench_blr_mpi_pivot_ooc_compare.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/bench_blr_mpi_pivot_ooc_compare.f90 -o $(LIBDIR)/bench_blr_mpi_pivot_ooc_compare.o

$(BINDIR)/bench_blr_mpi_pivot_ooc_geometric: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/bench_blr_mpi_pivot_ooc_geometric.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/bench_blr_mpi_pivot_ooc_geometric.o $(MPILIBS)

$(LIBDIR)/bench_blr_mpi_pivot_ooc_geometric.o: $(TEST)/bench_blr_mpi_pivot_ooc_geometric.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/bench_blr_mpi_pivot_ooc_geometric.f90 -o $(LIBDIR)/bench_blr_mpi_pivot_ooc_geometric.o

$(BINDIR)/bench_blr_mpi_pivot_ooc_massive_rhs: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/bench_blr_mpi_pivot_ooc_massive_rhs.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/bench_blr_mpi_pivot_ooc_massive_rhs.o $(MPILIBS)

$(LIBDIR)/bench_blr_mpi_pivot_ooc_massive_rhs.o: $(TEST)/bench_blr_mpi_pivot_ooc_massive_rhs.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/bench_blr_mpi_pivot_ooc_massive_rhs.f90 -o $(LIBDIR)/bench_blr_mpi_pivot_ooc_massive_rhs.o

$(BINDIR)/bench_blr_mpi_pivot_compressibility: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/bench_blr_mpi_pivot_compressibility.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/bench_blr_mpi_pivot_compressibility.o $(MPILIBS)

$(LIBDIR)/bench_blr_mpi_pivot_compressibility.o: $(TEST)/bench_blr_mpi_pivot_compressibility.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/bench_blr_mpi_pivot_compressibility.f90 -o $(LIBDIR)/bench_blr_mpi_pivot_compressibility.o

$(BINDIR)/bench_blr_mpi_pivot_ratio: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/bench_blr_mpi_pivot_ratio.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/bench_blr_mpi_pivot_ratio.o $(MPILIBS)

$(LIBDIR)/bench_blr_mpi_pivot_ratio.o: $(TEST)/bench_blr_mpi_pivot_ratio.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/bench_blr_mpi_pivot_ratio.f90 -o $(LIBDIR)/bench_blr_mpi_pivot_ratio.o

$(BINDIR)/bench_blr_mpi_pivot_arena: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/bench_blr_mpi_pivot_arena.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/bench_blr_mpi_pivot_arena.o $(MPILIBS)

$(LIBDIR)/bench_blr_mpi_pivot_arena.o: $(TEST)/bench_blr_mpi_pivot_arena.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/bench_blr_mpi_pivot_arena.f90 -o $(LIBDIR)/bench_blr_mpi_pivot_arena.o

# --- References A/B : pipeline "matrice compressee de bout en bout" -----
$(LIBDIR)/progress_bar_mod.o: $(SRC)/progress_bar_mod.f90
	$(FC) $(OPTF) -c $(SRC)/progress_bar_mod.f90 -o $(LIBDIR)/progress_bar_mod.o

$(BINDIR)/bench_blr_pipeline: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/bench_blr_pipeline.o
	$(FC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/bench_blr_pipeline.o $(LIBS)

$(LIBDIR)/thread_control_mod.o: $(SRC)/thread_control_mkl.f90
	$(FC) $(OPTF) -c $(SRC)/thread_control_mkl.f90 -o $(LIBDIR)/thread_control_mod.o

$(LIBDIR)/blr_matrix_mod.o: $(SRC)/blr_matrix_mod.f90 $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o
	$(FC) $(OPTF) -c $(SRC)/blr_matrix_mod.f90 -o $(LIBDIR)/blr_matrix_mod.o

$(LIBDIR)/bench_blr_pipeline.o: $(TEST)/bench_blr_pipeline.f90 $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/bench_blr_pipeline.f90 -o $(LIBDIR)/bench_blr_pipeline.o

$(BINDIR)/test_blr_matrix_progress: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/test_blr_matrix_progress.o
	$(FC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/test_blr_matrix_progress.o $(LIBS)

$(LIBDIR)/test_blr_matrix_progress.o: $(TEST)/test_blr_matrix_progress.f90 $(LIBDIR)/blr_matrix_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/test_blr_matrix_progress.f90 -o $(LIBDIR)/test_blr_matrix_progress.o

$(BINDIR)/test_blr_workspace: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/test_blr_workspace.o
	$(FC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/test_blr_workspace.o $(LIBS)

$(LIBDIR)/test_blr_workspace.o: $(TEST)/test_blr_workspace.f90 $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/test_blr_workspace.f90 -o $(LIBDIR)/test_blr_workspace.o

# --- Exemples canoniques (examples/) : programmes simplifies, un par
# paradigme (voir examples/README.md), pour prise en main rapide -- PAS
# dans "all" (inutiles a run_all_tests.sh), a construire explicitement
# avec `make examples`. Reutilisent directement les .o/.mod deja
# construits ci-dessus (dans lib/), comme legacy/Makefile. Les .o propres
# a chaque exemple restent dans examples/ (deja a l'ecart de la racine) ;
# seuls leurs executables rejoignent bin/.
examples: $(BINDIR)/example_blr_sequentiel $(BINDIR)/example_dense_sequentiel $(BINDIR)/example_mpi_dense $(BINDIR)/example_blr_mpi_sans_pivot $(BINDIR)/example_blr_mpi_pivot $(BINDIR)/example_blr_mpi_pivot_arena $(BINDIR)/example_blr_mpi_pivot_fixedcap $(BINDIR)/example_blr_mpi_pivot_ooc $(BINDIR)/example_user_matrix $(BINDIR)/example_user_matrix_ooc $(BINDIR)/example_blr_mpi_pivot_user

$(BINDIR)/example_blr_sequentiel: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(EX)/example_blr_sequentiel.o
	$(FC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(EX)/example_blr_sequentiel.o $(LIBS)

$(EX)/example_blr_sequentiel.o: $(EX)/example_blr_sequentiel.f90 $(LIBDIR)/blr_matrix_mod.o
	$(FC) $(OPTF) $(INCS) -c $(EX)/example_blr_sequentiel.f90 -o $(EX)/example_blr_sequentiel.o

$(BINDIR)/example_dense_sequentiel: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/block_lu_mod.o $(EX)/example_dense_sequentiel.o
	$(FC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/block_lu_mod.o $(EX)/example_dense_sequentiel.o $(LIBS)

$(EX)/example_dense_sequentiel.o: $(EX)/example_dense_sequentiel.f90 $(LIBDIR)/block_lu_mod.o
	$(FC) $(OPTF) $(INCS) -c $(EX)/example_dense_sequentiel.f90 -o $(EX)/example_dense_sequentiel.o

$(BINDIR)/example_mpi_dense: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(EX)/example_mpi_dense.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(EX)/example_mpi_dense.o $(MPILIBS)

$(EX)/example_mpi_dense.o: $(EX)/example_mpi_dense.f90 $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(EX)/example_mpi_dense.f90 -o $(EX)/example_mpi_dense.o

$(BINDIR)/example_blr_mpi_sans_pivot: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(EX)/example_blr_mpi_sans_pivot.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(EX)/example_blr_mpi_sans_pivot.o $(MPILIBS)

$(EX)/example_blr_mpi_sans_pivot.o: $(EX)/example_blr_mpi_sans_pivot.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(EX)/example_blr_mpi_sans_pivot.f90 -o $(EX)/example_blr_mpi_sans_pivot.o

$(BINDIR)/example_blr_mpi_pivot: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(EX)/example_blr_mpi_pivot.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(EX)/example_blr_mpi_pivot.o $(MPILIBS)

$(EX)/example_blr_mpi_pivot.o: $(EX)/example_blr_mpi_pivot.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(EX)/example_blr_mpi_pivot.f90 -o $(EX)/example_blr_mpi_pivot.o

$(BINDIR)/example_user_matrix: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(EX)/example_user_matrix.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(EX)/example_user_matrix.o $(MPILIBS)

$(EX)/example_user_matrix.o: $(EX)/example_user_matrix.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(EX)/example_user_matrix.f90 -o $(EX)/example_user_matrix.o

$(BINDIR)/example_blr_mpi_pivot_arena: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(EX)/example_blr_mpi_pivot_arena.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(EX)/example_blr_mpi_pivot_arena.o $(MPILIBS)

$(EX)/example_blr_mpi_pivot_arena.o: $(EX)/example_blr_mpi_pivot_arena.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(EX)/example_blr_mpi_pivot_arena.f90 -o $(EX)/example_blr_mpi_pivot_arena.o

$(BINDIR)/example_blr_mpi_pivot_fixedcap: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(EX)/example_blr_mpi_pivot_fixedcap.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(EX)/example_blr_mpi_pivot_fixedcap.o $(MPILIBS)

$(EX)/example_blr_mpi_pivot_fixedcap.o: $(EX)/example_blr_mpi_pivot_fixedcap.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(EX)/example_blr_mpi_pivot_fixedcap.f90 -o $(EX)/example_blr_mpi_pivot_fixedcap.o

$(BINDIR)/example_blr_mpi_pivot_ooc: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(EX)/example_blr_mpi_pivot_ooc.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(EX)/example_blr_mpi_pivot_ooc.o $(MPILIBS)

$(EX)/example_blr_mpi_pivot_ooc.o: $(EX)/example_blr_mpi_pivot_ooc.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(EX)/example_blr_mpi_pivot_ooc.f90 -o $(EX)/example_blr_mpi_pivot_ooc.o

$(BINDIR)/example_user_matrix_ooc: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(EX)/example_user_matrix_ooc.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(EX)/example_user_matrix_ooc.o $(MPILIBS)

$(EX)/example_user_matrix_ooc.o: $(EX)/example_user_matrix_ooc.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(EX)/example_user_matrix_ooc.f90 -o $(EX)/example_user_matrix_ooc.o

# --- Version utilisateur (user_case/, voir user_case/tools/README.md) :
# meme structure que les exemples ci-dessus, mais assemblage/second
# membre decoupes dans des modules dedies (cli_options_mod/checkpoint_
# io_mod/profiling_mod/residual_mod/user_matrix_mod, ce dernier a
# personnaliser) plutot qu'un seul fichier -- leurs .o/.mod rejoignent
# lib/ comme le reste de la bibliotheque, seul le .o du programme
# principal reste dans user_case/ (meme convention que les EX/*.o
# ci-dessus).
$(LIBDIR)/cli_options_mod.o: $(USERDIR)/cli_options_mod.f90
	$(MPIFC) $(OPTF) -c $(USERDIR)/cli_options_mod.f90 -o $(LIBDIR)/cli_options_mod.o

$(LIBDIR)/checkpoint_io_mod.o: $(USERDIR)/checkpoint_io_mod.f90 $(LIBDIR)/mpi_blr_matrix_mod.o
	$(MPIFC) $(OPTF) -c $(USERDIR)/checkpoint_io_mod.f90 -o $(LIBDIR)/checkpoint_io_mod.o

$(LIBDIR)/profiling_mod.o: $(USERDIR)/profiling_mod.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o
	$(MPIFC) $(OPTF) -c $(USERDIR)/profiling_mod.f90 -o $(LIBDIR)/profiling_mod.o

$(LIBDIR)/residual_mod.o: $(USERDIR)/residual_mod.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/blr_matrix_mod.o
	$(MPIFC) $(OPTF) -c $(USERDIR)/residual_mod.f90 -o $(LIBDIR)/residual_mod.o

$(LIBDIR)/user_matrix_mod.o: $(USERDIR)/user_matrix_mod.f90 $(LIBDIR)/mpi_blr_matrix_mod.o
	$(MPIFC) $(OPTF) -c $(USERDIR)/user_matrix_mod.f90 -o $(LIBDIR)/user_matrix_mod.o

$(BINDIR)/example_blr_mpi_pivot_user: $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/cli_options_mod.o $(LIBDIR)/checkpoint_io_mod.o $(LIBDIR)/profiling_mod.o $(LIBDIR)/residual_mod.o $(LIBDIR)/user_matrix_mod.o $(USERDIR)/example_blr_mpi_pivot_user.o
	$(MPIFC) -o $@ $(OPTF) $(LIBDIR)/blr_kinds_mod.o $(LIBDIR)/mpi_grid_mod.o $(LIBDIR)/blr_compress_mod.o $(LIBDIR)/blr_arena_mod.o $(LIBDIR)/blr_ooc_mod.o $(LIBDIR)/progress_bar_mod.o $(LIBDIR)/thread_control_mod.o $(LIBDIR)/blr_matrix_mod.o $(LIBDIR)/mpi_blr_lu_pivot_bcast_mod.o $(LIBDIR)/mpi_lu_pivot_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/cli_options_mod.o $(LIBDIR)/checkpoint_io_mod.o $(LIBDIR)/profiling_mod.o $(LIBDIR)/residual_mod.o $(LIBDIR)/user_matrix_mod.o $(USERDIR)/example_blr_mpi_pivot_user.o $(MPILIBS)

$(USERDIR)/example_blr_mpi_pivot_user.o: $(USERDIR)/example_blr_mpi_pivot_user.f90 $(LIBDIR)/mpi_blr_matrix_mod.o $(LIBDIR)/mpi_lu_solve_2d_mod.o $(LIBDIR)/bench_utils_mod.o $(LIBDIR)/cli_options_mod.o $(LIBDIR)/checkpoint_io_mod.o $(LIBDIR)/profiling_mod.o $(LIBDIR)/residual_mod.o $(LIBDIR)/user_matrix_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(USERDIR)/example_blr_mpi_pivot_user.f90 -o $(USERDIR)/example_blr_mpi_pivot_user.o

clean:
	$(RM) -r $(LIBDIR) $(BINDIR)
	$(RM) $(EX)/*.o
	$(RM) $(USERDIR)/*.o
