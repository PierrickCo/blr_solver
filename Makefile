# Solveur BLR "from scratch" -- ligne VIVANTE (5 solveurs de reference :
# A/B/MPI/C/BLR+MPI, voir README.md). Le code mis de cote (LDLt, ancienne
# lignee MPI a stockage dense) vit dans legacy/, avec son propre Makefile
# (voir legacy/README.md) -- construire ce Makefile-ci EN PREMIER, legacy/
# reutilise directement ses .o/.mod.
#
# Prerequis : source /opt/intel/oneapi/setvars.sh (voir INSTALL.md pour le
# chemin alternatif gcc+OpenBLAS+OpenMPI).

MKLROOT = /opt/intel/oneapi/mkl/2026.1
IMPI_LIBDIR = /opt/intel/oneapi/mpi/2021.18/lib

FC   = ifx
OPTF = -O2 -qopenmp -warn all -check bounds

# -L$(IMPI_LIBDIR) EN PREMIER dans MPILIBS : sans ca, le linker peut
# resoudre "-lmpi" vers un OpenMPI systeme non lie (/usr/lib/x86_64-
# linux-gnu/libmpi.so.40, installe par un paquet apt sans rapport) au lieu
# de celui d'Intel MPI (libmpi.so.12) -- plantage immediat au premier
# appel MPI. Piege rencontre et corrige lors de l'exploration STRUMPACK
# (voir la memoire du projet), applique ici de maniere preventive.
MPIFC = mpiifx

SRC  = src
TEST = test

INCS = -I$(SRC)

LIBMKL = -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
LIBS   = $(LIBMKL) -lpthread -lm -ldl
MPILIBS = -L$(IMPI_LIBDIR) $(LIBS)

all: test_block_lu test_blr_compress test_mpi_grid test_mpi_lu_pivot test_openmp_dense test_mpi_lu_solve_2d test_blr_matrix_progress test_blr_workspace test_mpi_blr_matrix test_mpi_blr_matrix_pivot bench_seq bench_mpi bench_blr_pipeline

blr_kinds_mod.o: $(SRC)/blr_kinds_mod.f90
	$(FC) $(OPTF) -c $(SRC)/blr_kinds_mod.f90

# --- Reference C : LU par blocs dense, sans compression (NB=N degenere en
# un seul CGETRF monolithique, voir bench_seq.f90) -----------------------
test_block_lu: blr_kinds_mod.o block_lu_mod.o test_block_lu.o
	$(FC) -o $@ $(OPTF) blr_kinds_mod.o block_lu_mod.o test_block_lu.o $(LIBS)

block_lu_mod.o: $(SRC)/block_lu_mod.f90 blr_kinds_mod.o
	$(FC) $(OPTF) -c $(SRC)/block_lu_mod.f90

test_block_lu.o: $(TEST)/test_block_lu.f90 block_lu_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/test_block_lu.f90

# --- Bibliotheque d'algebre par blocs full/low-rank (base de A/B et de
# BLR+MPI) ----------------------------------------------------------------
test_blr_compress: blr_kinds_mod.o blr_compress_mod.o test_blr_compress.o
	$(FC) -o $@ $(OPTF) blr_kinds_mod.o blr_compress_mod.o test_blr_compress.o $(LIBS)

blr_compress_mod.o: $(SRC)/blr_compress_mod.f90 blr_kinds_mod.o
	$(FC) $(OPTF) -c $(SRC)/blr_compress_mod.f90

test_blr_compress.o: $(TEST)/test_blr_compress.f90 blr_compress_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/test_blr_compress.f90

# --- Grille de processus + distribution bloc-cyclique 2D (base MPI) -----
test_mpi_grid: blr_kinds_mod.o mpi_grid_mod.o test_mpi_grid.o
	$(MPIFC) -o $@ $(OPTF) blr_kinds_mod.o mpi_grid_mod.o test_mpi_grid.o $(MPILIBS)

mpi_grid_mod.o: $(SRC)/mpi_grid_mod.f90 blr_kinds_mod.o
	$(MPIFC) $(OPTF) -c $(SRC)/mpi_grid_mod.f90

test_mpi_grid.o: $(TEST)/test_mpi_grid.f90 mpi_grid_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/test_mpi_grid.f90

# --- Reference MPI (base) : factorisation LU dense distribuee, AVEC pivot
# Reutilise block_lu_mod.o (compile via FC=ifx) tel quel a l'edition de
# liens : ifx et mpiifx sont le meme compilateur sous-jacent (mpiifx
# ajoute seulement les includes/libs MPI), .o/.mod directement compatibles.
test_mpi_lu_pivot: blr_kinds_mod.o mpi_grid_mod.o mpi_lu_pivot_mod.o block_lu_mod.o test_mpi_lu_pivot.o
	$(MPIFC) -o $@ $(OPTF) blr_kinds_mod.o mpi_grid_mod.o mpi_lu_pivot_mod.o block_lu_mod.o test_mpi_lu_pivot.o $(MPILIBS)

mpi_lu_pivot_mod.o: $(SRC)/mpi_lu_pivot_mod.f90 mpi_grid_mod.o
	$(MPIFC) $(OPTF) -c $(SRC)/mpi_lu_pivot_mod.f90

test_mpi_lu_pivot.o: $(TEST)/test_mpi_lu_pivot.f90 mpi_lu_pivot_mod.o block_lu_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/test_mpi_lu_pivot.f90

# --- mpi_blr_lu_pivot_bcast_mod.o : garde UNIQUEMENT pour BCAST_BLR_BLOCK,
# reutilisee par mpi_blr_matrix_mod.f90 (diffusion de panneaux DEJA
# compresses) -- sa propre routine de factorisation top-level (proposition
# 3 historique) n'est plus testee ici, voir legacy/ pour son test dedie.
mpi_blr_lu_pivot_bcast_mod.o: $(SRC)/mpi_blr_lu_pivot_bcast_mod.f90 mpi_grid_mod.o blr_compress_mod.o mpi_lu_pivot_mod.o
	$(MPIFC) $(OPTF) -c $(SRC)/mpi_blr_lu_pivot_bcast_mod.f90

# --- Case "dense + OpenMP" : pas de directive OMP a nous (un seul gros
# CGEMM par panneau, deja thread par MKL) -- ce test MESURE que
# l'acceleration existe quand meme plutot que de l'affirmer.
test_openmp_dense: blr_kinds_mod.o block_lu_mod.o test_openmp_dense.o
	$(FC) -o $@ $(OPTF) blr_kinds_mod.o block_lu_mod.o test_openmp_dense.o $(LIBS)

test_openmp_dense.o: $(TEST)/test_openmp_dense.f90 block_lu_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/test_openmp_dense.f90

# --- Benchmarks de reference (A/B/MPI/C, N=10000/NRHS=100 par defaut) ---
bench_utils_mod.o: $(SRC)/bench_utils_mod.f90
	$(FC) $(OPTF) -c $(SRC)/bench_utils_mod.f90

# reference C
bench_seq: blr_kinds_mod.o block_lu_mod.o bench_utils_mod.o bench_seq.o
	$(FC) -o $@ $(OPTF) blr_kinds_mod.o block_lu_mod.o bench_utils_mod.o bench_seq.o $(LIBS)

bench_seq.o: $(TEST)/bench_seq.f90 block_lu_mod.o bench_utils_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/bench_seq.f90

# reference MPI (dense+pivot + resolution 2D)
bench_mpi: blr_kinds_mod.o mpi_grid_mod.o mpi_lu_pivot_mod.o mpi_lu_solve_2d_mod.o block_lu_mod.o bench_utils_mod.o bench_mpi.o
	$(MPIFC) -o $@ $(OPTF) blr_kinds_mod.o mpi_grid_mod.o mpi_lu_pivot_mod.o mpi_lu_solve_2d_mod.o block_lu_mod.o bench_utils_mod.o bench_mpi.o $(MPILIBS)

bench_mpi.o: $(TEST)/bench_mpi.f90 mpi_lu_pivot_mod.o mpi_lu_solve_2d_mod.o block_lu_mod.o bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/bench_mpi.f90

# --- Resolution MPI distribuee, vraie distribution 2D du second membre --
test_mpi_lu_solve_2d: blr_kinds_mod.o mpi_grid_mod.o mpi_lu_pivot_mod.o mpi_lu_solve_2d_mod.o test_mpi_lu_solve_2d.o
	$(MPIFC) -o $@ $(OPTF) blr_kinds_mod.o mpi_grid_mod.o mpi_lu_pivot_mod.o mpi_lu_solve_2d_mod.o test_mpi_lu_solve_2d.o $(MPILIBS)

mpi_lu_solve_2d_mod.o: $(SRC)/mpi_lu_solve_2d_mod.f90 mpi_grid_mod.o mpi_lu_pivot_mod.o
	$(MPIFC) $(OPTF) -c $(SRC)/mpi_lu_solve_2d_mod.f90

test_mpi_lu_solve_2d.o: $(TEST)/test_mpi_lu_solve_2d.f90 mpi_lu_pivot_mod.o mpi_lu_solve_2d_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/test_mpi_lu_solve_2d.f90

# --- Solveur "BLR+MPI" (frontiere du projet) : stockage local en
# BLR_BLOCK_T (jamais dense), panneaux L/U diffuses deja compresses,
# resolution 2D complete -- variantes SANS et AVEC pivot confine au bloc
# diagonal --------------------------------------------------------------
test_mpi_blr_matrix: blr_kinds_mod.o mpi_grid_mod.o blr_compress_mod.o progress_bar_mod.o thread_control_mod.o blr_matrix_mod.o mpi_blr_lu_pivot_bcast_mod.o mpi_lu_pivot_mod.o mpi_lu_solve_2d_mod.o mpi_blr_matrix_mod.o bench_utils_mod.o test_mpi_blr_matrix.o
	$(MPIFC) -o $@ $(OPTF) blr_kinds_mod.o mpi_grid_mod.o blr_compress_mod.o progress_bar_mod.o thread_control_mod.o blr_matrix_mod.o mpi_blr_lu_pivot_bcast_mod.o mpi_lu_pivot_mod.o mpi_lu_solve_2d_mod.o mpi_blr_matrix_mod.o bench_utils_mod.o test_mpi_blr_matrix.o $(MPILIBS)

mpi_blr_matrix_mod.o: $(SRC)/mpi_blr_matrix_mod.f90 mpi_grid_mod.o blr_compress_mod.o blr_matrix_mod.o mpi_blr_lu_pivot_bcast_mod.o progress_bar_mod.o
	$(MPIFC) $(OPTF) -c $(SRC)/mpi_blr_matrix_mod.f90

test_mpi_blr_matrix.o: $(TEST)/test_mpi_blr_matrix.f90 mpi_blr_matrix_mod.o mpi_lu_solve_2d_mod.o bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/test_mpi_blr_matrix.f90

test_mpi_blr_matrix_pivot: blr_kinds_mod.o mpi_grid_mod.o blr_compress_mod.o progress_bar_mod.o thread_control_mod.o blr_matrix_mod.o mpi_blr_lu_pivot_bcast_mod.o mpi_lu_pivot_mod.o mpi_lu_solve_2d_mod.o mpi_blr_matrix_mod.o bench_utils_mod.o test_mpi_blr_matrix_pivot.o
	$(MPIFC) -o $@ $(OPTF) blr_kinds_mod.o mpi_grid_mod.o blr_compress_mod.o progress_bar_mod.o thread_control_mod.o blr_matrix_mod.o mpi_blr_lu_pivot_bcast_mod.o mpi_lu_pivot_mod.o mpi_lu_solve_2d_mod.o mpi_blr_matrix_mod.o bench_utils_mod.o test_mpi_blr_matrix_pivot.o $(MPILIBS)

test_mpi_blr_matrix_pivot.o: $(TEST)/test_mpi_blr_matrix_pivot.f90 mpi_blr_matrix_mod.o mpi_lu_solve_2d_mod.o bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/test_mpi_blr_matrix_pivot.f90

# --- References A/B : pipeline "matrice compressee de bout en bout" -----
progress_bar_mod.o: $(SRC)/progress_bar_mod.f90
	$(FC) $(OPTF) -c $(SRC)/progress_bar_mod.f90

bench_blr_pipeline: blr_kinds_mod.o blr_compress_mod.o progress_bar_mod.o thread_control_mod.o blr_matrix_mod.o bench_utils_mod.o bench_blr_pipeline.o
	$(FC) -o $@ $(OPTF) blr_kinds_mod.o blr_compress_mod.o progress_bar_mod.o thread_control_mod.o blr_matrix_mod.o bench_utils_mod.o bench_blr_pipeline.o $(LIBS)

thread_control_mod.o: $(SRC)/thread_control_mkl.f90
	$(FC) $(OPTF) -c $(SRC)/thread_control_mkl.f90 -o thread_control_mod.o

blr_matrix_mod.o: $(SRC)/blr_matrix_mod.f90 blr_compress_mod.o progress_bar_mod.o thread_control_mod.o
	$(FC) $(OPTF) -c $(SRC)/blr_matrix_mod.f90

bench_blr_pipeline.o: $(TEST)/bench_blr_pipeline.f90 blr_matrix_mod.o bench_utils_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/bench_blr_pipeline.f90

test_blr_matrix_progress: blr_kinds_mod.o blr_compress_mod.o progress_bar_mod.o thread_control_mod.o blr_matrix_mod.o test_blr_matrix_progress.o
	$(FC) -o $@ $(OPTF) blr_kinds_mod.o blr_compress_mod.o progress_bar_mod.o thread_control_mod.o blr_matrix_mod.o test_blr_matrix_progress.o $(LIBS)

test_blr_matrix_progress.o: $(TEST)/test_blr_matrix_progress.f90 blr_matrix_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/test_blr_matrix_progress.f90

test_blr_workspace: blr_kinds_mod.o blr_compress_mod.o test_blr_workspace.o
	$(FC) -o $@ $(OPTF) blr_kinds_mod.o blr_compress_mod.o test_blr_workspace.o $(LIBS)

test_blr_workspace.o: $(TEST)/test_blr_workspace.f90 blr_compress_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/test_blr_workspace.f90

# --- Exemples canoniques (examples/) : programmes simplifies, un par
# paradigme (voir examples/README.md), pour prise en main rapide -- PAS
# dans "all" (inutiles a run_all_tests.sh), a construire explicitement
# avec `make examples`. Reutilisent directement les .o/.mod deja
# construits ci-dessus, comme legacy/Makefile.
EX = examples

examples: example_blr_sequentiel example_dense_sequentiel example_mpi_dense example_blr_mpi_sans_pivot example_blr_mpi_pivot

example_blr_sequentiel: blr_kinds_mod.o blr_compress_mod.o progress_bar_mod.o thread_control_mod.o blr_matrix_mod.o $(EX)/example_blr_sequentiel.o
	$(FC) -o $@ $(OPTF) blr_kinds_mod.o blr_compress_mod.o progress_bar_mod.o thread_control_mod.o blr_matrix_mod.o $(EX)/example_blr_sequentiel.o $(LIBS)

$(EX)/example_blr_sequentiel.o: $(EX)/example_blr_sequentiel.f90 blr_matrix_mod.o
	$(FC) $(OPTF) $(INCS) -c $(EX)/example_blr_sequentiel.f90 -o $(EX)/example_blr_sequentiel.o

example_dense_sequentiel: blr_kinds_mod.o block_lu_mod.o $(EX)/example_dense_sequentiel.o
	$(FC) -o $@ $(OPTF) blr_kinds_mod.o block_lu_mod.o $(EX)/example_dense_sequentiel.o $(LIBS)

$(EX)/example_dense_sequentiel.o: $(EX)/example_dense_sequentiel.f90 block_lu_mod.o
	$(FC) $(OPTF) $(INCS) -c $(EX)/example_dense_sequentiel.f90 -o $(EX)/example_dense_sequentiel.o

example_mpi_dense: blr_kinds_mod.o mpi_grid_mod.o mpi_lu_pivot_mod.o mpi_lu_solve_2d_mod.o $(EX)/example_mpi_dense.o
	$(MPIFC) -o $@ $(OPTF) blr_kinds_mod.o mpi_grid_mod.o mpi_lu_pivot_mod.o mpi_lu_solve_2d_mod.o $(EX)/example_mpi_dense.o $(MPILIBS)

$(EX)/example_mpi_dense.o: $(EX)/example_mpi_dense.f90 mpi_lu_pivot_mod.o mpi_lu_solve_2d_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(EX)/example_mpi_dense.f90 -o $(EX)/example_mpi_dense.o

example_blr_mpi_sans_pivot: blr_kinds_mod.o mpi_grid_mod.o blr_compress_mod.o progress_bar_mod.o thread_control_mod.o blr_matrix_mod.o mpi_blr_lu_pivot_bcast_mod.o mpi_lu_pivot_mod.o mpi_lu_solve_2d_mod.o mpi_blr_matrix_mod.o $(EX)/example_blr_mpi_sans_pivot.o
	$(MPIFC) -o $@ $(OPTF) blr_kinds_mod.o mpi_grid_mod.o blr_compress_mod.o progress_bar_mod.o thread_control_mod.o blr_matrix_mod.o mpi_blr_lu_pivot_bcast_mod.o mpi_lu_pivot_mod.o mpi_lu_solve_2d_mod.o mpi_blr_matrix_mod.o $(EX)/example_blr_mpi_sans_pivot.o $(MPILIBS)

$(EX)/example_blr_mpi_sans_pivot.o: $(EX)/example_blr_mpi_sans_pivot.f90 mpi_blr_matrix_mod.o mpi_lu_solve_2d_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(EX)/example_blr_mpi_sans_pivot.f90 -o $(EX)/example_blr_mpi_sans_pivot.o

example_blr_mpi_pivot: blr_kinds_mod.o mpi_grid_mod.o blr_compress_mod.o progress_bar_mod.o thread_control_mod.o blr_matrix_mod.o mpi_blr_lu_pivot_bcast_mod.o mpi_lu_pivot_mod.o mpi_lu_solve_2d_mod.o mpi_blr_matrix_mod.o $(EX)/example_blr_mpi_pivot.o
	$(MPIFC) -o $@ $(OPTF) blr_kinds_mod.o mpi_grid_mod.o blr_compress_mod.o progress_bar_mod.o thread_control_mod.o blr_matrix_mod.o mpi_blr_lu_pivot_bcast_mod.o mpi_lu_pivot_mod.o mpi_lu_solve_2d_mod.o mpi_blr_matrix_mod.o $(EX)/example_blr_mpi_pivot.o $(MPILIBS)

$(EX)/example_blr_mpi_pivot.o: $(EX)/example_blr_mpi_pivot.f90 mpi_blr_matrix_mod.o mpi_lu_solve_2d_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(EX)/example_blr_mpi_pivot.f90 -o $(EX)/example_blr_mpi_pivot.o

clean:
	$(RM) *.o *.mod test_block_lu test_blr_compress test_mpi_grid test_mpi_lu_pivot test_openmp_dense test_mpi_lu_solve_2d bench_seq bench_mpi bench_blr_pipeline test_blr_matrix_progress test_blr_workspace test_mpi_blr_matrix test_mpi_blr_matrix_pivot
	$(RM) $(EX)/*.o example_blr_sequentiel example_dense_sequentiel example_mpi_dense example_blr_mpi_sans_pivot example_blr_mpi_pivot
