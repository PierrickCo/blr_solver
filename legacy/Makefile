# Build des modules/tests LEGACY (mis de cote, voir legacy/README.md).
#
# PREREQUIS : compiler d'abord le projet PRINCIPAL (`make` depuis
# blr_solver/) -- ce Makefile reutilise directement les .o/.mod DEJA
# construits la-bas (../blr_kinds_mod.o, ../mpi_grid_mod.o, etc.) plutot
# que de les recompiler, pour ne jamais dupliquer le code source des
# modules encore vivants. Meme compilateur/options que le projet
# principal (doivent rester synchronises -- copies volontairement, pas
# generees, pour que ce fichier reste lisible seul).

MKLROOT = /opt/intel/oneapi/mkl/2026.1
IMPI_LIBDIR = /opt/intel/oneapi/mpi/2021.18/lib

FC    = ifx
MPIFC = mpiifx
OPTF  = -O2 -qopenmp -warn all -check bounds

MAIN  = ..
SRC   = src
TEST  = test

INCS = -I$(SRC) -I$(MAIN)

LIBMKL  = -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
LIBS    = $(LIBMKL) -lpthread -lm -ldl
MPILIBS = -L$(IMPI_LIBDIR) $(LIBS)

all: test_ldlt test_blr_ldlt test_mpi_matvec test_mpi_lu_nopiv test_mpi_lu_solve test_mpi_blr_lu_pivot test_mpi_blr_lu_nopiv test_blr_lu test_blr_lu_storage test_pivoting test_openmp test_mpi_blr_lu_pivot_bcast bench_seq_full bench_mpi_full

# --- LDLt (abandonne, voir blr_ldlt_mod.f90 pour le detail de l'echec) --
block_ldlt_mod.o: $(SRC)/block_ldlt_mod.f90 $(MAIN)/blr_kinds_mod.o
	$(FC) $(OPTF) $(INCS) -c $(SRC)/block_ldlt_mod.f90

test_ldlt: $(MAIN)/blr_kinds_mod.o block_ldlt_mod.o test_ldlt.o
	$(FC) -o $@ $(OPTF) $(MAIN)/blr_kinds_mod.o block_ldlt_mod.o test_ldlt.o $(LIBS)

test_ldlt.o: $(TEST)/test_ldlt.f90 block_ldlt_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/test_ldlt.f90

blr_ldlt_mod.o: $(SRC)/blr_ldlt_mod.f90 $(MAIN)/blr_compress_mod.o
	$(FC) $(OPTF) $(INCS) -c $(SRC)/blr_ldlt_mod.f90

test_blr_ldlt: $(MAIN)/blr_kinds_mod.o block_ldlt_mod.o $(MAIN)/blr_compress_mod.o blr_ldlt_mod.o test_blr_ldlt.o
	$(FC) -o $@ $(OPTF) $(MAIN)/blr_kinds_mod.o block_ldlt_mod.o $(MAIN)/blr_compress_mod.o blr_ldlt_mod.o test_blr_ldlt.o $(LIBS)

test_blr_ldlt.o: $(TEST)/test_blr_ldlt.f90 block_ldlt_mod.o blr_ldlt_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/test_blr_ldlt.f90

# --- Etape 2/5 + proposition 2 : BLR a stockage dense/compresse, ---------
# --- remplacees par le pipeline "jamais dense" de blr_matrix_mod.f90 ----
blr_lu_mod.o: $(SRC)/blr_lu_mod.f90 $(MAIN)/blr_compress_mod.o
	$(FC) $(OPTF) $(INCS) -c $(SRC)/blr_lu_mod.f90

test_blr_lu: $(MAIN)/blr_kinds_mod.o $(MAIN)/block_lu_mod.o $(MAIN)/blr_compress_mod.o blr_lu_mod.o test_blr_lu.o
	$(FC) -o $@ $(OPTF) $(MAIN)/blr_kinds_mod.o $(MAIN)/block_lu_mod.o $(MAIN)/blr_compress_mod.o blr_lu_mod.o test_blr_lu.o $(LIBS)

test_blr_lu.o: $(TEST)/test_blr_lu.f90 $(MAIN)/block_lu_mod.o blr_lu_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/test_blr_lu.f90

test_pivoting: $(MAIN)/blr_kinds_mod.o $(MAIN)/block_lu_mod.o $(MAIN)/blr_compress_mod.o blr_lu_mod.o test_pivoting.o
	$(FC) -o $@ $(OPTF) $(MAIN)/blr_kinds_mod.o $(MAIN)/block_lu_mod.o $(MAIN)/blr_compress_mod.o blr_lu_mod.o test_pivoting.o $(LIBS)

test_pivoting.o: $(TEST)/test_pivoting.f90 $(MAIN)/block_lu_mod.o blr_lu_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/test_pivoting.f90

test_openmp: $(MAIN)/blr_kinds_mod.o $(MAIN)/block_lu_mod.o $(MAIN)/blr_compress_mod.o blr_lu_mod.o test_openmp.o
	$(FC) -o $@ $(OPTF) $(MAIN)/blr_kinds_mod.o $(MAIN)/block_lu_mod.o $(MAIN)/blr_compress_mod.o blr_lu_mod.o test_openmp.o $(LIBS)

test_openmp.o: $(TEST)/test_openmp.f90 blr_lu_mod.o $(MAIN)/block_lu_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/test_openmp.f90

blr_lu_storage_mod.o: $(SRC)/blr_lu_storage_mod.f90 $(MAIN)/block_lu_mod.o blr_lu_mod.o
	$(FC) $(OPTF) $(INCS) -c $(SRC)/blr_lu_storage_mod.f90

test_blr_lu_storage: $(MAIN)/blr_kinds_mod.o $(MAIN)/block_lu_mod.o $(MAIN)/blr_compress_mod.o blr_lu_mod.o blr_lu_storage_mod.o test_blr_lu_storage.o
	$(FC) -o $@ $(OPTF) $(MAIN)/blr_kinds_mod.o $(MAIN)/block_lu_mod.o $(MAIN)/blr_compress_mod.o blr_lu_mod.o blr_lu_storage_mod.o test_blr_lu_storage.o $(LIBS)

test_blr_lu_storage.o: $(TEST)/test_blr_lu_storage.f90 blr_lu_storage_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/test_blr_lu_storage.f90

bench_seq_full: $(MAIN)/blr_kinds_mod.o $(MAIN)/block_lu_mod.o $(MAIN)/blr_compress_mod.o blr_lu_mod.o blr_lu_storage_mod.o $(MAIN)/bench_utils_mod.o bench_seq_full.o
	$(FC) -o $@ $(OPTF) $(MAIN)/blr_kinds_mod.o $(MAIN)/block_lu_mod.o $(MAIN)/blr_compress_mod.o blr_lu_mod.o blr_lu_storage_mod.o $(MAIN)/bench_utils_mod.o bench_seq_full.o $(LIBS)

bench_seq_full.o: $(TEST)/bench_seq_full.f90 blr_lu_storage_mod.o $(MAIN)/bench_utils_mod.o
	$(FC) $(OPTF) $(INCS) -c $(TEST)/bench_seq_full.f90

# --- 6b/6c : matvec + LU sans pivot distribues (etapes vers 6d) ---------
mpi_matvec_mod.o: $(SRC)/mpi_matvec_mod.f90 $(MAIN)/mpi_grid_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(SRC)/mpi_matvec_mod.f90

test_mpi_matvec: $(MAIN)/blr_kinds_mod.o $(MAIN)/mpi_grid_mod.o mpi_matvec_mod.o test_mpi_matvec.o
	$(MPIFC) -o $@ $(OPTF) $(MAIN)/blr_kinds_mod.o $(MAIN)/mpi_grid_mod.o mpi_matvec_mod.o test_mpi_matvec.o $(MPILIBS)

test_mpi_matvec.o: $(TEST)/test_mpi_matvec.f90 mpi_matvec_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/test_mpi_matvec.f90

mpi_lu_nopiv_mod.o: $(SRC)/mpi_lu_nopiv_mod.f90 $(MAIN)/mpi_grid_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(SRC)/mpi_lu_nopiv_mod.f90

test_mpi_lu_nopiv: $(MAIN)/blr_kinds_mod.o $(MAIN)/mpi_grid_mod.o mpi_lu_nopiv_mod.o test_mpi_lu_nopiv.o
	$(MPIFC) -o $@ $(OPTF) $(MAIN)/blr_kinds_mod.o $(MAIN)/mpi_grid_mod.o mpi_lu_nopiv_mod.o test_mpi_lu_nopiv.o $(MPILIBS)

test_mpi_lu_nopiv.o: $(TEST)/test_mpi_lu_nopiv.f90 mpi_lu_nopiv_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/test_mpi_lu_nopiv.f90

# --- 1D distributed solve (remplacee par la resolution 2D) --------------
mpi_lu_solve_mod.o: $(SRC)/mpi_lu_solve_mod.f90 $(MAIN)/mpi_grid_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(SRC)/mpi_lu_solve_mod.f90

test_mpi_lu_solve: $(MAIN)/blr_kinds_mod.o $(MAIN)/mpi_grid_mod.o $(MAIN)/mpi_lu_pivot_mod.o mpi_lu_solve_mod.o test_mpi_lu_solve.o
	$(MPIFC) -o $@ $(OPTF) $(MAIN)/blr_kinds_mod.o $(MAIN)/mpi_grid_mod.o $(MAIN)/mpi_lu_pivot_mod.o mpi_lu_solve_mod.o test_mpi_lu_solve.o $(MPILIBS)

test_mpi_lu_solve.o: $(TEST)/test_mpi_lu_solve.f90 $(MAIN)/mpi_lu_pivot_mod.o mpi_lu_solve_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/test_mpi_lu_solve.f90

bench_mpi_full: $(MAIN)/blr_kinds_mod.o $(MAIN)/mpi_grid_mod.o $(MAIN)/blr_compress_mod.o $(MAIN)/mpi_lu_pivot_mod.o mpi_blr_lu_pivot_mod.o $(MAIN)/mpi_blr_lu_pivot_bcast_mod.o mpi_lu_solve_mod.o $(MAIN)/mpi_lu_solve_2d_mod.o $(MAIN)/block_lu_mod.o $(MAIN)/bench_utils_mod.o bench_mpi_full.o
	$(MPIFC) -o $@ $(OPTF) $(MAIN)/blr_kinds_mod.o $(MAIN)/mpi_grid_mod.o $(MAIN)/blr_compress_mod.o $(MAIN)/mpi_lu_pivot_mod.o mpi_blr_lu_pivot_mod.o $(MAIN)/mpi_blr_lu_pivot_bcast_mod.o mpi_lu_solve_mod.o $(MAIN)/mpi_lu_solve_2d_mod.o $(MAIN)/block_lu_mod.o $(MAIN)/bench_utils_mod.o bench_mpi_full.o $(MPILIBS)

bench_mpi_full.o: $(TEST)/bench_mpi_full.f90 $(MAIN)/mpi_lu_pivot_mod.o mpi_blr_lu_pivot_mod.o $(MAIN)/mpi_blr_lu_pivot_bcast_mod.o mpi_lu_solve_mod.o $(MAIN)/mpi_lu_solve_2d_mod.o $(MAIN)/block_lu_mod.o $(MAIN)/bench_utils_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/bench_mpi_full.f90

# --- 6e + proposition 3 : BLR+pivot MPI, dense-puis-compression / --------
# --- compresse-puis-diffusion (BCAST_BLR_BLOCK reste vivant dans le  -----
# --- projet principal -- ne PAS le recompiler ici, reutiliser ../mpi_ ---
# --- blr_lu_pivot_bcast_mod.o tel quel) ----------------------------------
mpi_blr_lu_pivot_mod.o: $(SRC)/mpi_blr_lu_pivot_mod.f90 $(MAIN)/mpi_grid_mod.o $(MAIN)/blr_compress_mod.o $(MAIN)/mpi_lu_pivot_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(SRC)/mpi_blr_lu_pivot_mod.f90

test_mpi_blr_lu_pivot: $(MAIN)/blr_kinds_mod.o $(MAIN)/mpi_grid_mod.o $(MAIN)/blr_compress_mod.o $(MAIN)/mpi_lu_pivot_mod.o mpi_blr_lu_pivot_mod.o $(MAIN)/block_lu_mod.o test_mpi_blr_lu_pivot.o
	$(MPIFC) -o $@ $(OPTF) $(MAIN)/blr_kinds_mod.o $(MAIN)/mpi_grid_mod.o $(MAIN)/blr_compress_mod.o $(MAIN)/mpi_lu_pivot_mod.o mpi_blr_lu_pivot_mod.o $(MAIN)/block_lu_mod.o test_mpi_blr_lu_pivot.o $(MPILIBS)

test_mpi_blr_lu_pivot.o: $(TEST)/test_mpi_blr_lu_pivot.f90 mpi_blr_lu_pivot_mod.o $(MAIN)/block_lu_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/test_mpi_blr_lu_pivot.f90

test_mpi_blr_lu_pivot_bcast: $(MAIN)/blr_kinds_mod.o $(MAIN)/mpi_grid_mod.o $(MAIN)/blr_compress_mod.o $(MAIN)/mpi_lu_pivot_mod.o $(MAIN)/mpi_blr_lu_pivot_bcast_mod.o $(MAIN)/block_lu_mod.o test_mpi_blr_lu_pivot_bcast.o
	$(MPIFC) -o $@ $(OPTF) $(MAIN)/blr_kinds_mod.o $(MAIN)/mpi_grid_mod.o $(MAIN)/blr_compress_mod.o $(MAIN)/mpi_lu_pivot_mod.o $(MAIN)/mpi_blr_lu_pivot_bcast_mod.o $(MAIN)/block_lu_mod.o test_mpi_blr_lu_pivot_bcast.o $(MPILIBS)

test_mpi_blr_lu_pivot_bcast.o: $(TEST)/test_mpi_blr_lu_pivot_bcast.f90 $(MAIN)/mpi_blr_lu_pivot_bcast_mod.o $(MAIN)/block_lu_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/test_mpi_blr_lu_pivot_bcast.f90

# --- Combinaison BLR+MPI sans pivot, avant mpi_blr_matrix_mod.f90 --------
mpi_blr_lu_nopiv_mod.o: $(SRC)/mpi_blr_lu_nopiv_mod.f90 $(MAIN)/mpi_grid_mod.o $(MAIN)/blr_compress_mod.o mpi_lu_nopiv_mod.o $(MAIN)/mpi_blr_lu_pivot_bcast_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(SRC)/mpi_blr_lu_nopiv_mod.f90

test_mpi_blr_lu_nopiv: $(MAIN)/blr_kinds_mod.o $(MAIN)/mpi_grid_mod.o $(MAIN)/blr_compress_mod.o mpi_lu_nopiv_mod.o $(MAIN)/mpi_lu_pivot_mod.o $(MAIN)/mpi_blr_lu_pivot_bcast_mod.o mpi_blr_lu_nopiv_mod.o test_mpi_blr_lu_nopiv.o
	$(MPIFC) -o $@ $(OPTF) $(MAIN)/blr_kinds_mod.o $(MAIN)/mpi_grid_mod.o $(MAIN)/blr_compress_mod.o mpi_lu_nopiv_mod.o $(MAIN)/mpi_lu_pivot_mod.o $(MAIN)/mpi_blr_lu_pivot_bcast_mod.o mpi_blr_lu_nopiv_mod.o test_mpi_blr_lu_nopiv.o $(MPILIBS)

test_mpi_blr_lu_nopiv.o: $(TEST)/test_mpi_blr_lu_nopiv.f90 mpi_blr_lu_nopiv_mod.o
	$(MPIFC) $(OPTF) $(INCS) -c $(TEST)/test_mpi_blr_lu_nopiv.f90

clean:
	$(RM) *.o *.mod test_ldlt test_blr_ldlt test_mpi_matvec test_mpi_lu_nopiv test_mpi_lu_solve test_mpi_blr_lu_pivot test_mpi_blr_lu_pivot_bcast test_mpi_blr_lu_nopiv test_blr_lu test_blr_lu_storage test_pivoting test_openmp bench_seq_full bench_mpi_full
