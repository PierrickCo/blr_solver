# LDL^T BLR+MPI pour matrices complexes symetriques (voir README.md).
#
# PREREQUIS : construire d'abord le projet PRINCIPAL (`make` depuis
# blr_solver/) -- ce Makefile reutilise directement ses .o/.mod (../lib/)
# sans les recompiler, comme legacy/Makefile. Memes compilateur/options que
# le projet principal (copies volontairement, pas generees).
#
# Sorties : .o/.mod propres au LDL^T dans lib/, executables dans bin/.

MKLROOT = /opt/intel/oneapi/mkl/2026.1
IMPI_LIBDIR = /opt/intel/oneapi/mpi/2021.18/lib

MPIFC  = mpiifx
LIBDIR = lib
BINDIR = bin
MAIN   = ../lib
SRC    = src
TEST   = test
EX     = examples

# -module $(LIBDIR) : .mod du LDL^T ecrits ici ; -I$(MAIN) : .mod du projet
# principal relus la-bas.
OPTF = -O2 -qopenmp -warn all -check bounds -module $(LIBDIR) -I$(MAIN)

$(shell mkdir -p $(BINDIR) $(LIBDIR))

LIBMKL  = -L$(MKLROOT)/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
LIBS    = $(LIBMKL) -lpthread -lm -ldl
MPILIBS = -L$(IMPI_LIBDIR) $(LIBS)

# Objets du projet principal necessaires (meme liste que test_mpi_blr_matrix).
MAIN_OBJS = $(MAIN)/blr_kinds_mod.o $(MAIN)/mpi_grid_mod.o $(MAIN)/blr_compress_mod.o \
            $(MAIN)/blr_arena_mod.o $(MAIN)/blr_ooc_mod.o $(MAIN)/progress_bar_mod.o \
            $(MAIN)/thread_control_mod.o $(MAIN)/blr_matrix_mod.o $(MAIN)/mpi_blr_lu_pivot_bcast_mod.o \
            $(MAIN)/mpi_lu_pivot_mod.o $(MAIN)/mpi_lu_solve_2d_mod.o $(MAIN)/mpi_blr_matrix_mod.o \
            $(MAIN)/bench_utils_mod.o

LDLT_OBJS = $(LIBDIR)/mpi_blr_ldlt_mod.o $(LIBDIR)/mpi_blr_ldlt_arena_mod.o $(LIBDIR)/mpi_blr_ldlt_ooc_mod.o $(LIBDIR)/sym_testmat_mod.o

all: $(BINDIR)/test_mpi_blr_ldlt $(BINDIR)/test_mpi_blr_ldlt_arena $(BINDIR)/test_mpi_blr_ldlt_ooc $(BINDIR)/example_blr_mpi_ldlt $(BINDIR)/bench_ldlt_vs_lu

$(LIBDIR)/mpi_blr_ldlt_mod.o: $(SRC)/mpi_blr_ldlt_mod.f90 $(MAIN)/mpi_blr_matrix_mod.o
	$(MPIFC) $(OPTF) -c $(SRC)/mpi_blr_ldlt_mod.f90 -o $@

$(LIBDIR)/mpi_blr_ldlt_arena_mod.o: $(SRC)/mpi_blr_ldlt_arena_mod.f90 $(LIBDIR)/mpi_blr_ldlt_mod.o
	$(MPIFC) $(OPTF) -c $(SRC)/mpi_blr_ldlt_arena_mod.f90 -o $@

$(LIBDIR)/test_mpi_blr_ldlt_arena.o: $(TEST)/test_mpi_blr_ldlt_arena.f90 $(LDLT_OBJS)
	$(MPIFC) $(OPTF) -c $(TEST)/test_mpi_blr_ldlt_arena.f90 -o $@

$(BINDIR)/test_mpi_blr_ldlt_arena: $(MAIN_OBJS) $(LDLT_OBJS) $(LIBDIR)/test_mpi_blr_ldlt_arena.o
	$(MPIFC) -o $@ $(OPTF) $^ $(MPILIBS)

$(LIBDIR)/mpi_blr_ldlt_ooc_mod.o: $(SRC)/mpi_blr_ldlt_ooc_mod.f90 $(LIBDIR)/mpi_blr_ldlt_mod.o
	$(MPIFC) $(OPTF) -c $(SRC)/mpi_blr_ldlt_ooc_mod.f90 -o $@

$(LIBDIR)/test_mpi_blr_ldlt_ooc.o: $(TEST)/test_mpi_blr_ldlt_ooc.f90 $(LDLT_OBJS)
	$(MPIFC) $(OPTF) -c $(TEST)/test_mpi_blr_ldlt_ooc.f90 -o $@

$(BINDIR)/test_mpi_blr_ldlt_ooc: $(MAIN_OBJS) $(LDLT_OBJS) $(LIBDIR)/test_mpi_blr_ldlt_ooc.o
	$(MPIFC) -o $@ $(OPTF) $^ $(MPILIBS)

$(LIBDIR)/sym_testmat_mod.o: $(SRC)/sym_testmat_mod.f90 $(MAIN)/blr_kinds_mod.o
	$(MPIFC) $(OPTF) -c $(SRC)/sym_testmat_mod.f90 -o $@

$(LIBDIR)/test_mpi_blr_ldlt.o: $(TEST)/test_mpi_blr_ldlt.f90 $(LDLT_OBJS)
	$(MPIFC) $(OPTF) -c $(TEST)/test_mpi_blr_ldlt.f90 -o $@

$(BINDIR)/test_mpi_blr_ldlt: $(MAIN_OBJS) $(LDLT_OBJS) $(LIBDIR)/test_mpi_blr_ldlt.o
	$(MPIFC) -o $@ $(OPTF) $^ $(MPILIBS)

$(LIBDIR)/example_blr_mpi_ldlt.o: $(EX)/example_blr_mpi_ldlt.f90 $(LDLT_OBJS)
	$(MPIFC) $(OPTF) -c $(EX)/example_blr_mpi_ldlt.f90 -o $@

$(BINDIR)/example_blr_mpi_ldlt: $(MAIN_OBJS) $(LDLT_OBJS) $(LIBDIR)/example_blr_mpi_ldlt.o
	$(MPIFC) -o $@ $(OPTF) $^ $(MPILIBS)

$(LIBDIR)/bench_ldlt_vs_lu.o: $(TEST)/bench_ldlt_vs_lu.f90 $(LDLT_OBJS)
	$(MPIFC) $(OPTF) -c $(TEST)/bench_ldlt_vs_lu.f90 -o $@

$(BINDIR)/bench_ldlt_vs_lu: $(MAIN_OBJS) $(LDLT_OBJS) $(LIBDIR)/bench_ldlt_vs_lu.o
	$(MPIFC) -o $@ $(OPTF) $^ $(MPILIBS)

clean:
	rm -rf $(LIBDIR) $(BINDIR)

.PHONY: all clean
