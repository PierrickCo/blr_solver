#===============================================================================
# export_solution -- relit la solution de solveur_blr (--save-solution) et
# l'ecrit en texte, binaire et/ou NetCDF. Programme autonome et sequentiel :
# ni MPI, ni bibliotheque du solveur.
#
#   make                         ifx, avec NetCDF
#   make FC=gfortran             autre compilateur Fortran
#   make NETCDF=0                sans NetCDF (formats txt et bin seulement)
#   make NETCDF_LIBS="-L/chemin/lib -lnetcdf"   bibliotheque NetCDF ailleurs
#
# NetCDF : seule la bibliotheque C (libnetcdf) est necessaire, appelee
# directement depuis Fortran -- pas besoin de l'interface Fortran (netcdff).
#===============================================================================
FC ?= ifx
ifeq ($(FC),f77)
  FC = ifx
endif

NETCDF ?= 1
NETCDF_LIBS ?= -lnetcdf

ifneq (,$(findstring gfortran,$(FC)))
  FFLAGS ?= -O2
  FPP = -cpp
  MODF = -J build
else
  FFLAGS ?= -O2
  FPP = -fpp
  MODF = -module build
endif

OBJS = build/solution_reader_mod.o
ifeq ($(NETCDF),1)
  OBJS += build/netcdf_c_mod.o
  DEFS = -DWITH_NETCDF
  LIBS = $(NETCDF_LIBS)
endif

all: export_solution

export_solution: $(OBJS) build/export_solution.o
	$(FC) $(FFLAGS) -o $@ $^ $(LIBS)

build/%.o: src/%.f90 | build
	$(FC) $(FFLAGS) $(MODF) -c $< -o $@

build/export_solution.o: src/export_solution.F90 $(OBJS) | build
	$(FC) $(FFLAGS) $(FPP) $(DEFS) $(MODF) -c $< -o $@

build:
	mkdir -p build

clean:
	rm -rf build export_solution

.PHONY: all clean
