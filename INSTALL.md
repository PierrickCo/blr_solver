# Installation

## Prérequis (chemin validé : Intel oneAPI)

Ce projet est développé et testé avec :

- **ifx** (Intel Fortran Compiler) — version testée : 2026.1.1
- **mpiifx** (wrapper MPI d'ifx, fourni par Intel oneAPI) — pour les modules
  MPI uniquement, aucune dépendance supplémentaire au runtime MPI tant que
  les modules purement séquentiels/OpenMP ne sont pas utilisés
- **Intel MKL** — version testée : 2026.1 (BLAS/LAPACK : `CGETRF`, `CGESVD`,
  `CGEMM`, `CTRSM`, etc.)
- **Intel MPI** — version testée : 2021.18

## Chemin alternatif : gcc + BLAS/LAPACK + OpenMPI

Pour qui n'a pas de licence Intel oneAPI, un second chemin de build existe
(`Makefile.gcc`), **validé fonctionnellement** (les 13 binaires de la ligne
vivante compilent et passent leurs tests, erreurs ~1e-7 comme sous Intel) :

```bash
make -f Makefile.gcc
# ou, si libopenblas-dev est installe (sudo apt install libopenblas-dev) :
make -f Makefile.gcc THREAD_IMPL=openblas
```

Testé avec `gfortran` 9.4.0, BLAS/LAPACK de référence Netlib
(`libblas-dev`/`liblapack-dev`, déjà présents sur la plupart des
distributions Debian/Ubuntu) et OpenMPI 4.0.3 (`libopenmpi-dev`,
`mpif90`). Tous les artefacts vont dans `build-gcc/` — jamais mélangés
avec les `.o`/`.mod` du build Intel (formats incompatibles entre
compilateurs).

Par défaut (`THREAD_IMPL=noop`), la bibliothèque BLAS/LAPACK de référence
est mono-thread — pas de garde-fou anti-sursouscription nécessaire (rien à
limiter). Avec `THREAD_IMPL=openblas` (si `libopenblas-dev` est installé),
le comportement redevient équivalent à MKL (garde-fou actif, threads
contrôlables).

**Piège rencontré et documenté, pas anodin** : sur cette machine de
développement (Windows + WSL2), compiler directement sur un chemin
`/mnt/c/...` (disque Windows monté via DrvFs) **échoue** avec `gfortran`
9.4.0 — le compilateur n'arrive pas à relire ses propres fichiers `.mod`
compressés (`f951: Fatal Error: ... Unexpected EOF`), alors que le fichier
est un flux gzip parfaitement valide (`gunzip -t` le confirme). C'est une
incompatibilité entre `gfortran`/DrvFs, PAS un bug du projet — confirmé en
recompilant à l'identique sur le système de fichiers natif WSL2 (`~/` ou
`/tmp/`), où tout fonctionne. **Un vrai cluster Linux n'a pas ce problème**
(système de fichiers natif) ; seul un développement local sous WSL2 avec
le dépôt sur le disque Windows peut le rencontrer — dans ce cas, travailler
depuis une copie sur le système de fichiers Linux natif (`~/...`) pour la
chaîne gcc spécifiquement (la chaîne Intel, elle, n'a pas ce problème).

## Chemin alternatif : ifort 2019 (Intel Parallel Studio XE, compilateur classique)

Pour un cluster n'ayant que l'ancienne generation de toolkit Intel
(`ifort`/`mpiifort`, pas `ifx`/`mpiifx` de oneAPI). **NON valide en
conditions reelles** (aucune installation Parallel Studio XE 2019
disponible dans cet environnement de developpement) — construit par
analogie stricte avec `Makefile.gcc` (deja valide) et `Makefile` (deja
valide), avec dry-run (`make -n`) verifie sur le graphe de dependances,
mais jamais reellement compile avec `ifort`. A tester avant usage en
production :

```bash
source /opt/intel/compilers_and_libraries_2019/linux/bin/compilervars.sh intel64
source /opt/intel/compilers_and_libraries_2019/linux/mkl/bin/mklvars.sh intel64
make -f Makefile.ifort2019
# ou, si MKLROOT/IMPI_LIBDIR different sur votre install :
make -f Makefile.ifort2019 MKLROOT=/chemin/vers/mkl IMPI_LIBDIR=/chemin/vers/mpi/intel64/lib
```

**Precautions specifiques a `ifort` 2019** (voir l'en-tete de
`Makefile.ifort2019` pour le detail) :
- Pas de `setvars.sh` (introduit avec oneAPI/2020) — sourcer
  `compilervars.sh` et `mklvars.sh` separement (ou `psxevars.sh` si
  disponible).
- Structure de dossiers MKL/Intel MPI differente de oneAPI : les
  bibliotheques sont sous `lib/intel64/`, pas `lib/` directement — deja
  pris en compte dans `Makefile.ifort2019`, mais a verifier si votre
  installation a une disposition non standard (surcharger `MKLROOT`/
  `IMPI_LIBDIR`).
- `mpi_blr_matrix_mod.f90` contient quelques lignes en forme libre
  depassant 132 colonnes (limite standard Fortran). `ifx` les tolere
  silencieusement ; le comportement d'`ifort` 2019 sur ce point precis
  n'est **pas verifie**. Si la compilation de ce fichier echoue avec une
  erreur de troncature/longueur de ligne, c'est la cause la plus probable
  — reformater les lignes incriminees avec une continuation `&` resout le
  probleme (aucun flag `ifort` connu ne desactive cette limite en forme
  libre, contrairement a `-ffree-line-length-none` sous `gfortran`).
- Tous les artefacts vont dans `build-ifort2019/`, jamais melanges avec
  ceux de `Makefile`/`Makefile.gcc` (format `.mod` non garanti compatible
  entre generations de compilateur).

## Compiler

```bash
source /opt/intel/oneapi/setvars.sh
cd blr_solver
make
```

Ceci construit les binaires de la ligne vivante (5 solveurs de référence +
BLR+MPI) dans `bin/`, à partir des `.o`/`.mod` compilés dans `lib/` (voir
le Makefile). Le code mis de côté (`legacy/`) se construit séparément,
**après** le projet principal (il réutilise directement les `.o`/`.mod`
déjà construits dans `lib/` plutôt que de les dupliquer) :

```bash
cd legacy
make
```

`make clean` nettoie les `.o`/`.mod`/binaires de chaque dossier
indépendamment.

Pour compiler les [exemples canoniques](examples/README.md) (un
programme court par paradigme, en plus de la ligne vivante) :

```bash
make examples             # ou : make -f Makefile.gcc examples
```

## Vérifier l'installation

Le plus simple : `./run_all_tests.sh` — compile puis lance les 13 binaires
de la ligne vivante (séquentiel, OpenMP, MPI, BLR+MPI) à petite échelle,
avec un résumé PASS/FAIL clair à la fin (pas besoin de connaître chaque
binaire/syntaxe individuellement). Logs détaillés dans
`run_all_tests_logs/` en cas d'échec.

Pour un test unitaire manuel, un cas minuscule, séquentiel, sans MPI :

```bash
./bin/test_block_lu 200 32
```

Doit se terminer par `=== TEST OK ===` (ou équivalent) avec une erreur
relative de l'ordre de `1e-7`.

Un cas MPI minimal (nécessite `mpirun` fonctionnel) :

```bash
mpirun -np 4 ./bin/test_mpi_blr_matrix 2 2 500 50 20 1e-4
```

Doit afficher un diagnostic d'affinité CPU (voir [`CLUSTER.md`](CLUSTER.md))
puis `=== TEST MPI_BLR_MATRIX : OK ===`.

## Référence des binaires et de leur syntaxe

Les binaires vivent dans `bin/` (`build-gcc/` pour la chaîne gcc, qui ne
sépare pas bin/lib) -- préfixer chaque commande ci-dessous par `bin/`.

| Binaire | Usage | Référence |
|---|---|---|
| `test_block_lu` | `bin/test_block_lu N NB` | A (unitaire) |
| `bench_blr_pipeline` | `bin/bench_blr_pipeline [N] [NRHS] [NB] [TOL]` (défauts 10000 100 100 1e-4 ; `TOL=0` pour B) | A / B |
| `bench_seq` | `bin/bench_seq [N] [NRHS] [NB] [TOL]` (`NB=N` pour la référence C standard) | C |
| `test_mpi_lu_pivot` | `mpirun -np <PR*PC> bin/test_mpi_lu_pivot PR PC N NB` | MPI (unitaire) |
| `bench_mpi` | `mpirun -np <PR*PC> bin/bench_mpi PR PC [N] [NRHS] [NB] [TOL]` | MPI |
| `test_mpi_blr_matrix` | `mpirun -np <PR*PC> bin/test_mpi_blr_matrix PR PC N NB NRHS TOL` | BLR+MPI, sans pivot |
| `test_mpi_blr_matrix_pivot` | `mpirun -np <PR*PC> bin/test_mpi_blr_matrix_pivot PR PC N NB NRHS TOL` | BLR+MPI, avec pivot confiné |
| `test_mpi_grid` | `mpirun -np <PR*PC> bin/test_mpi_grid PR PC N NB` | infrastructure de grille |
| `test_mpi_lu_solve_2d` | `mpirun -np <PR*PC> bin/test_mpi_lu_solve_2d PR PC N NB NRHS` | résolution 2D distribuée |
| `test_openmp_dense` | `bin/test_openmp_dense N NB` | mesure le threading MKL interne (référence C) |
| `test_blr_compress` | `bin/test_blr_compress` (sans argument) | bibliothèque d'algèbre par blocs |

`PR`/`PC` : dimensions de la grille de processus MPI (`PR*PC` doit égaler le
nombre de processus lancés par `mpirun -np`). `N` : taille de la matrice.
`NB` : taille de bloc (voir `CLUSTER.md` pour son réglage). `NRHS` : nombre
de seconds membres. `TOL` : tolérance de compression BLR (`0` = pas de
compression).
