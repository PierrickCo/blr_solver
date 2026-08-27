# Code mis de côté (`legacy/`)

Ce dossier contient les étapes intermédiaires du projet qui ne font plus
partie de la ligne de référence actuelle (voir le `README.md` racine),
mais qui restent compilables et fonctionnelles — gardées pour leur valeur
pédagogique (elles montrent le cheminement qui a mené à la conception
actuelle) et comme base de départ si l'une de ces pistes est reprise un
jour.

**Aucun de ces modules n'est utilisé par le code vivant** (`src/`,
`test/`), à une seule exception près : `mpi_blr_lu_pivot_bcast_mod.f90`
(proposition 3) est resté dans `src/` — pas ici — car sa routine
`BCAST_BLR_BLOCK` est réutilisée par `mpi_blr_matrix_mod.f90`. Seule sa
propre routine de factorisation historique (`MPI_BLR_LU_PIVOT_BCAST_FACTOR`)
et son test dédié sont ici.

## Compiler

Compiler d'abord le projet principal (`make` depuis `blr_solver/`), **puis**
ce dossier (`make` depuis `blr_solver/legacy/`) — ce Makefile réutilise
directement les `.o`/`.mod` déjà construits par le projet principal pour
les modules encore vivants (`blr_kinds_mod`, `mpi_grid_mod`,
`blr_compress_mod`, `block_lu_mod`, `mpi_lu_pivot_mod`,
`mpi_blr_lu_pivot_bcast_mod`, `mpi_lu_solve_2d_mod`, `bench_utils_mod`),
plutôt que de dupliquer leur code source.

```bash
cd blr_solver && make
cd legacy && make
```

## Contenu, par étape historique

- **LDLt (abandonné)** — `block_ldlt_mod.f90` (LDLt dense, Bunch-Kaufman
  via `CLASYF`, fonctionnel) et `blr_ldlt_mod.f90` (tentative
  LDLt+compression BLR, **non fonctionnelle** — `CLASYF` utilise un schéma
  de mise à jour différée qui ne se décompose pas comme pour LU ; aurait
  demandé de réécrire Bunch-Kaufman à la main pour exposer les blocs à
  compresser. Voir l'en-tête du fichier pour le détail de l'échec).
- **Étape 2/5 + proposition 2** — `blr_lu_mod.f90` (BLR séquentielle,
  facteurs stockés denses, compression n'accélère que la factorisation) et
  `blr_lu_storage_mod.f90` (facteurs stockés compressés, en deux phases).
  Remplacées par le pipeline « jamais dense » de `blr_matrix_mod.f90`
  (compression dès l'assemblage, pas de phase de recompression séparée).
- **6a-6e (première lignée MPI, stockage local toujours dense)** —
  `mpi_matvec_mod.f90` (6b, produit matrice-vecteur distribué),
  `mpi_lu_nopiv_mod.f90` (6c, LU distribué sans pivot),
  `mpi_blr_lu_pivot_mod.f90` (6e, BLR+pivot, diffusion dense-puis-
  compression), `mpi_blr_lu_nopiv_mod.f90` (BLR+MPI sans pivot, complète
  la grille des 6 combinaisons dense/BLR × OpenMP/MPI historiquement
  demandée). **6a** (`mpi_grid_mod.f90`) et **6d**
  (`mpi_lu_pivot_mod.f90`) restent vivants — ce sont la base de la
  référence MPI actuelle, pas du code mis de côté.
- **Résolution 1D** — `mpi_lu_solve_mod.f90`, remplacée par la résolution
  2D complète (`mpi_lu_solve_2d_mod.f90`, vivante) qui distribue le second
  membre par lignes ET par colonnes sans réplication résiduelle.
- **Benchmarks complets** — `bench_seq_full.f90`/`bench_mpi_full.f90`
  (copies non tronquées de `bench_seq.f90`/`bench_mpi.f90` avant que ceux-ci
  ne soient recentrés sur les seules références C et MPI) — comparent en
  plus les variantes ci-dessus et l'ancienne résolution par rassemblement.

Chaque module garde ses commentaires d'en-tête d'origine, expliquant le
raisonnement et les décisions prises au moment de sa construction.
