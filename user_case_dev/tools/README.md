# tools/ -- profilage de `example_blr_mpi_pivot_user`

Deux scripts Python qui lisent les CSV produits par `profiling/profiling_mod.f90`
et en tirent des graphiques. Ils ne font aucun calcul BLR eux-memes -- ce
sont de purs outils de visualisation, a lancer APRES un run du solveur
avec les options `--profile-usage` et/ou `--dump-structure`.

## Prerequis

```bash
pip install matplotlib numpy
```

(`numpy` n'est necessaire que pour `plot_blr_structure.py`.)

## 1) Suivi memoire/CPU par etape -- `--profile-usage`

Chaque processus MPI ecrit son propre `usage_profile_rank<N>.csv`
(colonnes `rank,stage,t_wall_s,rss_mb,cpu_pct`) : un echantillon a
chaque etape du solveur (debut, assemblage, factorisation, second
membre, resolution, verification, fin), PLUS un echantillon
supplementaire tous les 5% de progression a l'interieur de
l'assemblage/la factorisation/la resolution (etiquettes `etape@NN%`).

```bash
# depuis user_case_dev/, apres `make` (voir ../Makefile)
mpirun -np 4 ./example_blr_mpi_pivot_user 2 2 4000 200 4 1e-4 --profile-usage
python3 tools/plot_usage_profile.py . usage_profile.png
```

Produit `usage_profile.png` : deux graphiques empiles contre le temps
ecoule -- RSS (Mo) par processus, et %CPU par processus (peut depasser
100% si plusieurs threads MKL/OpenMP sont actifs, meme convention que
`top`) -- avec les bornes d'etape (assemblage/factorisation/...) reperees
par des traits verticaux.

## 2) Structure de compression -- `--dump-structure`

Chaque processus ecrit la nature de ses blocs LOCAUX (diagonal / dense /
bas-rang, et le rang associe) a la fin de l'assemblage ET a la fin de la
factorisation, dans `blr_structure_assembly_rank<N>.csv` et
`blr_structure_factorization_rank<N>.csv`.

```bash
mpirun -np 4 ./example_blr_mpi_pivot_user 2 2 4000 200 4 1e-4 --dump-structure
python3 tools/plot_blr_structure.py . blr_structure.png
```

Produit `blr_structure.png` : une figure a deux panneaux (fin
assemblage / fin factorisation), un pixel par bloc GLOBAL de VOTRE
matrice (celle definie dans `matrix/user_matrix_mod.f90`), colore par
rang (blocs bas-rang) ou par categorie fixe (diagonal en gris fonce,
dense en rouge).

## Les deux en meme temps

```bash
mpirun -np 4 ./example_blr_mpi_pivot_user 2 2 4000 200 4 1e-4 --profile-usage --dump-structure
python3 tools/plot_usage_profile.py .
python3 tools/plot_blr_structure.py .
```

## Notes

- Les fichiers CSV/PNG sont ecrits dans le repertoire courant au moment
  du run/du plot -- passez un chemin en premier argument aux deux
  scripts pour pointer ailleurs.
- Jamais relus par le solveur -- a supprimer librement entre deux runs
  (`rm usage_profile_rank*.csv blr_structure_*.csv`).
