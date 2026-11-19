# export_solution

Relit la solution X écrite par `solveur_blr ... --save-solution` (un fichier
`saved_solution_rank*.dat` par processus) et l'écrit dans **un seul fichier**
par format choisi : texte, binaire et/ou NetCDF.

Programme autonome et séquentiel : il ne dépend ni de MPI ni de la
bibliothèque du solveur, et se lance sans `mpirun`, après le calcul.

## Compilation

```bash
source /opt/intel/oneapi/setvars.sh     # pour ifx
make                                    # avec NetCDF
make NETCDF=0                           # sans NetCDF (txt et bin seulement)
make FC=gfortran                        # autre compilateur
make NETCDF_LIBS="-L/chemin/lib -lnetcdf"   # NetCDF installé ailleurs
```

Seule la bibliothèque **C** de NetCDF (`libnetcdf`) est nécessaire : elle est
appelée directement depuis Fortran (`src/netcdf_c_mod.f90`), l'interface
Fortran de NetCDF (`netcdff`) n'est pas utilisée. Si seule
`libnetcdf.so.<version>` existe (paquet d'exécution sans `-dev`), donner son
chemin complet : `make NETCDF_LIBS=/usr/lib/x86_64-linux-gnu/libnetcdf.so.15`.

## Utilisation

```bash
# 1) calcul, en sauvegardant la solution
mpirun -np 4 ../solveur_pcmatgen/solveur_blr 2 2 4000 200 8 1e-4 \
       --save-solution --checkpoint-dir=resultats
# 2) export
./export_solution --dir=resultats --nb=200 --format=txt,netcdf --out=X
```

| Option | | Rôle |
|---|---|---|
| `--nb=NB` | **obligatoire** | taille de bloc NB du calcul (voir l'avertissement ci-dessous) |
| `--format=LISTE` | **obligatoire** | `txt`, `bin`, `netcdf` séparés par des virgules, ou `all` |
| `--dir=CHEMIN` | défaut `.` | dossier des `saved_solution_rank*.dat` (le `--checkpoint-dir` du calcul) |
| `--out=PREFIXE` | défaut `solution` | nom des fichiers produits, sans extension |
| `--rows=I1:I2` | défaut `1:N` | lignes exportées |
| `--cols=J1:J2` | défaut `1:NRHS` | seconds membres exportés |
| `--chunk=K` | défaut 16 | colonnes lues à la fois : la mémoire est d'environ lignes × K complexes |

N, NRHS et la grille PR×PC sont lus dans les fichiers.

**Attention à NB** : il n'est pas stocké dans les fichiers du solveur. Un NB
faux n'est détecté que s'il change le nombre de lignes par processus ; sinon
(par exemple N = 1000, PR = 2, NB = 50 au lieu de 100) X est exporté dans le
désordre **sans message d'erreur**.

## Formats produits

- **`PREFIXE.txt`** : quelques lignes d'en-tête commençant par `#`, puis une
  ligne par élément, colonne par colonne : `i  j  Re(X(i,j))  Im(X(i,j))`
  (indices globaux à partir de 1, 9 chiffres significatifs).
- **`PREFIXE.bin`** (accès direct, sans marqueurs Fortran) : 4 entiers 32 bits
  `NROWS NCOLS IMIN JMIN`, puis X en complexes simple précision (réel puis
  imaginaire, 2 × 32 bits), colonne par colonne.
  Lecture en Python : `h = np.fromfile(f, np.int32, 4)` puis
  `X = np.fromfile(f, np.complex64, offset=16).reshape(h[1], h[0]).T`.
- **`PREFIXE.nc`** (NetCDF-4) : dimensions `row` et `col` ; variables
  `x_re(col,row)` et `x_im(col,row)` (float), `row` et `col` (indices globaux) ;
  attributs globaux `N`, `NRHS`, `title`. En C/Python l'ordre est
  `x_re[col][row]` ; en Fortran, `x_re(row,col)`.

## Contenu

```
Makefile
src/export_solution.F90       programme et options
src/solution_reader_mod.f90   lecture des fichiers par processus (format et répartition)
src/netcdf_c_mod.f90          appels à la bibliothèque C de NetCDF
```
