# Exemples canoniques

Sept programmes courts, un par paradigme du projet (voir le tableau de
[`README.md`](../README.md)), montrant l'enchainement minimal d'appels
pour chaque solveur -- sans l'instrumentation de mesure ni la validation
croisee contre une reference independante des drivers de `test/`/`bench/`
(qui restent la source de verite pour la correction, voir
`../run_all_tests.sh`).

| Fichier | Reference | Paradigme |
|---|---|---|
| `example_blr_sequentiel.f90` | A/B | Sequentiel/OpenMP, matrice compressee (BLR) de bout en bout |
| `example_dense_sequentiel.f90` | C | Sequentiel/OpenMP, LU dense par blocs, pivot complet |
| `example_mpi_dense.f90` | MPI | MPI, LU dense distribue, pivot, resolution 2D |
| `example_blr_mpi_sans_pivot.f90` | BLR+MPI | MPI, stockage BLR distribue, sans pivot (diagonale dominante) |
| `example_blr_mpi_pivot.f90` | BLR+MPI | MPI, stockage BLR distribue, pivot confine au bloc diagonal |
| `example_blr_mpi_pivot_arena.f90` | BLR+MPI | Identique, stockage bas-rang carve dans un allocateur arena (voir `../README.md`) |
| `example_blr_mpi_pivot_ooc.f90` | BLR+MPI | Identique, out-of-core -- matrice/facteurs/second membre sur disque, RAM bornee par `MAX_MB_PER_PROCESS` (voir plus bas) |

## Compiler

```bash
cd ..                    # depuis examples/, remonter a blr_solver/
make examples             # chaine Intel (ifx+MKL+Intel MPI)
# ou
make -f Makefile.gcc examples   # chaine gcc+BLAS/LAPACK+OpenMPI
```

Chaque exemple se construit comme binaire dans `blr_solver/bin/`
(ou `build-gcc/` pour la chaine gcc), en reutilisant directement les
`.o`/`.mod` deja compiles par `make`/`make -f Makefile.gcc` dans
`blr_solver/lib/` -- compiler le projet principal d'abord.

## Lancer

```bash
./bin/example_blr_sequentiel [N] [NB] [NRHS] [TOL]
./bin/example_dense_sequentiel [N] [NB]
mpirun -np 4 ./bin/example_mpi_dense 2 2 [N] [NB] [NRHS]
mpirun -np 4 ./bin/example_blr_mpi_sans_pivot 2 2 [N] [NB] [NRHS] [TOL]
mpirun -np 4 ./bin/example_blr_mpi_pivot 2 2 [N] [NB] [NRHS] [TOL] [OPTIONS...]
mpirun -np 4 ./bin/example_blr_mpi_pivot_ooc 2 2 [N] [NB] [NRHS] [TOL] [MAX_MB_PER_PROCESS] [NRHS_BATCH] [SCRATCH_DIR] \
    [EPS] [RATIO] [USE_RSVD] [CHECK_ERROR] [CHECK_RESIDUAL] [PROFILE_USAGE] [DUMP_STRUCTURE] \
    [SAVE_MATRIX] [SAVE_RHS] [SAVE_FACTORED] [SAVE_SOLUTION] [LOAD_MATRIX] [LOAD_FACTORED] [LOAD_RHS] [CHECKPOINT_DIR]
```

`example_blr_mpi_pivot` seul : au-dela de `PR PC [N] [NB] [NRHS] [TOL]`
(positionnels, dans cet ordre), toutes les autres options sont NOMMEES
(`--cle=valeur`, ou `--cle` seul pour une booleenne = `--cle=1`), dans
n'importe quel ordre -- `--help`/`-h` affiche la liste complete. Voir
`## Sauvegarde et reprise (checkpoint/restart)` plus bas pour le detail
des options de checkpoint (communes aux deux exemples BLR+MPI avec pivot
-- `example_blr_mpi_pivot_ooc` garde, lui, la forme positionnelle
d'origine).

Chaque exemple assemble une matrice de demonstration (formule analytique,
pas de fichier d'entree), factorise, resout un second membre arbitraire,
et affiche les premieres valeurs de la solution -- point de depart a
adapter (callback `GET_BLOCK`, second membre) pour un probleme reel.

## Empreinte memoire theorique (solveurs MPI)

Les 3 exemples MPI (`example_mpi_dense.f90`, `example_blr_mpi_sans_pivot.f90`,
`example_blr_mpi_pivot.f90`) affichent, a 3 moments cles, l'empreinte
memoire THEORIQUE (calculee a partir des tailles de tableaux -- PAS
mesuree via l'OS, voir `GET_RSS_KB`/`CLUSTER.md` pour la RSS reelle) :

1. **Matrice** (dense ou compressee BLR) juste apres l'assemblage.
2. **Facteurs LU** juste apres la factorisation -- pour la reference MPI
   (dense), c'est EXACTEMENT la meme valeur qu'a l'assemblage (facteurs
   ecrasent la matrice en place, aucune compression ne peut changer la
   taille) ; pour BLR+MPI, cette valeur peut DIFFERER de celle de
   l'assemblage a cause du remplissage ("fill-in") introduit par la mise
   a jour de Schur -- observe concretement sur `example_blr_mpi_pivot`
   (matrice 1.9 Mo -> facteurs 2.0 Mo a N=1000/NB=100/grille 2x2).
3. **Second membre B** juste apres son allocation/remplissage.

Chaque ligne affiche DEUX nombres (`REPORT_MEMORY_MB`, nouvelle routine
partagee de `mpi_grid_mod.f90`) : le **total** agrege sur toute la grille
(`MPI_ALLREDUCE`/`MPI_SUM`, utile pour la capacite globale du cluster) et
le **maximum par-processus** (`MPI_ALLREDUCE`/`MPI_MAX`, c'est LUI qui
determine si le job tient dans la RAM de chaque noeud -- la contrainte
qui bloque reellement un job, pas le total). Pour BLR+MPI, la memoire de
la matrice/des facteurs reutilise directement `MPI_BLR_MATRIX_MEMORY_
BYTES` (deja existante) ; pour la reference MPI dense et le second
membre (les deux BLR+MPI inclus), c'est une simple formule
`lignes_locales * colonnes_locales * 8 octets` (`COMPLEX(WP)` simple
precision).

**Comparaison a la RSS reelle** (`example_blr_mpi_pivot.f90` uniquement,
via `GET_RSS_KB()`/`bench_utils_mod.f90`, meme mecanisme `REPORT_MEMORY_
MB`) : la RSS mesuree est TOUJOURS superieure a la valeur theorique
(elle inclut le runtime MPI/MKL charge en memoire, les workspaces
SVD/QR transitoires par thread, l'executable lui-meme...), mais le
RATIO n'est PAS constant -- mesure sur ce projet : ~65x a N=2000 (l'
overhead fixe du runtime domine a cette echelle), ~2.5-3x a N=20000
(les donnees elles-memes dominent). C'est ce second regime (grande
echelle) qui sert de reference pour detecter une vraie anomalie (fuite,
copie transitoire oubliee) -- pas un ratio absolu independant de N.

## Barres de progression et taux de compression

Les 5 exemples activent `SHOW_PROGRESS=.TRUE.` sur chacune de leurs etapes
couteuses (assemblage, factorisation, resolution) -- couverture complete
des 5 solveurs de reference, pas seulement BLR+MPI. Seule exception :
`BLOCK_LU_SOLVE` (reference C, dense sequentiel) n'a pas de barre --
elle resout en un seul appel `CTRSM` par triangle, sans boucle par bloc a
instrumenter (deja tres rapide, `O(N^2)`, pas `O(N^3)` comme la
factorisation).

Les deux exemples BLR+MPI (`example_blr_mpi_sans_pivot.f90`,
`example_blr_mpi_pivot.f90`) affichent en plus le taux de compression
(`N_COMPRESSED/N_OFFDIAG`) apres l'assemblage et a la fin de la
factorisation. Ces deux compteurs sont deja des sommes GLOBALES (agregees
par `MPI_ALLREDUCE` a l'interieur de `mpi_blr_matrix_mod.f90`) --
identiques sur tous les processus, affiches une seule fois (rang 0).

**Deux mecanismes d'affichage distincts** (voir `src/progress_bar_mod.f90`) :
- **Sequentiel/OpenMP** (`example_blr_sequentiel.f90`,
  `example_dense_sequentiel.f90`) : `PRINT_PROGRESS_BAR`, une barre a
  rafraichissement en place (retour chariot `CHAR(13)`) -- fonctionne
  correctement, attache a un vrai terminal.
- **MPI** (les 3 autres exemples, et les modules `mpi_lu_pivot_mod.f90`/
  `mpi_lu_solve_2d_mod.f90`/`mpi_blr_matrix_mod.f90` sous-jacents) :
  `PRINT_PROGRESS_MILESTONE`, un affichage par PALIERS (une ligne complete
  tous les 5%) -- `mpirun` relaie le stdout de chaque processus via un
  pipe, pas un pseudo-terminal, donc le retour chariot n'y est PAS
  interprete comme "ecraser la ligne precedente" : la barre en place
  finissait par polluer la sortie d'une ligne par mise a jour (bug
  constate et corrige) -- d'ou l'affichage par paliers, robuste quel que
  soit le lanceur.

Pour les 3 boucles de factorisation/resolution MPI (reference MPI et
BLR+MPI, avec ou sans pivot), le pourcentage affiche est EXACT : chaque
etape K est globalement synchronisee (un `MPI_BCAST`/diffusion a chaque
K), donc `K/NBLK` (ou l'equivalent pour la resolution, 2 sous-etapes)
represente le meme etat sur tous les processus. Seul l'assemblage BLR+MPI
fait exception : aucune synchronisation entre processus a cette etape
(chacun assemble ses blocs locaux de son cote) -- l'affichage du rang 0
n'est donc qu'une approximation (mais fidele, le decoupage bloc-cyclique
equilibrant la charge entre processus).

## Investigation de l'ecart RSS-vs-theorique (`example_blr_mpi_pivot.f90`)

A N=20000/NB=500/grille 2x2, la factorisation ajoute ~49 Mo/processus de
RSS reelle alors que le stockage BLR (theorique) ne grossit que de ~10 Mo
(remplissage/fill-in) -- un ecart de ~39 Mo non explique par les donnees
elles-memes. Quatre hypotheses ciblees ont ete testees, TOUTES ecartees
ou negligeables :

1. **Cache interne MKL** (`mkl_free_buffers()`, appele apres la
   factorisation) -- ~15% de l'ecart recupere, pas la cause principale.
2. **Fragmentation en haut de tas** (`malloc_trim(0)`, glibc) -- ~0%
   recupere (il ne peut agir que sur la memoire libre CONTIGUE en haut du
   tas, pas les trous au milieu).
3. **Reallocation implicite Fortran** (`ABLK = NEWBLK` declenche une
   reallocation+copie supplementaire meme si `NEWBLK` est deja a la bonne
   taille) -- corrige via `MOVE_ALLOC` (voir `SCHUR_UPDATE_BLR_BLOCK`,
   `blr_compress_mod.f90`, argument optionnel `USE_MOVE_ALLOC`,
   propage par `MPI_BLR_MATRIX_LU_FACTOR_PIVOT`) -- resultat numerique
   IDENTIQUE (verifie bit-a-bit), mais ~0% de gain memoire mesure.
4. **Seuil `mmap` glibc** (`mallopt(M_MMAP_THRESHOLD, ...)`, testé puis
   retire -- voir l'historique git) -- l'appel reussit (code retour 1)
   mais ~0% d'effet, tres probablement parce que le runtime Fortran
   d'Intel (`ifx`) ne route pas les `ALLOCATE`/`DEALLOCATE` par le
   `malloc()` de glibc que ce reglage controle (non confirme -- non
   teste sous gfortran, qui route typiquement `ALLOCATE` directement par
   `malloc()`).

**Conclusion** : un tracage RSS par etape K (`SHOW_RSS_TRACE`, argument
optionnel de `MPI_BLR_MATRIX_LU_FACTOR_PIVOT`, actif dans cet exemple)
montre une croissance CONTINUE qui RALENTIT et se stabilise vers la fin
de la factorisation -- une signature de fragmentation proportionnelle au
VOLUME de travail (le nombre de blocs mis a jour a chaque etape diminue
mecaniquement quand la sous-matrice traINante retrecit), pas d'une fuite
qui grossirait sans lien avec le travail restant. Aucun allocateur
standard (glibc, jemalloc, tcmalloc) ne peut compacter (deplacer) de la
memoire encore utilisee -- seule une reecriture manuelle (arena/pool
dedie sur toute la duree de la factorisation, pas juste une etape --
le stockage des blocs doit survivre au-dela d'une seule etape K)
offrirait une garantie reelle, avec son propre compromis (memoire
gaspillee par les anciennes versions de blocs jamais recuperees) --
chantier volontairement laisse de cote pour l'instant, le rapport
effort/gain n'etant pas demontre.

`USE_MOVE_ALLOC` est conserve (defaut `.FALSE.`, aucun changement de
comportement) car c'est une bonne pratique defensive independante du
gain memoire observe ici. `SHOW_RSS_TRACE` est conserve comme outil de
diagnostic permanent (defaut `.FALSE.`, cout nul).

## Solveur out-of-core (`example_blr_mpi_pivot_ooc.f90`)

Meme algorithme numerique EXACT que `example_blr_mpi_pivot.f90` (pivot
partiel confine au bloc diagonal, meme matrice de demonstration) -- mais
plus RIEN n'est suppose tenir en RAM en permanence : chaque processus MPI
possede un fichier scratch local (jamais partage) et RIEN d'autre que de
petits descripteurs (`M`/`N`/`RANG`/position sur fichier) ne reste
resident entre deux etapes. Voir `MPI_BLR_MATRIX_OOC_T` dans
`src/mpi_blr_matrix_mod.f90` et `src/blr_ooc_mod.f90` pour le detail
d'implementation (format de stockage, cache LRU, politique d'ecriture).

**4 etapes SEPAREES, chacune bornee en memoire par un parametre
utilisateur :**

1. **Assemblage de A** -> ecrit sur le fichier scratch, borne par
   `MAX_MB_PER_PROCESS` (cache LRU).
2. **Assemblage du second membre B** -> ECRIT SUR DISQUE LUI AUSSI, par
   paquets d'au plus `NRHS_BATCH` colonnes -- jamais tout B en RAM a la
   fois, meme dans un cas massivement multi-second-membre.
3. **Factorisation LU** -> les blocs "finis" (panneau/diagonale de
   l'etape K) repartent sur disque des qu'ils ne servent plus a la mise a
   jour de Schur ; la sous-matrice traineante reste bornee par
   `MAX_MB_PER_PROCESS` via le meme cache.
4. **Resolution** -> relit L/U depuis le disque, par paquets d'au plus
   `NRHS_BATCH` seconds membres.

### Choisir `MAX_MB_PER_PROCESS` -- NE PAS prendre le minimum

**Ne pas** choisir ce parametre "le plus petit possible" pour economiser
de la RAM -- un budget sous l'ensemble de travail REEL du processus ne
degrade pas gracieusement, il fait THRASHER le cache LRU, avec un cout
disproportionne par rapport a la RAM economisee.

**Cause** : un bloc hors-diagonale n'est retouche par une mise a jour de
Schur qu'une fois tous les ~`NBLK` pas de la factorisation, et la
TOTALITE de l'ensemble de travail local est parcourue entre deux touches
du MEME bloc -- un schema d'acces quasi-SEQUENTIEL ou aucune taille de
cache inferieure a l'ensemble de travail complet n'aide, quelle que soit
la politique d'eviction.

**Mesure concrete** (N=20000, NB=400, grille 2x2, voir
`test/bench_blr_mpi_pivot_ooc_compare.f90`) :

| `MAX_MB_PER_PROCESS` | RSS reelle | Temps cache (lookup+E/S disque) | Temps factorisation totale |
|---|---|---|---|
| 200 Mo | ~127 Mo | 10.1 s (68% du temps) | 14.9 s |
| 800 Mo | ~128 Mo | 0.75 s (÷13.5) | 3.1 s (÷4.9) |

La RSS reelle est restee QUASI IDENTIQUE entre les deux budgets -- la
difference de temps ne vient PAS d'un budget "un peu court avec un peu
plus de RAM utilisee en echange", elle vient du fait que 200 Mo etait
EN-DESSOUS de l'ensemble de travail reel (~130 Mo) et 800 Mo est
au-dessus : sous ce seuil, le cache thrashe ; au-dessus, le cout cache
s'effondre et le temps de factorisation converge vers le chemin
classique (calcul pur deja comparable : 0.97 s OOC contre 0.61 s
classique a budget=200 deja).

**Conclusion pratique** : choisir `MAX_MB_PER_PROCESS` au moins egal a
l'ensemble de travail reel du processus (verifiable en augmentant le
budget jusqu'a ce que la RSS reelle cesse de croitre, ou via
`MPI_BLR_MATRIX_OOC_RESIDENT_BYTES` a un budget genereux). L'interet de
l'out-of-core n'est PAS de tourner avec le budget minimal possible, mais
de rester CORRECT et de degrader la RAM au lieu de planter quand le
probleme est trop gros pour tenir en RAM du tout -- en-dessous du seuil
"ensemble de travail", le compromis temps/memoire devient tres
defavorable et s'aggrave a mesure que N grandit.

### Diagnostic fin (`T_CACHE_GET`/`T_CACHE_EVICT`/`T_COMPUTE`/`T_MPICALL`)

`MPI_BLR_MATRIX_LU_FACTOR_PIVOT_OOC` accepte 4 arguments de sortie
optionnels (aucun surcout si absents) qui decomposent `T_FACTOR` en 4
categories qui se somment exactement :
- `T_CACHE_GET` : temps dans TOUS les acces cache (lookup + lecture
  disque si le bloc n'est pas deja resident).
- `T_CACHE_EVICT` : temps des evictions de panneau (etape 3,
  TOUJOURS reelles, ecriture disque inconditionnelle quel que soit le
  budget).
- `T_COMPUTE` : calcul pur (`CGETRF`, TRSM, mises a jour de Schur, swap
  de ligne) -- si ce terme seul est deja comparable au chemin classique
  mais que `T_FACTOR` reste tres superieur, c'est le signe que
  `T_CACHE_GET` domine (voir ci-dessus).
- `T_MPICALL` : diffusions MPI (`BCAST` de `IPIV_K`/`DKK`, diffusion des
  panneaux compresses).

### Cache de resolution (`MPI_BLR_MATRIX_LU_SOLVE_PIVOT_OOC_2D`)

La resolution relit `DIAGBLK`/`OFFBLK` depuis le disque via `CACHE_A`/
`CACHE_D` -- les MEMES caches que la factorisation, reutilises tels quels
(vides mais structurellement valides apres le flush de fin de
factorisation). A l'interieur d'un SEUL paquet, un bloc L/U n'est touche
qu'une fois (pas de gain a esperer la) -- mais les MEMES blocs sont
necessaires a CHAQUE paquet (`NBATCHES` fois au total, ils ne dependent
pas des colonnes de second membre), donc un cache qui SURVIT au paquet en
cours evite de les relire depuis le disque a chaque fois.

**Mesure** (N=20000, NB=400, matrice quasi incompressible -- 1% compresse,
donc blocs denses volumineux --, `MAX_MB_PER_PROCESS`=900,
`NRHS_BATCH`=2 donc `NBATCHES`=3) : temps de resolution divise par **2.76**
(36.5 s -> 13.2 s), RSS quasi inchangee. Gain plus modeste que le facteur
`NBATCHES`=3 aurait pu suggerer -- coherent avec la meme regle que pour la
factorisation : `MAX_MB_PER_PROCESS`=900 est proche mais pas tres
au-dessus du besoin reel observe (~915-920 Mo), une partie du gain se
perd donc encore en eviction. Sur une matrice tres compressible (blocs
petits), l'effet est negligeable dans les deux sens -- ce cache n'apporte
un vrai gain que lorsque NBATCHES>1 ET que les blocs sont assez gros/
nombreux pour que relire NBATCHES fois coute reellement cher.

## Sauvegarde et reprise (checkpoint/restart)

`example_blr_mpi_pivot.f90` et `example_blr_mpi_pivot_ooc.f90` acceptent
7 options communes pour sauvegarder un objet a un instant donne et
reprendre un calcul plus tard sans repasser par les etapes deja faites --
nommees (`--save-matrix=1`, n'importe quel ordre) pour `example_blr_mpi_
pivot.f90`, positionnelles (`SAVE_MATRIX=1`, ordre fixe) pour la variante
out-of-core (voir `## Lancer` plus haut pour le detail des deux styles) :

| Option (nommee / positionnelle OOC) | Effet |
|---|---|
| `--save-matrix` / `SAVE_MATRIX=1` | Sauvegarde la matrice assemblee (avant factorisation) |
| `--save-rhs` / `SAVE_RHS=1` | Sauvegarde le second membre B (avant resolution) |
| `--save-factored` / `SAVE_FACTORED=1` | Sauvegarde la matrice factorisee (facteurs LU + `IPIV`) |
| `--save-solution` / `SAVE_SOLUTION=1` | Sauvegarde la solution X (apres resolution) |
| `--load-matrix` / `LOAD_MATRIX=1` | Relit la matrice sauvegardee, saute l'assemblage |
| `--load-factored` / `LOAD_FACTORED=1` | Relit la matrice factorisee, saute assemblage ET factorisation (prioritaire sur `--load-matrix`/`LOAD_MATRIX` si les deux sont positionnes) |
| `--load-rhs` / `LOAD_RHS=1` | Relit le second membre sauvegarde, saute son calcul |
| `--checkpoint-dir=.` / `CHECKPOINT_DIR=.` | Repertoire des fichiers de checkpoint (defaut `.`, doit deja exister) |

**Format** : UN FICHIER PAR PROCESSUS MPI (`STREAM` non formate,
`saved_<objet>_rank<N>.dat` pour `example_blr_mpi_pivot.f90`, ou
`saved_<objet>_ooc_rank<N>.dat` pour la variante out-of-core -- prefixes
DISTINCTS pour eviter tout chargement croise accidentel entre les deux
exemples). Chaque fichier commence par un en-tete (N/NB ou NRHS, grille
PR/PC/MYROW/MYCOL/MYID) verifie AVANT tout chargement : en cas
d'incoherence (grille differente, N/NB differents...) sur NE SERAIT-CE
QU'UN SEUL rang, TOUS les rangs s'arretent proprement (reduction
`MPI_ALLREDUCE`/`MPI_LAND` collective) -- jamais de blocage partiel ni de
calcul silencieusement corrompu.

**Variante out-of-core** : `SAVE_BLR_MATRIX_OOC`/`LOAD_BLR_MATRIX_OOC`
vidangent d'abord les caches (`OOC_CACHE_FLUSH_ALL`/`OOC_DIAG_CACHE_
FLUSH_ALL`) puis relisent/reecrivent chaque bloc UN A LA FOIS via
`OOC_READ_BLOCK_RAW`/`OOC_WRITE_BLOCK_FIXED` -- jamais plus d'un bloc
resident en RAM, meme pendant la sauvegarde/le chargement (l'avantage
memoire de l'OOC reste entier). `SAVE_SOLUTION` n'a pas d'equivalent RAM
a relire pour cette variante (X n'existe jamais en un seul tenant, voir
`CONSUME_X_BATCH`) : chaque paquet resolu est ecrit sur disque des qu'il
est disponible, au fil de l'eau, directement depuis le callback. Ces
fichiers de checkpoint sont DISTINCTS des fichiers scratch internes
(`a_rank_*.dat`/`b_rank_*.dat` dans `SCRATCH_DIR`) : ces derniers restent
du scratch pur, supprimes automatiquement en fin d'execution -- les
checkpoints, eux, sont PERSISTANTS.

**Exemple** (`example_blr_mpi_pivot.f90` -- assembler+factoriser une
fois, resoudre plus tard sans refactoriser) :

```bash
mkdir -p ckpt

# 1) assemble, factorise, sauvegarde la matrice factorisee
mpirun -np 4 ./bin/example_blr_mpi_pivot 2 2 4000 200 4 1e-4 --save-factored --checkpoint-dir=ckpt

# 2) plus tard : relit directement la matrice factorisee, saute
#    assemblage ET factorisation
mpirun -np 4 ./bin/example_blr_mpi_pivot 2 2 4000 200 4 1e-4 --load-factored --checkpoint-dir=ckpt
```

Meme principe pour la variante out-of-core (deux arguments
supplementaires avant `EPS`, `MAX_MB_PER_PROCESS` et `NRHS_BATCH`, plus
`SCRATCH_DIR` -- voir la ligne `Lancer` plus haut pour l'ordre complet) :

```bash
mkdir -p scratch ckpt
mpirun -np 4 ./bin/example_blr_mpi_pivot_ooc 2 2 4000 200 4 1e-4 64 2 scratch 0.0 0.333333 0 0 0 0 0 0 0 1 0 0 0 0 ckpt
mpirun -np 4 ./bin/example_blr_mpi_pivot_ooc 2 2 4000 200 4 1e-4 64 2 scratch 0.0 0.333333 0 0 0 0 0 0 0 0 0 0 1 0 ckpt
```
