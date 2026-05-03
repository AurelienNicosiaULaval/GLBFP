# Audit GLBFP avant préparation R Journal

Date: 2026-05-02

Référentiel audité: `main`, commit `7a66e63`

Sources consultées:

- CRAN Repository Policy, CRAN Repository Maintainers, consulté le 2026-05-02, https://stat.ethz.ch/CRAN/web/packages/policies.html
- The R Journal, "How to submit your article", consulté le 2026-05-02, https://journal.r-project.org/submissions.html
- rjtools documentation, consulté le 2026-05-02, https://rjournal.github.io/rjtools/
- Scott, D. W. (1992). Multivariate Density Estimation: Theory, Practice, and Visualization. Wiley. doi:10.1002/9780470316849
- Terrell, G. R., and Scott, D. W. (1985). Oversmoothed Nonparametric Density Estimates. Journal of the American Statistical Association, 80(389), 209-214. doi:10.1080/01621459.1985.10477163

## 1. État actuel du package

Le dépôt contient un package R standard avec `DESCRIPTION`, `NAMESPACE`, `R/`,
`man/`, `data/`, `vignettes/`, `README.md`, `NEWS.md`, `LICENSE.md`,
`.Rbuildignore` et un workflow GitHub Actions `R-CMD-check.yaml`.

Fonctions exportées auditées:

- `ASH()`, `ASH_estimate()`
- `LBFP()`, `LBFP_estimate()`
- `GLBFP()`, `GLBFP_estimate()`
- `compute_bi_optim()`
- `compute_G_star()`, `G_i()`, `K_mi()`

Données incluses:

- `ashua`, 4 389 lignes et 3 variables (`flow`, `level`, `day`), taille locale approximative 369 656 octets.

Check initial:

- Commande exécutée: `rcmdcheck::rcmdcheck(args = c("--as-cran"), build_args = c("--no-manual"), error_on = "never")`
- Résultat initial: 0 ERROR, 0 WARNING, 2 NOTE.
- Notes initiales: nouvelle soumission CRAN, validation HTML ignorée localement parce que `tidy` n'est pas assez récent.

## 2. Problèmes bloquants pour CRAN

Le check initial ne révèle pas de ERROR ni de WARNING. Les points suivants restent toutefois à corriger avant une soumission sérieuse:

- `DESCRIPTION` ne contient pas `URL` ni `BugReports`.
- Les exemples sont parfois plus longs que nécessaire pour un check CRAN.
- `tests/` est absent dans l'état initial.
- Les validations d'entrée ne couvrent pas clairement les données non numériques, les valeurs manquantes, les bornes dégénérées, ni les grilles invalides.
- La documentation des données ne fournit pas un script reproductible de création de `ashua`.
- La vignette contient des placeholders tels que `REF_ARTICLE` et `REFERENCE ARTICLE`.
- Le fichier `doc/` et `Meta/` sont présents dans le dépôt, mais exclus du build. Leur présence n'est pas bloquante, mais elle augmente le bruit du dépôt.

## 3. Problèmes bloquants pour The R Journal

- L'article méthodologique GLBFP original n'est pas cité avec une référence bibliographique complète vérifiable dans le dépôt.
- Les vignettes initiales ne suffisent pas à soutenir un article logiciel: elles sont longues, partiellement exploratoires, et contiennent des placeholders.
- Aucun benchmark structuré n'est disponible dans un dossier reproductible séparé.
- Aucun plan de reproductibilité spécifique à The R Journal n'est présent.
- Aucun plan de figures, tableaux, ou positionnement article logiciel n'est présent.
- Le package n'est pas encore accompagné d'une documentation publique de type pkgdown.

## 4. Problèmes de documentation

- Plusieurs pages roxygen décrivent les arguments de façon correcte mais trop générale.
- Les valeurs retournées ne sont pas toujours décrites avec assez de précision.
- Les exemples utilisent souvent `ashua` complet et des grilles relativement grandes.
- Les fonctions internes ne sont pas factorisées, ce qui rend les comportements d'erreur moins cohérents.
- La documentation de `ashua` ne donne pas l'identifiant hydrométrique ni le script de génération.

## 5. Problèmes de tests

- Aucun test `testthat` n'est présent dans l'état initial.
- Les comportements suivants ne sont pas couverts: données 1D, données 2D, valeurs manquantes, données constantes, entrées non numériques, grilles irrégulières, petits échantillons, sorties de tracé, non-négativité des densités, stabilité de résultats reproductibles.

## 6. Problèmes de reproductibilité

- La vignette initiale utilise un exemple de simulation mais ne sépare pas clairement code pédagogique et code de benchmark.
- Les données incluses ne sont pas accompagnées d'un script `data-raw/`.
- Les benchmarks ne sont pas séparés des vignettes.
- Les références méthodologiques ne sont pas toutes vérifiées.

## 7. Problèmes d'API

- Les noms historiques sont conservés et compréhensibles, mais il n'existe pas encore d'interface S3 commune.
- Les classes retournées sont propres aux fonctions (`"GLBFP"`, `"LBFP"`, `"ASH"`, etc.), mais aucun `summary()` ni `predict()` commun n'est disponible.
- Les comportements avec matrices, data frames et grilles irrégulières ne sont pas complètement harmonisés.
- Les messages d'erreur ne sont pas toujours spécifiques.

## 8. Problèmes de performance

- Les fonctions de grille appellent les estimateurs point par point avec `lapply()`. C'est simple et maintenable, mais peut devenir coûteux pour les grandes grilles.
- Les calculs 2D avec `GLBFP_estimate()` peuvent devenir lents lorsque `grid_size` ou `m` augmente.
- Aucun benchmark de temps ou de mémoire n'est disponible dans l'état initial.

## 9. Problèmes de style R

- Le style est globalement lisible, mais la validation est répétée entre fonctions.
- Certaines lignes longues et messages d'erreur pourraient être simplifiés.
- Quelques textes de documentation utilisent une terminologie imprécise ou promotionnelle.
- Les méthodes de tracé 2D ne gèrent pas clairement les grilles irrégulières.

## 10. Priorités

### Critique

- Ajouter des tests unitaires.
- Retirer les placeholders non vérifiés dans les vignettes.
- Ajouter une référence méthodologique vérifiée ou documenter explicitement son absence.
- Renforcer la validation d'entrée.
- Corriger les méthodes de tracé fragiles.
- Ajouter `URL` et `BugReports` dans `DESCRIPTION`.

### Importante

- Refaire le README.
- Ajouter des vignettes courtes, rapides et reproductibles.
- Ajouter des benchmarks séparés de `R CMD check`.
- Ajouter `inst/CITATION`.
- Ajouter un dossier de préparation The R Journal.
- Ajouter pkgdown et coverage GitHub Actions.

### Souhaitable

- Ajouter des méthodes S3 communes `summary()` et `predict()`.
- Améliorer progressivement les performances des évaluations sur grille.
- Ajouter un script complet `data-raw/` lorsque les fichiers sources de `ashua` sont disponibles.
- Ajouter une vérification URL avant soumission CRAN.
