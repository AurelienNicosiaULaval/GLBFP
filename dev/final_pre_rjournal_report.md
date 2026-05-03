# Rapport final de préparation pré-R Journal

Date: 2026-05-02

## 1. Résumé des modifications

Le package a été préparé pour un état nettement plus proche de CRAN et de The R
Journal:

- métadonnées `DESCRIPTION` complétées avec `URL`, `BugReports`, `Suggests`,
  `Config/testthat/edition` et roxygen markdown;
- validations d'entrée centralisées dans `R/internal_utils.R`;
- messages d'erreur clarifiés pour données non numériques, valeurs manquantes,
  bornes dégénérées, grilles invalides et covariance dégénérée;
- méthodes S3 ajoutées pour `summary()` et `predict()`;
- documentation roxygen enrichie et pages `man/` régénérées;
- tests `testthat` ajoutés;
- README réécrit;
- vignettes rapides et reproductibles ajoutées;
- benchmarks reproductibles séparés de `R CMD check`;
- fichiers de citation ajoutés;
- workflows GitHub Actions coverage et pkgdown ajoutés;
- préparation The R Journal ajoutée dans `rjournal-prep/`;
- audit initial et notes CRAN documentés dans `dev/`.

## 2. Liste des fichiers modifiés ou ajoutés

Métadonnées et configuration:

- `.Rbuildignore`
- `.lintr`
- `DESCRIPTION`
- `NAMESPACE`
- `NEWS.md`
- `_pkgdown.yml`
- `CITATION.cff`
- `cran-comments.md`
- `.github/workflows/R-CMD-check.yaml`
- `.github/workflows/coverage.yaml`
- `.github/workflows/pkgdown.yaml`

Code R:

- `R/ASH.R`
- `R/ASH_estimate.R`
- `R/GLBFP-package.R`
- `R/GLBFP.R`
- `R/GLBFP_estimate.R`
- `R/G_i.R`
- `R/G_star.R`
- `R/K(mi).R`
- `R/LBFP.R`
- `R/LBFP_estimate.R`
- `R/ashua.R`
- `R/compute_bi_optim.R`
- `R/globals.R`
- `R/imports.R`
- `R/internal_utils.R`
- `R/predict.R`
- `R/summary.R`

Documentation:

- `README.md`
- `inst/CITATION`
- pages `man/` régénérées par roxygen2
- `vignettes/GLBFP_estimation.Rmd`
- `vignettes/GLBFP_introduction.Rmd`
- `vignettes/GLBFP_validation.Rmd`
- `vignettes/getting-started.Rmd`
- `vignettes/glbfp-theory-brief.Rmd`
- `vignettes/two-dimensional-density.Rmd`

Tests, benchmarks et préparation article:

- `tests/testthat.R`
- `tests/testthat/test-edge-cases.R`
- `tests/testthat/test-invariants.R`
- `tests/testthat/test-regression.R`
- `tests/testthat/test-s3-api.R`
- `tests/testthat/test-stability.R`
- `benchmarks/README.md`
- `benchmarks/benchmark_glbfp_1d.R`
- `benchmarks/benchmark_glbfp_2d.R`
- `benchmarks/benchmark_loo_if_applicable.R`
- `inst/bench/run_validation.R`
- `inst/bench/sim_scenarios.R`
- `data-raw/README.md`
- `data-raw/prepare_ashua.R`
- `dev/audit_glbfp_pre_rjournal.md`
- `dev/cran_check_notes.md`
- `dev/references_to_verify.md`
- `rjournal-prep/rjournal_package_positioning.md`
- `rjournal-prep/rjournal_article_outline.md`
- `rjournal-prep/rjournal_reproducibility_plan.md`
- `rjournal-prep/rjournal_figures_and_tables_plan.md`
- `rjournal-prep/rjournal_checklist.md`

## 3. Résultats de R CMD check

Commande finale:

```r
rcmdcheck::rcmdcheck(
  args = c("--as-cran"),
  build_args = c("--no-manual"),
  error_on = "never"
)
```

Résultat final:

- 0 ERROR
- 0 WARNING
- 2 NOTE

Notes:

- `New submission`, attendu pour une première soumission CRAN.
- Validation HTML du manuel ignorée localement parce que HTML Tidy n'est pas
  assez récent dans l'environnement local.

Commande complémentaire:

```r
devtools::check(args = c("--as-cran"), manual = FALSE, error_on = "never")
```

Résultat:

- 0 ERROR
- 0 WARNING
- 0 NOTE

## 4. Résultats des tests

Commande:

```r
devtools::test()
```

Résultat:

- 50 tests passés;
- 0 échec;
- 0 warning;
- 0 skip.

Les tests couvrent les fonctions principales, données 1D et 2D, valeurs
manquantes, données constantes avec bornes explicites, valeurs extrêmes,
entrées invalides, sorties S3, grilles régulières et irrégulières, non-négativité
des densités, et stabilité de régression.

## 5. État des vignettes

Vignettes présentes:

- `getting-started.Rmd`
- `glbfp-theory-brief.Rmd`
- `two-dimensional-density.Rmd`
- `GLBFP_introduction.Rmd`
- `GLBFP_validation.Rmd`
- `GLBFP_estimation.Rmd`

Toutes les vignettes sont reconstruites avec succès pendant `R CMD check`.

La vignette leave-one-out n'a pas été ajoutée parce qu'aucune API leave-one-out
n'est actuellement exportée par le package. Le fichier
`benchmarks/benchmark_loo_if_applicable.R` documente explicitement cette absence.

## 6. État du README

Le README a été réécrit avec:

- statut du package;
- installation GitHub et installation CRAN future;
- exemple minimal;
- exemple 2D avec `ashua`;
- tableau des fonctions principales;
- vignettes;
- références générales vérifiées;
- citation du package;
- contribution;
- licence.

## 7. État de la documentation

La documentation roxygen a été régénérée avec roxygen2 7.3.3. Les fonctions
exportées disposent de titres, descriptions, paramètres, valeurs retournées,
exemples exécutables, détails courts et liens croisés.

La référence méthodologique GLBFP complète n'a pas été inventée. Elle est suivie
dans `dev/references_to_verify.md`.

## 8. État CRAN-readiness

État actuel: presque prêt pour une première soumission, sous réserve des points
restants ci-dessous.

Points positifs:

- `R CMD check --as-cran` sans ERROR ni WARNING;
- tests unitaires présents;
- exemples et vignettes rapides;
- URLs vérifiées avec `urlchecker::url_check()`;
- `lintr::lint_package()` sans lint;
- benchmarks exclus du build;
- `data-raw/`, `dev/`, `benchmarks/` et `rjournal-prep/` exclus du build.

Commandes complémentaires:

- `urlchecker::url_check()`: toutes les URLs sont correctes.
- `lintr::lint_package()`: aucun lint.
- `rhub::rhub_check()`: non exécuté, package `rhub` non disponible localement.
- `revdepcheck`: non exécuté, aucun reverse dependency pertinent pour ce package
  en développement.

## 9. État R Journal-readiness

État actuel: bonne base pour rédiger un article logiciel, mais le manuscrit n'est
pas encore écrit.

Ajouts réalisés:

- positionnement logiciel dans `rjournal-prep/rjournal_package_positioning.md`;
- plan d'article dans `rjournal-prep/rjournal_article_outline.md`;
- plan de reproductibilité;
- plan des figures et tableaux;
- checklist de soumission;
- benchmarks reproductibles;
- vignettes directement réutilisables comme base d'exemples.

## 10. Problèmes restants

- La référence complète de l'article méthodologique original GLBFP n'est pas
  vérifiée.
- La provenance complète de `ashua` doit être documentée dans `data-raw/`.
- Le package n'a pas encore été testé sur GitHub Actions après ces modifications.
- `rhub` n'est pas disponible localement.
- Les performances des fonctions de grille restent simples et robustes, mais pas
  encore fortement optimisées.

## 11. Recommandations avant soumission CRAN

- Vérifier et ajouter la référence méthodologique GLBFP complète.
- Compléter `data-raw/prepare_ashua.R` avec la source brute exacte.
- Lancer GitHub Actions sur macOS, Windows et Ubuntu.
- Exécuter un check externe avec R-hub ou une infrastructure équivalente.
- Relire `cran-comments.md` juste avant soumission.
- Vérifier que la version `0.5.1` correspond bien à la stratégie de release.

## 12. Recommandations avant rédaction de l'article R Journal

- Créer le manuscrit avec `rjtools`.
- Utiliser les vignettes comme base de code reproductible.
- Produire les tableaux à partir de `benchmarks/` ou `inst/bench/`.
- Ajouter du texte alternatif aux figures.
- Garder l'article centré sur la contribution logicielle.
- Ne pas présenter le package comme une nouvelle contribution théorique.
- Exécuter `rjtools::initial_check_article()` lorsque le manuscrit existe.
