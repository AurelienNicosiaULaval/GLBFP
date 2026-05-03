# Plan de reproductibilité

## Environnement

- Construire l'article avec `rjtools`, pas `rticles::rjournal_article`.
- Enregistrer `sessionInfo()` dans le matériel supplémentaire.
- Utiliser des graines fixes pour toutes les simulations.
- Éviter tout accès internet pendant la construction de l'article.

## Données

- Utiliser les données incluses `ashua` ou des données simulées.
- Compléter `data-raw/` avant soumission.
- Documenter la provenance de `ashua`.

## Code

- Les exemples principaux doivent provenir des vignettes ou de scripts dans
  `benchmarks/`.
- Les benchmarks lourds ne doivent pas tourner dans `R CMD check`.
- Les figures doivent être générées par code R évalué.

## Validation

- Exécuter `devtools::document()`.
- Exécuter `devtools::test()`.
- Exécuter `devtools::check()`.
- Exécuter `rcmdcheck::rcmdcheck(args = c("--as-cran"))`.
- Exécuter `rjtools::initial_check_article()` lorsque le manuscrit existe.
