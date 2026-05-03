# Checklist R Journal

## Package

- `R CMD check --as-cran` sans ERROR ni WARNING.
- Notes CRAN expliquées.
- Tests unitaires présents.
- Vignettes rapides.
- README à jour.
- `inst/CITATION` à jour.
- Site pkgdown préparé.

## Article

- Article créé avec `rjtools`.
- Moins de 20 pages, sauf justification éditoriale.
- Résumé de moins de 250 mots.
- Code et données reproductibles.
- Figures et tableaux générés par code.
- Bibliographie limitée aux références citées.
- Texte alternatif pour les figures.
- Référence méthodologique GLBFP complète vérifiée.

## Reproductibilité

- Scripts de benchmark séparés.
- Pas d'accès internet pendant la compilation.
- Graines fixées.
- Instructions de reconstruction incluses.
- Résultats régénérables en temps raisonnable.

## Avant soumission

- Vérifier URLs.
- Vérifier orthographe.
- Exécuter `rjtools::initial_check_article()`.
- Préparer lettre de motivation.
- Préparer archive de soumission complète.
