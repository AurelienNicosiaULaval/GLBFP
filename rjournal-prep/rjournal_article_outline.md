# Plan d'article R Journal

## 1. Introduction

- Motiver l'intérêt logiciel: disposer d'une implémentation R reproductible pour
  ASH, LBFP et GLBFP.
- Positionner l'article comme contribution logicielle.
- Décrire brièvement le contenu du package.

## 2. Background minimal

- Rappeler le lien entre histogrammes, frequency polygons et estimations lissées.
- Présenter ASH, LBFP et GLBFP au niveau nécessaire pour comprendre l'API.
- Citer les références vérifiées.
- Ajouter la référence méthodologique GLBFP complète une fois vérifiée.

## 3. Package design

- Structure du package.
- Types d'entrées acceptés.
- Objets retournés.
- Méthodes S3 disponibles.
- Gestion des erreurs et cas limites.

## 4. Core functions

- Fonctions pointwise: `ASH()`, `LBFP()`, `GLBFP()`.
- Fonctions de grille: `ASH_estimate()`, `LBFP_estimate()`, `GLBFP_estimate()`.
- Aide au choix de bande passante: `compute_bi_optim()`.

## 5. Examples

- Exemple 1D simulé.
- Exemple 2D avec `ashua`.
- Visualisations reproductibles.
- Interprétation prudente des sorties.

## 6. Computational performance

- Benchmarks 1D et 2D.
- Effet de `grid_size`.
- Effet de `m`.
- Comparaison de temps avec KDE de base seulement si l'objectif est bien défini.

## 7. Reproducibility

- Versions R et packages.
- Données et scripts.
- Tests unitaires.
- GitHub Actions.
- Instructions pour régénérer figures et tableaux.

## 8. Conclusion

- Résumer la contribution logicielle.
- Mentionner les limites.
- Indiquer les améliorations prévues.
