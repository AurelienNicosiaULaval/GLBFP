# Positionnement pour The R Journal

## Type d'article

Article logiciel portant sur un package R.

## Contribution principale

La contribution visée est l'implémentation R, la documentation, les exemples
reproductibles, les tests, les benchmarks et l'intégration dans l'écosystème R
pour les estimateurs ASH, LBFP et GLBFP.

L'article ne doit pas être présenté comme une nouvelle contribution théorique
principale.

## Public visé

- utilisatrices et utilisateurs de R intéressés par l'estimation de densité;
- statisticiens et statisticiennes appliqués;
- personnes préparant des analyses reproductibles ou de l'enseignement sur les
  estimateurs de densité non paramétriques.

## Message central

`GLBFP` fournit une interface R documentée et testée pour appliquer des
estimateurs de densité histogram-based, avec des exemples 1D et 2D, des méthodes
de visualisation, et une base reproductible pour comparer les résultats.

## Limites à déclarer

- Les performances peuvent se dégrader lorsque la dimension, `grid_size`, ou `m`
  augmente.
- Les données manquantes doivent être traitées avant estimation.
- Les exemples de benchmark ne remplacent pas une étude exhaustive de performance.
- La référence méthodologique GLBFP complète doit être vérifiée avant rédaction.
