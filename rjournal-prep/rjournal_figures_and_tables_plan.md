# Plan des figures et tableaux

## Figures possibles

1. Exemple 1D simulé
   - données simulées normales;
   - courbe GLBFP estimée;
   - comparaison visuelle avec `stats::density()` si utile.

2. Exemple 2D
   - contour GLBFP sur `ashua` ou sur simulation contrôlée;
   - axes clairement libellés;
   - figure statique reproductible.

3. Sensibilité
   - effet de différentes valeurs de `m`;
   - effet d'un facteur multiplicatif sur `b`.

## Tableaux possibles

1. Résumé des fonctions exportées
   - fonction;
   - rôle;
   - entrée principale;
   - sortie principale.

2. Benchmarks 1D
   - taille d'échantillon;
   - taille de grille;
   - méthode;
   - temps médian;
   - mémoire si mesurée avec `bench`.

3. Benchmarks 2D
   - même structure que le tableau 1D.

4. Résultats de validation
   - scénario;
   - méthode;
   - erreur intégrée approximative;
   - temps.

## Accessibilité

- Ajouter du texte alternatif aux figures dans le manuscrit R Journal.
- Préférer des figures statiques pour la version PDF.
- Garder les objets interactifs optionnels dans la version HTML seulement.
