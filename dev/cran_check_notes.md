# Notes R CMD check

## Check initial du 2026-05-02

Commande:

```r
rcmdcheck::rcmdcheck(
  args = c("--as-cran"),
  build_args = c("--no-manual"),
  error_on = "never"
)
```

Résultat:

- 0 ERROR
- 0 WARNING
- 2 NOTE

Notes:

- `New submission`: attendu pour un package non encore soumis à CRAN.
- HTML manual validation skipped: dépend de l'installation locale de HTML Tidy et ne reflète pas un problème du package.

## Check final du 2026-05-02

Commande:

```r
rcmdcheck::rcmdcheck(
  args = c("--as-cran"),
  build_args = c("--no-manual"),
  error_on = "never"
)
```

Résultat:

- 0 ERROR
- 0 WARNING
- 2 NOTE

Notes:

- `New submission`: attendu pour une première soumission CRAN.
- HTML manual validation skipped: dépend de l'installation locale de HTML Tidy.

Ces deux notes ne nécessitent pas de modification du package.
