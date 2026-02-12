.PHONY: all data rjournal jss

all: data rjournal jss

data:
	Rscript reproduce.R

rjournal:
	quarto render rjournal_outline.qmd

jss:
	quarto render jss_outline.qmd
