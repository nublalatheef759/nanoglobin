.PHONY: preview html typst all clean

preview:
	quarto preview

html:
	quarto render --to html

typst:
	quarto render --to typst

all: html typst

clean:
	rm -rf _book .quarto
