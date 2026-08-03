.PHONY: install test evaluation readme review thesis clean

install:
	python -m pip install -e ".[review]"

test:
	python -m pytest

evaluation:
	python scripts/check_evaluation.py

readme:
	quarto render README.qmd

# README.qmd executes the tests and evaluation checker itself.
review: install readme

thesis:
	quarto render thesis --to html
	quarto render thesis --to typst

clean:
	rm -rf README_files thesis/_book thesis/.quarto .pytest_cache
