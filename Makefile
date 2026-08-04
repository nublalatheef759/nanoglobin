.PHONY: install kernel test evaluation readme review thesis \
	public-hbb-smoke-fetch public-hbb-smoke-dry-run public-hbb-smoke-run clean

PUBLIC_HBB_RUNS ?= SRR37686269 SRR37686273
PUBLIC_HBB_MAX_SPOTS ?= 5000
PUBLIC_HBB_CORES ?= 8

install:
	python -m pip install -e ".[review]"

kernel:
	python setup.py build_ext --inplace

test: kernel
	python -m pytest -q

evaluation:
	python scripts/check_evaluation.py

readme:
	quarto render README.qmd

review: test evaluation readme

thesis:
	quarto render thesis --to html
	quarto render thesis --to typst

# Download bounded prefixes from the two PRJNA1439314 runs configured in
# config.yml. `prefetch` and `fasterq-dump` are provided by envs/public_hbb.yaml.
public-hbb-smoke-fetch:
	mkdir -p fastq
	@for run in $(PUBLIC_HBB_RUNS); do \
	  if [ -s "fastq/$$run.fastq" ]; then \
	    echo "fastq/$$run.fastq already exists"; \
	    continue; \
	  fi; \
	  echo "Fetching $$run (first $(PUBLIC_HBB_MAX_SPOTS) spots)"; \
	  prefetch --max-size 100G "$$run"; \
	  fasterq-dump -X $(PUBLIC_HBB_MAX_SPOTS) -e 4 -O fastq "$$run"; \
	done

public-hbb-smoke-dry-run:
	snakemake --use-conda --cores 1 --dry-run --rerun-incomplete

public-hbb-smoke-run:
	snakemake --use-conda --cores $(PUBLIC_HBB_CORES) \
	  --printshellcmds --rerun-incomplete

clean:
	rm -rf README_files thesis/_book thesis/.quarto .pytest_cache build \
	  nanoglobin/_genotype_kernel.c nanoglobin/_genotype_kernel*.so
