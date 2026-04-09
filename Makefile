# =============================================================================
# Root Makefile — orchestrates build, run, and development workflows
# =============================================================================

PYTHON      ?= python3
UV          ?= uv
FORTRAN_DIR  = src/fortran
SCRIPTS_DIR  = scripts/python
OUTPUT_DIR   = output

# --- High-level targets ---

.PHONY: all build run clean install lint format check help

help: ## Show this help message
	@grep -E '^[a-zA-Z_-]+:.*?##' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

all: build ## Build everything (Fortran solver)

# --- Fortran ---

build: ## Compile the Fortran solver
	$(MAKE) -C $(FORTRAN_DIR) all

run: build ## Build and run the simulation
	@mkdir -p $(OUTPUT_DIR)
	cd $(FORTRAN_DIR) && ./fkrplk.x

clean: ## Remove all build artifacts
	$(MAKE) -C $(FORTRAN_DIR) clean

# --- Python ---

install: ## Install Python deps with uv
	$(UV) sync

install-dev: ## Install Python deps including dev tools
	$(UV) sync --extra dev

lint: ## Run ruff linter on Python code
	$(UV) run ruff check $(SCRIPTS_DIR)

format: ## Auto-format Python code with ruff
	$(UV) run ruff format $(SCRIPTS_DIR)
	$(UV) run ruff check --fix $(SCRIPTS_DIR)

check: lint ## Run all checks (lint + type-check)
	$(UV) run mypy $(SCRIPTS_DIR)

# --- Plotting ---

plot: ## Generate standard plots from latest output
	$(UV) run $(PYTHON) $(SCRIPTS_DIR)/tdfp_plots.py \
		$(word 1, $(wildcard experiments/conf2_E10keV_10MeV/fkrplk.* $(OUTPUT_DIR)/fkrplk.*)) \
		-o $(OUTPUT_DIR)/plots
