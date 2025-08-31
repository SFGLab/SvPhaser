# Contributing to SvPhaser

Thanks for your interest in improving SvPhaser! 🎉  
This doc explains how to set up a dev environment, our coding standards, and how to submit changes.

## Table of contents
- [Development setup](#development-setup)
- [Project layout](#project-layout)
- [Coding standards](#coding-standards)
- [Testing](#testing)
- [Performance checks](#performance-checks)
- [Git & PR workflow](#git--pr-workflow)
- [Issue reports & feature requests](#issue-reports--feature-requests)
- [Security](#security)

## Development setup

```bash
# 1) Create a virtual environment (Python 3.10+ recommended)
python -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate

# 2) Install runtime + dev deps
pip install -U pip
pip install -e .        # installs svphaser in editable mode
pip install -r requirements-dev.txt

# (optional) Install pre-commit hooks
pre-commit install


## Project Layout

```bash
src/svphaser/
  phasing/
    algorithms.py   # math & classification (pure)
    _workers.py     # per-chromosome worker logic
    io.py           # orchestration + CSV/VCF writers
    types.py        # small shared types
    __init__.py     # public API exports


## Coding standards

- Type hints everywhere (PEP 484). Run mypy.
- Style: black for formatting, ruff for linting/imports.
- Docstrings: short, precise, with argument and return descriptions when non-trivial.
- Logging: use the module logger (e.g., logging.getLogger("svphaser.io")), no prints.