# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository is a cancer phylogenetics research toolkit for simulating clonal evolution data and running the fastppm (Perfect Phylogeny Mixture) algorithm to infer clone trees from sequencing data.

**Core workflow:** Simulate data → Run fastppm → Visualize results

## Common Commands

### Run full pipeline (simulate + fastppm + visualize)
```bash
python run.py -m 100 -s 10 -c 100
```
- `-m/--mutations`: Number of mutations (creates that many clones by default)
- `-s/--samples`: Number of sequenced samples
- `-c/--coverage`: Expected sequencing coverage
- `-l/--loss`: Loss function for fastppm (default: `l2`; can also use binomial-based losses)
- `-r/--seed`: Random seed (optional)

### Run just the visualization (after a run)
```bash
python visual.py
```
Output is saved to `out/fastppm_visual.png`

### Simulate data only
```bash
python scripts/simulate.py --mutations 100 --samples 10 --coverage 100 --output sim/test --seed 42
```

### Run fastppm directly (on existing simulated data)
```bash
./fastppm-cli -v sim/sim_variant_matrix.txt -d sim/sim_total_matrix.txt -t sim/sim_tree.txt -o out/results.json -f verbose -l l2
```

## Architecture

### Core Components

**`fastppm-cli`**: Pre-compiled binary that implements the Perfect Phylogeny Mixture regression algorithm. Takes variant counts, total depth, and a tree topology; outputs inferred usage matrix and objective score.

**`run.py`**: Orchestration script that:
1. Cleans `./sim` and `./out` directories
2. Runs `scripts/simulate.py` to generate synthetic data
3. Runs `fastppm-cli` on the simulated data
4. Copies results to `./out`
5. Runs `visual.py` for visualization

**`visual.py`**: Creates a 4-panel dashboard showing:
- Clone tree structure (networkx + graphviz layout)
- Observed VAF heatmap
- Latent frequencies (F = UB)
- Clone usage proportions (U)

### Simulation Model (`scripts/simulate.py`)

Simulates cancer phylogeny by:
1. Generating a random spanning tree as the clonal tree (root=0)
2. Assigning mutations to tree nodes (each node gets ≥1 mutation)
3. Creating a usage matrix (random subset of clones per sample, Dirichlet-distributed proportions)
4. Computing latent frequencies: F = U @ B (B = clonal matrix, B[j,i]=1 if path exists from i→j)
5. Simulating sequencing reads: Poisson(coverage) total counts, then Binomial(f) variant counts

Output files (all in `./sim/` with `sim_` prefix):
- `clonal_matrix.txt`: Binary matrix B
- `usage_matrix.txt`: Sample-clone proportions U
- `frequency_matrix.txt`: Observed clone-level frequencies
- `variant_matrix.txt`: Variant read counts (collapsed to clone level)
- `total_matrix.txt`: Total read counts
- `tree.txt`: Adjacency list format of true tree
- `citup_tree.txt`: CITUP-format tree string

### Evaluation Scripts

- `scripts/nextflow/*.nf`: Nextflow workflows for running algorithm comparisons (citup, orchard, sapling) on HPC clusters
- `scripts/reference_regression.py`: CVXPY-based reference implementation for comparison
- `scripts/summarize_evaluations.py`: Computes F1 scores, ancestor-descendant distances between inferred and true trees
- `scripts/edge_diff.py`: Utility to compute symmetric differences between tree edge sets

## Data Flow

```
┌─────────────────┐     ┌──────────────────┐     ┌─────────────────┐
│ simulate.py     │────→│ fastppm-cli      │────→│ visual.py       │
│ (creates ./sim/)│     │ (reads ./sim/)   │     │ (reads ./out/)  │
└─────────────────┘     │ (writes ./out/)  │     └─────────────────┘
                        └──────────────────┘
```

After `run.py` completes, `./out/` contains:
- `fastppm_results.json`: fastppm output with objective score and runtime
- `fastppm_visual.png`: Visualization dashboard
- `sim_*.txt`: Copied simulation inputs for reference

## Key Dependencies

Python:
- `numpy`, `pandas`: Matrix operations
- `networkx`: Tree representation and algorithms
- `matplotlib`, `seaborn`: Visualization
- `scipy`: Statistical distributions (Poisson, Binomial)

For Nextflow workflows (scripts/nextflow/):
- Requires access to institutional HPC paths (hardcoded in params)
- CVXPY with solvers (ECOS, CLARABEL, MOSEK) for reference regression
