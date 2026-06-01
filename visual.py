import os
import json
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns

# Set path configurations
OUT_DIR = "./out"
TREE_FILE = os.path.join(OUT_DIR, "sim_tree.txt")
V_FILE = os.path.join(OUT_DIR, "sim_variant_matrix.txt")
D_FILE = os.path.join(OUT_DIR, "sim_total_matrix.txt")
F_FILE = os.path.join(OUT_DIR, "sim_frequency_matrix.txt")
U_FILE = os.path.join(OUT_DIR, "sim_usage_matrix.txt")
JSON_FILE = os.path.join(OUT_DIR, "fastppm_results.json")


def load_matrix(filepath):
    """Loads a space-separated matrix file into a numpy array."""
    if os.path.exists(filepath):
        return np.loadtxt(filepath)
    return None


def parse_tree(filepath):
    """Parses the tree topology into a directed networkx graph."""
    G = nx.DiGraph()
    if not os.path.exists(filepath):
        return G
    with open(filepath, 'r') as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.split()
            parent = int(parts[0])
            G.add_node(parent)
            for child in parts[1:]:
                G.add_edge(parent, int(child))
    return G


def set_boundary_ticks(ax, num_rows, num_cols):
    """Sets ticks only at the beginning (0) and end (max) of the axes."""
    if num_cols > 0:
        ax.set_xticks([0.5, num_cols - 0.5])
        ax.set_xticklabels(['0', str(num_cols - 1)])
    if num_rows > 0:
        ax.set_yticks([0.5, num_rows - 0.5])
        ax.set_yticklabels(['0', str(num_rows - 1)], rotation=0)
    # Hide the little tick line marks for a cleaner look
    ax.tick_params(axis='both', which='both', length=0)


# Load data structures
G = parse_tree(TREE_FILE)
V = load_matrix(V_FILE)
D = load_matrix(D_FILE)
F = load_matrix(F_FILE)
U = load_matrix(U_FILE)

# Calculate Observed VAF (safely handling division by zero if any D=0)
if V is not None and D is not None:
    observed_VAF = np.divide(V, D, out=np.zeros_like(V), where=D != 0)
else:
    observed_VAF = None

# Extract dynamic dimensions to display in the summary
num_samples = V.shape[0] if V is not None else 0
num_mutations = V.shape[1] if V is not None else len(G.nodes)
num_clones = U.shape[1] if U is not None else 0

# Load final optimization scores from JSON
obj_score, runtime = "N/A", "N/A"
if os.path.exists(JSON_FILE):
    with open(JSON_FILE, 'r') as f:
        res = json.load(f)
        obj_score = f"{res.get('objective', 0):.5f}"
        runtime = f"{res.get('runtime', 0):.4f}s"

# Setup the Multi-Panel Figure Layout
fig = plt.figure(figsize=(28, 16))
grid = plt.GridSpec(2, 3, wspace=0.15, hspace=0.3)

# ----------------------------------------------------
# PANEL 1: Evolutionary Clone Tree
# ----------------------------------------------------
ax_tree = fig.add_subplot(grid[0:2, 0])

try:
    from networkx.drawing.nx_pydot import graphviz_layout
    pos = graphviz_layout(G, prog="dot")
except (ImportError, FileNotFoundError, Exception):
    pos = nx.kamada_kawai_layout(G)

# Draw nodes, edges, and small indices inside the nodes
nx.draw_networkx_nodes(G, pos, node_size=540, node_color='#3b82f6',
                       edgecolors='#1e3a8a', linewidths=0.5, ax=ax_tree)
nx.draw_networkx_edges(G, pos, arrows=True, arrowsize=8,
                       width=0.8, edge_color='#64748b', ax=ax_tree)
nx.draw_networkx_labels(G, pos, font_size=14, font_color='white', ax=ax_tree)

ax_tree.set_title("Inferred Clone Tree Structure ($T$)",
                  fontsize=16, fontweight='bold', pad=10)
ax_tree.axis('off')

# ----------------------------------------------------
# PANEL 2: Observed Variant Allele Frequency (VAF)
# ----------------------------------------------------
ax_vaf = fig.add_subplot(grid[0, 1])
if observed_VAF is not None:
    sns.heatmap(observed_VAF, cmap="YlOrRd", vmin=0, vmax=1, ax=ax_vaf,
                xticklabels=False, yticklabels=False, cbar_kws={'label': 'Frequency'})
    set_boundary_ticks(ax_vaf, num_samples, num_mutations)

ax_vaf.set_title("Observed Data: VAF ($V / D$)",
                 fontsize=15, fontweight='bold')
ax_vaf.set_xlabel("Mutations", fontsize=12)
ax_vaf.set_ylabel("Samples", fontsize=12)

# ----------------------------------------------------
# PANEL 3: Latent Frequencies (F)
# ----------------------------------------------------
ax_f = fig.add_subplot(grid[0, 2])
if F is not None:
    sns.heatmap(F, cmap="YlOrRd", vmin=0, vmax=1, ax=ax_f,
                xticklabels=False, yticklabels=False, cbar_kws={'label': 'Frequency'})
    set_boundary_ticks(ax_f, num_samples, num_mutations)

ax_f.set_title("Regression Fit: Latent Frequencies ($F = UB$)",
               fontsize=15, fontweight='bold')
ax_f.set_xlabel("Mutations", fontsize=12)

# ----------------------------------------------------
# PANEL 4: Clone Usage Proportions (U)
# ----------------------------------------------------
ax_u = fig.add_subplot(grid[1, 1:3])
if U is not None:
    sns.heatmap(U, cmap="Purples", vmin=0, ax=ax_u,
                xticklabels=False, yticklabels=False, cbar_kws={'label': 'Proportion'})
    set_boundary_ticks(ax_u, num_samples, num_clones)

ax_u.set_title("Inferred Clone Usage / Proportions Matrix ($U$)",
               fontsize=15, fontweight='bold')
ax_u.set_xlabel("Clones", fontsize=12)
ax_u.set_ylabel("Samples", fontsize=12)

# ----------------------------------------------------
# Summary Box: Final Regression Results
# ----------------------------------------------------
summary_text = (
    f"--- fastppm Run Summary ---\n"
    f"Optimization Loss: {obj_score}\n"
    f"Compute Runtime:   {runtime}\n"
    f"Dimensions:        {num_samples} Samples\n"
    f"                   {num_mutations} Mutations\n"
    f"                   {num_clones} Clones"
)
fig.text(0.13, 0.05, summary_text, fontsize=14, family='monospace',
         bbox=dict(boxstyle='round,pad=0.8', facecolor='#f8fafc', edgecolor='#cbd5e1'))

plt.suptitle("Perfect Phylogeny Mixture Regression Dashboard",
             fontsize=24, fontweight='bold', y=0.95)

# Save dashboard image
output_dashboard = os.path.join(OUT_DIR, "fastppm_visual.png")
plt.savefig(output_dashboard, dpi=300, bbox_inches='tight')
print(f"Visualization generated: {output_dashboard}")
