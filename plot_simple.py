#!/usr/bin/env python3
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from Bio import Phylo

tree = Phylo.read('STEP3_tree.nwk', 'newick')
leaves = list(tree.get_terminals())
n = len(leaves)

print(f'Creating circular tree with {n} sequences...')

fig, ax = plt.subplots(figsize=(16, 16), subplot_kw=dict(projection='polar'))

for i, leaf in enumerate(leaves):
    path = tree.get_path(leaf)
    dist = sum(node.branch_length for node in path if node.branch_length)
    angle = 2 * np.pi * i / n
    study = 'Study_' in leaf.name
    color = '#E74C3C' if study else '#3498DB'
    
    ax.plot([0, angle], [0, dist], 'b-', linewidth=2)
    ax.scatter(angle, dist, c=color, s=200, edgecolors='black', linewidths=2, zorder=5)
    
    label = leaf.name.replace('Study_', '').replace('Reference_', '').replace('_16S', '').replace('_', ' ')[:30]
    r = dist + 0.05
    rotation = np.degrees(angle) - 90
    if not (-90 < rotation < 90):
        rotation += 180
    ax.text(angle, r, label, fontsize=10, ha='center', va='center', 
           fontweight='bold' if study else 'normal', rotation=rotation, rotation_mode='anchor')

ax.set_ylim(0, max([sum(node.branch_length for node in tree.get_path(l) if node.branch_length) for l in leaves]) * 1.5)
ax.set_yticks([])
ax.set_xticks([])
ax.spines['polar'].set_visible(False)
ax.grid(False)
plt.title('Circular NJ Tree - 16S rRNA (10 Taxa)', fontsize=16, fontweight='bold', pad=20)

from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#E74C3C', edgecolor='black', label='Study Isolates (5)'),
    Patch(facecolor='#3498DB', edgecolor='black', label='Reference Strains (5)')
]
ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.3, 1.05), fontsize=12)

plt.savefig('FINAL_CIRCULAR_TREE.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('FINAL_CIRCULAR_TREE.svg', format='svg', bbox_inches='tight', facecolor='white')
print('SAVED: FINAL_CIRCULAR_TREE.png and .svg')

for i, leaf in enumerate(leaves, 1):
    t = 'S' if 'Study_' in leaf.name else 'R'
    print(f'{i}. [{t}] {leaf.name[:40]}')
