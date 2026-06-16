#!/usr/bin/env python3
from Bio import Phylo
import matplotlib.pyplot as plt
import numpy as np

tree = Phylo.read('STEP3_tree.nwk', 'newick')
leaves = list(tree.get_terminals())
n = len(leaves)

leaf_data = []
for leaf in leaves:
    path = tree.get_path(leaf)
    dist = sum(node.branch_length for node in path if node.branch_length)
    leaf_data.append((leaf.name, dist))
leaf_data.sort(key=lambda x: x[1])
max_dist = max(d for _, d in leaf_data)

fig, ax = plt.subplots(figsize=(14, 14), subplot_kw=dict(projection='polar'))

for i, (name, dist) in enumerate(leaf_data):
    angle = 2 * np.pi * i / n
    study = 'Study_' in name
    ax.plot([0, angle], [0, dist], 'b-', lw=2)
    ax.scatter(angle, dist, c='red' if study else 'blue', s=120, edgecolors='black', zorder=5)
    
    label = name.replace('Study_', '').replace('Reference_', '').replace('_16S', '').replace('_', ' ')[:35]
    r = dist + max_dist * 0.2
    rot = np.degrees(angle) - 90
    if not (-90 < rot < 90):
        rot += 180
    ax.text(angle, r, label, fontsize=9, ha='center', va='center', 
           fontweight='bold' if study else 'normal', rotation=rot, rotation_mode='anchor')

ax.set_ylim(0, max_dist * 1.4)
ax.set_yticks([])
ax.set_xticks([])
ax.spines['polar'].set_visible(False)
plt.title('Circular NJ Phylogenetic Tree - 16S rRNA', fontsize=14, fontweight='bold', pad=20)

plt.savefig('circular_tree.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('circular_tree.svg', format='svg', bbox_inches='tight', facecolor='white')

print('✓ circular_tree.png created')
print('✓ circular_tree.svg created')
print(f'\nTree has {n} sequences:')
for name, _ in leaf_data:
    t = 'S' if 'Study_' in name else 'R'
    print(f'  [{t}] {name[:40]}')
