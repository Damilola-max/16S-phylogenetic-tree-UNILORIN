#!/usr/bin/env python3
"""
Phylogenetic Tree Builder with MAFFT Alignment and Circular Visualization
Aligns 16S rRNA sequences and generates a circular Neighbor-Joining tree
"""

import subprocess
import sys
from pathlib import Path
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Phylo import write as phylo_write
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np

# File paths
INPUT_FASTA = "/Users/user/Documents/Grit/reference_genes/phylogenetic_tree_sequences_supplementary.fasta"
ALIGNED_FASTA = "/Users/user/Documents/Grit/reference_genes/ALIGNED_sequences.fasta"
TREE_NEWICK = "/Users/user/Documents/Grit/reference_genes/phylogenetic_tree.nwk"
CIRCULAR_PNG = "/Users/user/Documents/Grit/reference_genes/circular_phylogenetic_tree.png"
CIRCULAR_SVG = "/Users/user/Documents/Grit/reference_genes/circular_phylogenetic_tree.svg"

def run_mafft_alignment():
    """Run MAFFT alignment on input sequences"""
    print("=" * 70)
    print("STEP 1: Running MAFFT Multiple Sequence Alignment")
    print("=" * 70)
    
    cmd = [
        "mafft",
        "--auto",  # Automatic algorithm selection
        "--maxiterate", "1000",  # High accuracy
        "--thread", "4",  # Use 4 threads
        INPUT_FASTA
    ]
    
    print(f"Command: {' '.join(cmd)}")
    print(f"Input: {INPUT_FASTA}")
    print(f"Output: {ALIGNED_FASTA}")
    print("-" * 70)
    
    try:
        with open(ALIGNED_FASTA, 'w') as outfile:
            result = subprocess.run(cmd, stdout=outfile, stderr=subprocess.PIPE, text=True)
        
        if result.returncode == 0:
            print("✓ MAFFT alignment completed successfully!")
            
            # Count sequences in alignment
            with open(ALIGNED_FASTA) as f:
                content = f.read()
                seq_count = content.count('>')
            print(f"✓ Aligned {seq_count} sequences")
            return True
        else:
            print(f"✗ MAFFT failed with return code {result.returncode}")
            print(f"Error: {result.stderr}")
            return False
            
    except Exception as e:
        print(f"✗ Error running MAFFT: {e}")
        return False

def build_neighbor_joining_tree():
    """Build Neighbor-Joining tree from alignment"""
    print("\n" + "=" * 70)
    print("STEP 2: Building Neighbor-Joining Phylogenetic Tree")
    print("=" * 70)
    
    try:
        # Read alignment
        print(f"Reading alignment from: {ALIGNED_FASTA}")
        alignment = AlignIO.read(ALIGNED_FASTA, "fasta")
        print(f"✓ Loaded alignment with {len(alignment)} sequences")
        print(f"✓ Alignment length: {alignment.get_alignment_length()} bp")
        
        # Calculate distance matrix
        print("\nCalculating distance matrix (Kimura 2-parameter)...")
        calculator = DistanceCalculator('k2p')  # Kimura 2-parameter
        dm = calculator.get_distance(alignment)
        print("✓ Distance matrix calculated")
        
        # Build NJ tree
        print("\nConstructing Neighbor-Joining tree...")
        constructor = DistanceTreeConstructor()
        tree = constructor.nj(dm)
        print("✓ NJ tree constructed")
        
        # Root the tree at midpoint
        tree.root_at_midpoint()
        print("✓ Tree rooted at midpoint")
        
        # Save tree in Newick format
        phylo_write(tree, TREE_NEWICK, "newick")
        print(f"✓ Tree saved to: {TREE_NEWICK}")
        
        return tree, dm
        
    except Exception as e:
        print(f"✗ Error building tree: {e}")
        import traceback
        traceback.print_exc()
        return None, None

def create_circular_tree_visualization(tree):
    """Create a circular (radial) tree visualization"""
    print("\n" + "=" * 70)
    print("STEP 3: Creating Circular Tree Visualization")
    print("=" * 70)
    
    try:
        fig, ax = plt.subplots(figsize=(14, 14), subplot_kw=dict(projection='polar'))
        
        # Get all terminal nodes (leaves)
        terminals = list(tree.get_terminals())
        n_leaves = len(terminals)
        
        print(f"Total taxa: {n_leaves}")
        
        # Assign angles to leaves (evenly spaced around circle)
        leaf_angles = {}
        leaf_distances = {}
        
        for i, leaf in enumerate(terminals):
            angle = 2 * np.pi * i / n_leaves
            leaf_angles[leaf.name] = angle
            leaf_distances[leaf.name] = leaf.branch_length or 0
        
        # Function to get node position
        def get_node_position(node):
            if node.is_terminal():
                return leaf_angles[node.name], leaf_distances.get(node.name, 0)
            
            # For internal nodes, average the angles of children
            child_angles = []
            child_distances = []
            for child in node.clades:
                a, d = get_node_position(child)
                child_angles.append(a)
                child_distances.append(d + (child.branch_length or 0))
            
            # Handle circular mean
            avg_angle = np.arctan2(
                np.mean([np.sin(a) for a in child_angles]),
                np.mean([np.cos(a) for a in child_angles])
            )
            avg_dist = np.mean(child_distances)
            
            return avg_angle, avg_dist
        
        # Draw tree branches
        def draw_clade(clade, parent_angle=None, parent_radius=0):
            angle, radius = get_node_position(clade)
            
            if parent_angle is not None:
                # Draw branch from parent to current node
                # Use arc for curved branches
                r_start = parent_radius
                r_end = radius
                
                # Draw line
                ax.plot([parent_angle, angle], [r_start, r_end], 
                       'b-', linewidth=1.5, alpha=0.7)
            
            # Draw node
            if clade.is_terminal():
                # Color code: Study = red, Reference = blue
                color = '#E74C3C' if 'Study_' in clade.name else '#3498DB'
                marker_size = 80 if 'Study_' in clade.name else 60
                
                ax.scatter(angle, radius, c=color, s=marker_size, 
                          edgecolors='black', linewidths=1, zorder=5)
                
                # Add label
                label = clade.name.replace('Study_', '').replace('Reference_', '')
                label = label.replace('_16S', '').replace('_', ' ')
                
                # Position label outside the circle
                label_radius = radius + 0.15
                ha = 'left' if -np.pi/2 < angle < np.pi/2 else 'right'
                
                ax.text(angle, label_radius, label, 
                       fontsize=8, ha=ha, va='center',
                       fontweight='bold' if 'Study_' in clade.name else 'normal',
                       rotation=np.degrees(angle) - 90 if -np.pi/2 < angle < np.pi/2 else np.degrees(angle) + 90,
                       rotation_mode='anchor')
            
            # Recurse to children
            for child in clade.clades:
                draw_clade(child, angle, radius)
        
        # Draw the tree
        root = tree.root
        draw_clade(root)
        
        # Configure plot
        ax.set_ylim(0, max(leaf_distances.values()) * 1.3)
        ax.set_yticks([])
        ax.set_xticks([])
        ax.spines['polar'].set_visible(False)
        ax.grid(False)
        
        # Add title
        plt.title('Circular Neighbor-Joining Phylogenetic Tree\n16S rRNA Gene Sequences', 
                 fontsize=14, fontweight='bold', pad=20)
        
        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#E74C3C', edgecolor='black', label='Study Isolates'),
            Patch(facecolor='#3498DB', edgecolor='black', label='Reference Strains')
        ]
        ax.legend(handles=legend_elements, loc='upper right', 
                 bbox_to_anchor=(1.3, 1.1), fontsize=10)
        
        # Save figure
        plt.tight_layout()
        plt.savefig(CIRCULAR_PNG, dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        plt.savefig(CIRCULAR_SVG, format='svg', bbox_inches='tight',
                   facecolor='white', edgecolor='none')
        
        print(f"✓ Circular tree saved to: {CIRCULAR_PNG}")
        print(f"✓ Vector version saved to: {CIRCULAR_SVG}")
        
        plt.close()
        return True
        
    except Exception as e:
        print(f"✗ Error creating visualization: {e}")
        import traceback
        traceback.print_exc()
        return False

def print_summary(tree, dm):
    """Print summary of the analysis"""
    print("\n" + "=" * 70)
    print("ANALYSIS SUMMARY")
    print("=" * 70)
    
    # List all sequences
    print("\nSequences in analysis:")
    print("-" * 70)
    
    terminals = sorted(tree.get_terminals(), key=lambda x: x.name)
    for i, leaf in enumerate(terminals, 1):
        seq_type = "STUDY" if 'Study_' in leaf.name else "REF "
        short_name = leaf.name.replace('Study_', '').replace('Reference_', '').replace('_16S', '')
        print(f"{i:2d}. [{seq_type}] {short_name}")
    
    print("-" * 70)
    print(f"\nTotal sequences: {len(terminals)}")
    print(f"Study isolates: {sum(1 for t in terminals if 'Study_' in t.name)}")
    print(f"Reference strains: {sum(1 for t in terminals if 'Reference_' in t.name)}")
    
    print("\nOutput files:")
    print(f"  • Alignment: {ALIGNED_FASTA}")
    print(f"  • Tree (Newick): {TREE_NEWICK}")
    print(f"  • Tree image (PNG): {CIRCULAR_PNG}")
    print(f"  • Tree image (SVG): {CIRCULAR_SVG}")
    
    print("\n" + "=" * 70)
    print("✓ ANALYSIS COMPLETE")
    print("=" * 70)

def main():
    print("\n" + "=" * 70)
    print("  PHYLOGENETIC TREE BUILDER - 16S rRNA ANALYSIS")
    print("=" * 70)
    print()
    
    # Check input file exists
    if not Path(INPUT_FASTA).exists():
        print(f"✗ Error: Input file not found: {INPUT_FASTA}")
        sys.exit(1)
    
    # Step 1: Run MAFFT alignment
    if not run_mafft_alignment():
        print("\n✗ Alignment failed. Exiting.")
        sys.exit(1)
    
    # Step 2: Build NJ tree
    tree, dm = build_neighbor_joining_tree()
    if tree is None:
        print("\n✗ Tree building failed. Exiting.")
        sys.exit(1)
    
    # Step 3: Create circular visualization
    if not create_circular_tree_visualization(tree):
        print("\n⚠ Warning: Visualization failed, but tree was built successfully.")
    
    # Print summary
    print_summary(tree, dm)

if __name__ == "__main__":
    main()
