from Bio import Phylo, SeqIO
import sys
import os
import re

# Check input arguments
if len(sys.argv) != 3:
    print("Usage: python pruning_tree.py <fasta_file> <tree_file>")
    sys.exit(1)

# Retrieve input files from terminal arguments
fasta_file = sys.argv[1]
tree_file = sys.argv[2]

# Set output filename based on fasta filename
tree_outfile = os.path.basename(fasta_file).replace(".fasta", ".treefile")

# Step 1: Extract taxa names (sequence IDs) from fasta file
subsample_taxa = [record.id for record in SeqIO.parse(fasta_file, "fasta")]

# Step 2: Read the full tree from the provided tree file
full_tree = Phylo.read(tree_file, "newick")

# Step 3: Identify taxa that need to be pruned (not present in the subsample)
all_taxa_in_tree = [term.name for term in full_tree.get_terminals()]
prune_taxa = [taxon for taxon in all_taxa_in_tree if taxon not in subsample_taxa]

# Step 4: Prune unwanted taxa from the tree
for taxon in prune_taxa:
    full_tree.prune(taxon)

# Step 5: Save the pruned subtree and remove unnecessary parentheses
Phylo.write(full_tree, tree_outfile, "newick")

# Remove unnecessary parentheses from the tree file
with open(tree_outfile, 'r') as file:
    tree_data = file.read().strip()

cleaned_tree = re.sub(r'^\((.*)\);$', r'\1;', tree_data)

with open(tree_outfile, 'w') as file:
    file.write(cleaned_tree)

print(f"Pruned and cleaned subtree saved as {tree_outfile}")
