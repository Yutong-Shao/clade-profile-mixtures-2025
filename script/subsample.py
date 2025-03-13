import sys
import random
import os

# Retrieve input arguments from command line
input_fasta = sys.argv[1]
second_argument = sys.argv[2]

# Extract base filename without extension
base_name = os.path.splitext(os.path.basename(input_fasta))[0]

# Read sequences from fasta file
with open(input_fasta, 'r') as file:
    lines = file.readlines()

# Pair taxa names with their sequences
sequences = []
for i in range(0, len(lines), 2):
    if lines[i].startswith('>'):
        sequences.append((lines[i].strip(), lines[i + 1].strip()))

# Check if second argument is numeric (random selection) or a taxa file (selection from list)
if second_argument.isdigit():
    num_taxa_to_select = int(second_argument)
    # Check if requested number of taxa is valid
    if num_taxa_to_select > len(sequences):
        print(f"Error: Requested number of taxa ({num_taxa_to_select}) exceeds available taxa ({len(sequences)}).")
        sys.exit(1)

    # Randomly select specified number of taxa
    selected_sequences = random.sample(sequences, num_taxa_to_select)

    # Output filename for random selection
    output_fasta = f"{base_name}_selected_{num_taxa_to_select}.fasta"
else:
    taxa_file = second_argument
    if not os.path.exists(taxa_file):
        print(f"Error: Taxa file '{taxa_file}' not found.")
        sys.exit(1)

    # Read taxa names from taxa file
    with open(taxa_file, 'r') as tf:
        taxa_list = [line.strip() for line in tf.readlines()]

    # Select sequences matching taxa from taxa file
    selected_sequences = [seq for seq in sequences if seq[0][1:] in taxa_list]

    # Check if all taxa were found
    found_taxa = [seq[0][1:] for seq in selected_sequences]
    missing_taxa = set(taxa_list) - set(found_taxa)
    if missing_taxa:
        print(f"Warning: These taxa were not found in fasta file: {', '.join(missing_taxa)}")

    # Output filename for taxa list selection
    output_fasta = f"{base_name}_selected_taxa.fasta"

# Write selected sequences to output fasta file
with open(output_fasta, 'w') as output_handle:
    for name, seq in selected_sequences:
        output_handle.write(f"{name}\n{seq}\n")

print(f"Successfully selected taxa, saved to {output_fasta}.")
print(f"Total number of taxa selected: {len(selected_sequences)}")
