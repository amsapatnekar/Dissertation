import os
from Bio import SeqIO
import yaml

# Input file
input_file = "filtered_psychrophile_proteins.fasta"

# Output folder
output_folder = "Psychrophile_yamls"

# Minimum sequence length to exclude
min_length = 1500

# Create the output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Function to process sequences and generate YAML files
def generate_yaml_files(input_file, output_folder, min_length):
    # Parse the input FASTA file
    fasta_sequences = SeqIO.parse(input_file, "fasta")
    
    count = 0

    for record in fasta_sequences:
        # Process only sequences shorter than or equal to the min_length
        if len(record.seq) <= min_length:
            yaml_data = {
                "version": 1,
                "sequences": [
                    {
                        "protein": {
                            "id": "A",
                            "sequence": str(record.seq)
                        }
                    }
                ]
            }

            # Generate YAML file name
            yaml_filename = os.path.join(output_folder, f"{record.id}.yaml")

            # Write the YAML file
            with open(yaml_filename, "w") as yaml_file:
                yaml.dump(yaml_data, yaml_file, default_flow_style=False, sort_keys=False)
            
            count += 1
    
    print(f"Generated {count} YAML files in the folder '{output_folder}'.")

# Run the function
generate_yaml_files(input_file, output_folder, min_length)
