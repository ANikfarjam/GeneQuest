from flask import Blueprint, request, jsonify
from Bio import Entrez, AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
import tempfile
import subprocess
from dotenv import load_dotenv
import json
import hashlib

# Blueprint setup
phylo_tree_bp = Blueprint('phylo_tree_generator', __name__)
load_dotenv()

# Set Entrez email
Entrez.email = os.getenv("NCBI_EMAIL")

# Define the cache file
CACHE_FILE = 'phylo_tree_cache.json'

# Load the cache from the file if it exists
def load_cache():
    if os.path.exists(CACHE_FILE):
        with open(CACHE_FILE, 'r') as cache_file:
            return json.load(cache_file)
    else:
        return {}

# Save the cache to the file
def save_cache(cache_data):
    with open(CACHE_FILE, 'w') as cache_file:
        json.dump(cache_data, cache_file, indent=4)

# Initialize the cache (loaded once when the server starts)
cache = load_cache()

@phylo_tree_bp.route('/api/phylo_tree', methods=['POST'])
def phylogenetic_tree():
    try:
        # Get accession numbers from request
        data = request.get_json()
        accession_numbers = data.get('accession_numbers', [])

        if not accession_numbers:
            return jsonify({"error": "No accession numbers provided"}), 400

        # Generate a unique key for the accession numbers
        accession_key = hashlib.md5(json.dumps(sorted(accession_numbers)).encode()).hexdigest()

        # Check if the result is already in the cache
        if accession_key in cache:
            print("Returning cached phylogenetic tree")
            return jsonify({"phylogenetic_tree": cache[accession_key]}), 200

        # Fetch sequences from GenBank
        sequences = []
        for acc in accession_numbers:
            if acc in cache:  # Check if the sequence is already cached
                print(f"Using cached sequence for {acc}")
                sequences.append(SeqRecord(Seq(cache[acc]), id=acc))
            else:
                with Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text") as handle:
                    fasta_data = handle.read()
                    sequence = "".join(fasta_data.split('\n')[1:])
                    sequences.append(SeqRecord(Seq(sequence), id=acc))
                    cache[acc] = sequence  # Cache the fetched sequence

        # Save updated cache for sequences
        save_cache(cache)

        # Write sequences to a temporary file
        with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fasta') as fasta_file:
            for seq in sequences:
                fasta_file.write(f">{seq.id}\n{seq.seq}\n")
            input_file = fasta_file.name

        # Perform sequence alignment using MAFFT
        aligned_file = tempfile.NamedTemporaryFile(delete=False, suffix='.fasta').name
        command = [
            "../mafft-mac/mafft.bat",  # MAFFT batch file for Windows
            "--auto",  # Use automatic alignment settings
            input_file  # Input file
        ]

        with open(aligned_file, "w") as output:
            subprocess.run(command, stdout=output, check=True)

        # Read aligned sequences
        alignment = AlignIO.read(aligned_file, "fasta")

        # Calculate distance matrix
        calculator = DistanceCalculator("identity")
        distance_matrix = calculator.get_distance(alignment)

        # Build phylogenetic tree using UPGMA
        constructor = DistanceTreeConstructor(calculator, method="upgma")
        tree = constructor.build_tree(alignment)

        # Convert tree to JSON
        tree_json = tree_to_json(tree)

        # Cache the phylogenetic tree
        cache[accession_key] = tree_json
        save_cache(cache)

        return jsonify({"phylogenetic_tree": tree_json}), 200

    except Exception as e:
        print(f"Error: {e}")
        return jsonify({"error": str(e)}), 500


def tree_to_json(tree):
    """
    Convert a Phylo tree to a JSON structure.
    """
    def node_to_dict(node):
        children = [node_to_dict(child) for child in node.clades] if node.clades else []
        return {
            "name": str(node.name) if node.name else "Unnamed",
            "branch_length": node.branch_length,
            "children": children
        }
    
    return node_to_dict(tree.root)
