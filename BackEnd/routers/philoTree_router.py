from flask import Blueprint, request, jsonify
from Bio import Entrez, AlignIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
import tempfile
import subprocess
from dotenv import load_dotenv

# Blueprint setup
phylo_tree_bp = Blueprint('phylo_tree_generator', __name__)

load_dotenv()

# Set Entrez email
Entrez.email = os.getenv("NCBI_EMAIL")


@phylo_tree_bp.route('/phylo_tree', methods=['POST'])
def phylogenetic_tree():
    try:
        # Get accession numbers from request
        data = request.get_json()
        accession_numbers = data.get('accession_numbers', [])
        
        if not accession_numbers:
            return jsonify({"error": "No accession numbers provided"}), 400
        
        # Fetch sequences from GenBank
        sequences = []
        for acc in accession_numbers:
            with Entrez.efetch(db="nucleotide", id=acc, rettype="fasta", retmode="text") as handle:
                fasta_data = handle.read()
                seq_record = SeqRecord(Seq("".join(fasta_data.split('\n')[1:])), id=acc)
                sequences.append(seq_record)
        
        # Write sequences to a temporary file
        with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.fasta') as fasta_file:
            for seq in sequences:
                fasta_file.write(f">{seq.id}\n{seq.seq}\n")
            input_file = fasta_file.name

        # Perform sequence alignment using MAFFT
        aligned_file = tempfile.NamedTemporaryFile(delete=False, suffix='.fasta').name
        command = [
            "mafft.bat",  # MAFFT batch file for Windows
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

        return jsonify({"phylogenetic_tree": tree_json}), 200

    except Exception as e:
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


