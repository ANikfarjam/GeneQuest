from flask import Blueprint, request, jsonify, send_file
from Bio import Phylo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MafftCommandline
from Bio.Align import MultipleSeqAlignment
from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import tempfile
import os
import matplotlib
import matplotlib.pyplot as plt

# Use the Agg backend for non-GUI rendering
matplotlib.use('Agg')

data_anl_bp = Blueprint('data_anl', __name__)

def align_sequences_mafft(sequences):
    """
    Align sequences using MAFFT.
    """
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as input_file, \
            tempfile.NamedTemporaryFile(mode='r', delete=False) as output_file:

        # Write sequences to input FASTA file
        SeqIO.write(sequences, input_file.name, "fasta")
        input_file.close()

        # Run MAFFT alignment
        mafft_cline = MafftCommandline(input=input_file.name)
        mafft_cline.set_parameter('--auto', True)  # Let MAFFT choose the best strategy
        mafft_cline.set_parameter('--quiet', True)  # Suppress console output

        stdout, stderr = mafft_cline()
        with open(output_file.name, "w") as f:
            f.write(stdout)

        # Parse the aligned sequences
        alignment = SeqIO.parse(output_file.name, "fasta")
        return list(alignment)

def generate_tree(sequences):
    """
    Generate a phylogenetic tree from aligned sequences.
    """
    # Align sequences using MAFFT
    aligned_sequences = align_sequences_mafft(sequences)

    # Convert aligned sequences to Biopython's MultipleSeqAlignment
    alignment = MultipleSeqAlignment(aligned_sequences)

    # Calculate distance matrix
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)

    # Construct a phylogenetic tree using UPGMA algorithm
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(distance_matrix)

    return tree

def save_tree_to_image(tree):
    """
    Save the phylogenetic tree to a temporary PNG image file.
    """
    with tempfile.NamedTemporaryFile(delete=False, suffix=".png") as img_file:
        plt.figure(figsize=(10, 8))  # Adjust figure size as needed
        Phylo.draw(tree, do_show=False)  # Use Biopython's Phylo.draw function
        plt.savefig(img_file.name, format="png")  # Save the figure to the file
        plt.close()  # Close the figure to free memory
        return img_file.name

@data_anl_bp.route("/gb_dat_anl", methods=["POST"])
def gb_dat_anl():
    try:
        # Parse the input JSON
        data = request.json
        if not data or "sequences" not in data:
            return jsonify({"error": "No sequences provided"}), 400

        sequences = data["sequences"]

        # Validate input sequences
        if not all("hitSequence" in seq and seq["hitSequence"] for seq in sequences):
            return jsonify({"error": "All sequences must have a non-empty 'hitSequence' field."}), 400

        # Extract `requestedQuerySequence` if available
        requested_query_sequence = sequences[0].get("requestedQuerySequence", None)
        if requested_query_sequence:
            sequences.append({"accession": "requestedQuery", "hitSequence": requested_query_sequence})

        # Convert to Biopython SeqRecord objects
        seq_records = [SeqRecord(Seq(seq["hitSequence"]), id=seq["accession"]) for seq in sequences]

        # Generate phylogenetic tree
        tree = generate_tree(seq_records)

        # Save tree to an image file
        img_path = save_tree_to_image(tree)

        # Return the image as a response
        return send_file(img_path, mimetype="image/png")

    except Exception as e:
        print(f"Error during tree generation: {e}")
        return jsonify({"error": str(e)}), 500
    finally:
        # Clean up temporary files if needed
        if "img_path" in locals() and os.path.exists(img_path):
            os.unlink(img_path)
