from flask import Blueprint, request, jsonify  # Flask modules for routing and JSON handling
from Bio.Align import MultipleSeqAlignment  # BioPython module for alignments
from Bio.Seq import Seq  # Module for sequence representation
from Bio.SeqRecord import SeqRecord  # Module for representing sequences with metadata
from Bio import SeqIO  # For sequence file handling
from Bio.Align.Applications import ClustalOmegaCommandline  # Alignment tool
import tempfile  # For creating temporary files

# Create a Flask Blueprint for organizing routes
data_anl_bp = Blueprint('data_anl', __name__)

def align_sequences(sequences):
    """
    Align sequences using Clustal Omega and return aligned sequences.
    """
    try:
        with tempfile.NamedTemporaryFile(delete=False, mode='w') as input_file, \
                tempfile.NamedTemporaryFile(delete=False) as output_file:
            # Write sequences to a FASTA file
            SeqIO.write(sequences, input_file.name, "fasta")
            input_file.close()

            # Run Clustal Omega
            clustal_cline = ClustalOmegaCommandline(
                infile=input_file.name,
                outfile=output_file.name,
                verbose=True,
                auto=True
            )
            clustal_cline()

            # Parse and return aligned sequences
            return list(SeqIO.parse(output_file.name, "fasta"))
    except Exception as e:
        raise RuntimeError(f"Sequence alignment failed: {e}")

@data_anl_bp.route("/alignment_conservation/seq_logo", methods=["POST"])
def fetch_seq_logo():
    """
    Endpoint to perform sequence alignment and conservation analysis using Biopython.
    """
    try:
        # Step 1: Parse the incoming JSON request
        input_data = request.json  # Expecting a list of {"id", "description", "sequence"}
        if not input_data or not isinstance(input_data, list):
            return jsonify({"error": "Invalid input data. Must be a list of sequences."}), 400

        # Validate sequence data structure
        if not all("id" in item and "sequence" in item for item in input_data):
            return jsonify({"error": "Each sequence must have 'id' and 'sequence' fields."}), 400

        # Step 2: Convert JSON input into Biopython SeqRecord objects
        sequences = [
            SeqRecord(
                Seq(item["sequence"]),
                id=item["id"],
                description=item.get("description", "")
            )
            for item in input_data
        ]

        # Step 3: Align sequences if they are not pre-aligned
        try:
            aligned_sequences = align_sequences(sequences)
        except RuntimeError as align_error:
            return jsonify({"error": str(align_error)}), 500

        # Create a MultipleSeqAlignment object
        alignment = MultipleSeqAlignment(aligned_sequences)

        # Step 4: Identify conserved regions in the alignment
        conserved_regions = [
            i for i in range(alignment.get_alignment_length())
            if len(set(alignment[:, i])) == 1
        ]

        # Step 5: Construct the response data
        response_data = {
            "alignment": [
                {
                    "id": record.id,
                    "sequence": str(record.seq),
                    "description": record.description,
                }
                for record in alignment
            ],
            "conserved_regions": conserved_regions,
        }

        # Step 6: Return the alignment and conservation analysis as JSON
        return jsonify(response_data), 200

    except Exception as e:
        # If any error occurs, log the exception and return an error response
        return jsonify({"error": f"Server error: {str(e)}"}), 500
