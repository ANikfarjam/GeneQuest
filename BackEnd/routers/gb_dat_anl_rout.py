from flask import Blueprint, request, jsonify  # Flask modules for routing and JSON handling
from Bio.Align import MultipleSeqAlignment  # BioPython module for alignments
from Bio.Seq import Seq  # Module for sequence representation
from Bio.SeqRecord import SeqRecord  # Module for representing sequences with metadata

# Create a Flask Blueprint for organizing routes
data_anl_bp = Blueprint('data_anl', __name__)

@data_anl_bp.route("/alignment_conservation/seq_logo", methods=["POST"])
def fetch_seq_logo():
    """
    Endpoint to perform sequence alignment and conservation analysis using Biopython.

    1. Accepts JSON input containing sequence data.
    2. Aligns sequences using MultipleSeqAlignment (manual or pre-aligned data).
    3. Identifies conserved regions in the alignment.
    4. Returns alignment and conserved region information in JSON format.
    """
    try:
        # Step 1: Parse the incoming JSON request
        input_data = request.json  # Expecting a list of {"id", "description", "sequence"}
        if not input_data or not isinstance(input_data, list):
            # Return an error if input data is invalid or not a list
            return jsonify({"error": "Invalid input data"}), 400

        # Step 2: Convert JSON input into Biopython SeqRecord objects
        sequences = [
            SeqRecord(
                Seq(item["sequence"]),  # Sequence string
                id=item["id"],  # Unique identifier for the sequence
                description=item["description"]  # Sequence description
            )
            for item in input_data
        ]

        # Step 3: Create a MultipleSeqAlignment object
        # Note: This example assumes sequences are pre-aligned. If not, you would need
        # an external alignment tool like MAFFT or Clustal Omega to align them first.
        alignment = MultipleSeqAlignment(sequences)

        # Step 4: Identify conserved regions in the alignment
        # A conserved region is a column where all characters are identical
        conserved_regions = []
        for i in range(alignment.get_alignment_length()):
            column = alignment[:, i]  # Extract the i-th column from the alignment
            if len(set(column)) == 1:  # Check if all characters in the column are the same
                conserved_regions.append(i)  # Add the conserved position to the list

        # Step 5: Construct the response data
        response_data = {
            "alignment": [
                {
                    "id": record.id,  # Sequence ID
                    "sequence": str(record.seq),  # Aligned sequence
                    "description": record.description,  # Sequence description
                }
                for record in alignment  # Iterate over each record in the alignment
            ],
            "conserved_regions": conserved_regions,  # List of conserved positions
        }

        # Step 6: Return the alignment and conservation analysis as JSON
        return jsonify(response_data)

    except Exception as e:
        # If any error occurs, log the exception and return an error response
        return jsonify({"error": str(e)}), 500
