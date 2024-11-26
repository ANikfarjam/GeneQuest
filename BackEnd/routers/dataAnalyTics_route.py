from flask import Blueprint, request, jsonify, send_file
from .ConservedRegionCalc import calculate_conserved_region
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
anlytics_bp = Blueprint("analytics_bp", __name__)
@anlytics_bp.route('/conservativeRegion', methods=["POST"])
def get_conserved():
    try:
        # Parse the input JSON
        data = request.json
        print(data)
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
        print(seq_records)
        #get conservative Records
        conserved_region = calculate_conserved_region(seq_records)
        return jsonify({"Conserved region": conserved_region})
    except Exception as e:
        print(f"Error: {e}")
        return jsonify({"error": str(e)}), 500
