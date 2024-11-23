from flask import Blueprint, request, jsonify
from Bio.Blast import NCBIWWW, NCBIXML
import os
import json
import hashlib
from dotenv import load_dotenv

# Blueprint setup
megablast_bp = Blueprint("megablast", __name__)

load_dotenv()

# Define the cache file
CACHE_FILE = 'megablast_cache.json'

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


@megablast_bp.route('/api/megablast', methods=['POST'])
def megablast():
    try:
        # Get the sequence from the request
        data = request.get_json()
        sequence = data.get("sequence", None)

        if not sequence:
            return jsonify({"error": "No sequence provided"}), 400

        print(f"Received sequence: {sequence}")  # Debugging

        # Generate a unique cache key for the sequence
        sequence_key = hashlib.md5(sequence.encode()).hexdigest()

        # Check if the result is already in the cache
        if sequence_key in cache:
            print("Returning cached MegaBLAST results")
            return jsonify({"megablast_results": cache[sequence_key]}), 200

        # Perform MegaBLAST
        result_handle = NCBIWWW.qblast(
            program="blastn",
            database="nt",
            sequence=sequence,
            megablast=True
        )

        # Parse the BLAST result
        blast_records = NCBIXML.parse(result_handle)
        blast_results = []
        for record in blast_records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    blast_results.append({
                        "title": alignment.title,
                        "length": alignment.length,
                        "score": hsp.score,
                        "e_value": hsp.expect,
                        "identity": hsp.identities,
                        "query": hsp.query,
                        "match": hsp.match,
                        "subject": hsp.sbjct,
                    })

        # Cache the MegaBLAST results
        cache[sequence_key] = blast_results
        save_cache(cache)

        print(f"BLAST Results: {blast_results}")  # Debugging

        return jsonify({"megablast_results": blast_results}), 200

    except Exception as e:
        print(f"Error: {e}")  # Debugging
        return jsonify({"error": str(e)}), 500
