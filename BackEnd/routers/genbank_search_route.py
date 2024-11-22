from dotenv import load_dotenv
import os
from flask import Flask, Blueprint, request, jsonify
#from flask_cors import CORS
from Bio import Entrez
import ssl
import json

# Load environment variables
load_dotenv()

# Disable SSL verification (for development purposes)
ssl._create_default_https_context = ssl._create_unverified_context

# Blueprint setup
gen_bank_bp = Blueprint('gen_bank_search', __name__)

# Set Entrez email
Entrez.email = os.getenv("NCBI_EMAIL")

# Define the cache file
CACHE_FILE = 'query_cache.json'

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

@gen_bank_bp.route('/api/genbank_search', methods=['POST'])
def genbank_search():
    try:
        data = request.get_json()
        query = data.get('query', None)

        if not query:
            return jsonify({"error": "Query parameter is required"}), 400
        print(f"Received GenBank query: {query}")
        print(f"Received query: {query}")  # Log the query

        # Check cache for the query
        if query in cache:
            print("Returning cached result")
            return jsonify(cache[query]), 200

        # Search for IDs matching the query
        with Entrez.esearch(db="nucleotide", term=query, retmax=10) as handle:
            search_results = Entrez.read(handle)
            id_list = search_results.get("IdList", [])
            print(f"Search results: {search_results}")  # Log the search results
    

        id_list = search_results.get("IdList", [])

        # Fetch sequence data for the IDs
        if id_list:
            with Entrez.efetch(db="nucleotide", id=",".join(id_list), rettype="fasta", retmode="text") as fetch_handle:
                sequences = fetch_handle.read() #result of our search

            print(f"Fetched sequences:\n{sequences[:500]}")  # Log a snippet of the sequences

            # Save result to cache
            result_data = {
                "query": query,
                "sequences": sequences
            }

            cache[query] = result_data
            save_cache(cache)

            return jsonify(result_data), 200
        else:
            return jsonify({"error": "No results found for the query"}), 404

    except Exception as e:
        print(f"Internal Server Error: {e}")
        return jsonify({"error": "Internal server ermror"}), 500


