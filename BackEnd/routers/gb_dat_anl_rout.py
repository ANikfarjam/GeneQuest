from flask import Flask, Blueprint, request, jsonify
from Bio import AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Align import MultipleSeqAlignment


#Blueprint setup
data_anl_bp = Blueprint('data_anl', __name__)
#set up routes
@data_anl_bp.route("/alignment_conservation/seq_logo", methods=["POST"])
def fetch_seq_logo():
    try:
        # Load input JSON
        input_data = request.json
        input_file = "input_sequences.fasta"
        output_file = "aligned_output.aln"

        # Write input sequences to FASTA
        with open(input_file, "w") as fasta_file:
            for item in input_data:
                fasta_file.write(f">{item['id']} {item['description']}\n{item['sequence']}\n\n")

        # Run Clustal Omega
        clustal_cmd = ClustalOmegaCommandline(infile=input_file, outfile=output_file, verbose=True, auto=True)
        clustal_cmd()

        # Load alignment
        alignment = AlignIO.read(output_file, "clustal")

        # Identify conserved regions
        conserved_regions = []
        for i in range(alignment.get_alignment_length()):
            column = alignment[:, i]
            if len(set(column)) == 1:
                conserved_regions.append(i)

        # Clean up
        os.remove(input_file)
        os.remove(output_file)

        # Return JSON response
        return jsonify({"conserved_regions": conserved_regions, "alignment": str(alignment)})
    except Exception as e:
        return jsonify({"error": str(e)}), 500
@data_anl_bp.route('/alignment_conservation/mut_anl/var_map', methods=['GET'])
def var_map():
    pass
@data_anl_bp.route('/alignment_conservation/mut_anl/pi_chart', methods=['GET'])
def pie_chart():
    pass



