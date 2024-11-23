from flask import Flask, jsonify
from flask_cors import CORS
from routers.genbank_search_route import gen_bank_bp
from routers.philoTree_router import phylo_tree_bp
from routers.gb_dat_anl_rout  import data_anl_bp
from routers.megablast_route import megablast_bp
from dotenv import load_dotenv
import os

# Load environment variables
load_dotenv()

from routers.gb_dat_anl_rout import data_anl_bp
app = Flask(__name__)
# allow universal requests
CORS(app, resources={
    r"/api/*": {"origins": "http://localhost:3000"},
    r"/alignment_conservation/*": {"origins": "http://localhost:3000"}
})

#register genbank search
app.register_blueprint(gen_bank_bp) #gen-bank search
app.register_blueprint(phylo_tree_bp) #fetch phylo tree info
app.register_blueprint(data_anl_bp) #fetch genbank data analysis
app.register_blueprint(megablast_bp) #megablast route
if __name__ == '__main__':
    app.run(debug=True)