import React, { useState } from "react";
import axios from "axios";
import Tree from "react-d3-tree";

const PhylogeneticTree = () => {
  const [accessionNumbers, setAccessionNumbers] = useState("");
  const [treeData, setTreeData] = useState(null);
  const [error, setError] = useState(null);
  const [loading, setLoading] = useState(false); // New state for loading

  const handleGenerateTree = async () => {
    setLoading(true); // Set loading to true
    setError(null);
    setTreeData(null);

    try {
      const response = await axios.post("http://127.0.0.1:5000/api/phylo_tree", {
        accession_numbers: accessionNumbers.split(",").map((acc) => acc.trim()),
      });
      setTreeData(transformTreeData(response.data.phylogenetic_tree));
    } catch (err) {
      console.error("Error generating phylogenetic tree:", err);
      setError(err.response?.data?.error || "An unknown error occurred.");
    } finally {
      setLoading(false); // Set loading to false
    }
  };

  const transformTreeData = (node) => ({
    name: node.name || "Unnamed",
    attributes: {
      branch_length: node.branch_length ? node.branch_length.toFixed(3) : null,
    },
    children: node.children ? node.children.map(transformTreeData) : [],
  });

  return (
    <div>
      <h2>Generate Phylogenetic Tree</h2>
      <textarea
        rows="5"
        placeholder="Enter accession numbers separated by commas"
        value={accessionNumbers}
        onChange={(e) => setAccessionNumbers(e.target.value)}
      />
      <br />
      <button onClick={handleGenerateTree}>Generate Tree</button>
      {loading && <p style={{ color: "blue", marginTop: "10px" }}>Generating Tree. Please wait...</p>} {/* Loading Message */}
      {error && <p style={{ color: "red", marginTop: "10px" }}>{error}</p>}
      {treeData && (
        <div style={{ width: "100%", height: "500px", border: "1px solid black", marginTop: "20px" }}>
          <Tree data={treeData} orientation="vertical"
          textLayout={{
            textAnchor: "middle", // Centers the label horizontally
            x: 0, // Keeps the label aligned with the node
            y: 8, // Positions the label below the node
          }}
          />
        </div>
      )}
    </div>
  );
};

export default PhylogeneticTree;

