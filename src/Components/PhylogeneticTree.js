import React, { useState } from "react";
import axios from "axios";
import Tree from "react-d3-tree";

const PhylogeneticTree = () => {
  const [accessionNumbers, setAccessionNumbers] = useState("");
  const [treeData, setTreeData] = useState(null);
  const [error, setError] = useState(null);
  const [loading, setLoading] = useState(false);

  const handleGenerateTree = async () => {
    setLoading(true);
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
      setLoading(false);
    }
  };

  const transformTreeData = (node) => ({
    name: node.name || "Unnamed",
    attributes: {
      branch_length: node.branch_length ? node.branch_length.toFixed(3) : null, // Show branch length
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
      {loading && (
        <p style={{ color: "blue", marginTop: "10px" }}>
          Generating Tree. Please wait...
        </p>
      )}
      {error && <p style={{ color: "red", marginTop: "10px" }}>{error}</p>}
      {treeData && (
        <div
          style={{
            width: "100%",
            height: "500px", // Fixed height for the tree container
            marginTop: "20px", // Ensure some space below the button
            position: "relative", // Ensure proper positioning
          }}
        >
          <Tree
            data={treeData}
            orientation="horizontal" // Left-to-right rendering
            pathFunc="elbow" // Right-angled branches
            textLayout={{
              textAnchor: "start", // Align labels to the right of nodes
              x: 15, // Position labels slightly to the right
              y: 0, // Keep labels vertically aligned
            }}
            translate={{
              x: 700, // Static horizontal position
              y: 200, // Static vertical position (adjusted to be below the button)
            }}
            zoomable={false} // Disable drag/zoom functionality
            nodeSvgShape={{
              shape: "circle",
              shapeProps: {
                r: 0.5, // Small node size
                fill: "black", // Node color
              },
            }}
            styles={{
              links: {
                stroke: "black", // Branch color
                strokeWidth: 1, // Thin branch lines
              },
              nodes: {
                node: {
                  circle: { fill: "black" }, // Node circle color
                },
                name: {
                  stroke: "none",
                  fill: "black", // Label color
                },
              },
            }}
          />
        </div>
      )}
    </div>
  );
};

export default PhylogeneticTree;
