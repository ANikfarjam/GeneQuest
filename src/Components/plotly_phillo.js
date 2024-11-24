import React from "react";
import Plot from "react-plotly.js";

const PhyloTree = ({ Data }) => {
  if (!Data) return <p>No data available for tree generation.</p>;

  // Helper function to recursively process the tree
  const processTree = (node, x = 0, y = 0, stepY = 1) => {
    const result = [];
    const edges = [];

    if (node.name) {
      result.push({ x, y, name: node.name }); // Add the current node
    }

    if (node.children) {
      const totalChildren = node.children.length;
      node.children.forEach((child, i) => {
        const childY = y - stepY * (totalChildren - 1) / 2 + i * stepY; // Spread children vertically
        const childX = x + 1; // Move children to the right
        result.push(...processTree(child, childX, childY, stepY / 2)); // Recurse
        edges.push({ x: [x, childX], y: [y, childY] }); // Add edge from parent to child
      });
    }

    return { nodes: result, edges };
  };

  const { nodes, edges } = processTree(Data);

  return (
    <Plot
      data={[
        {
          x: nodes.map((n) => n.x),
          y: nodes.map((n) => n.y),
          text: nodes.map((n) => n.name),
          mode: "markers+text",
          type: "scatter",
          marker: { size: 10 },
          textposition: "top center",
        },
        ...edges.map((edge) => ({
          x: edge.x,
          y: edge.y,
          mode: "lines",
          type: "scatter",
          line: { color: "gray", width: 1 },
        })),
      ]}
      layout={{
        title: "Phylogenetic Tree",
        xaxis: { showgrid: false, zeroline: false, visible: false },
        yaxis: { showgrid: false, zeroline: false, visible: false },
        margin: { l: 0, r: 0, t: 40, b: 0 },
        height: 600,
      }}
    />
  );
};

export default PhyloTree;
