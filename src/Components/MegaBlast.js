import React, { useState } from "react";
import axios from "axios";
import "./DataTable";

const MegaBLAST = () => {
  const [sequence, setSequence] = useState("");
  const [results, setResults] = useState([]);
  const [error, setError] = useState(null);
  const [loading, setLoading] = useState(false);

  // Initial column widths
  const [columnWidths] = useState({
    title: "50%", // Wider Title column
    length: "12.5%",
    score: "12.5%",
    eValue: "12.5%",
    identity: "12.5%",
  });

  const handleMegaBLAST = async () => {
    setLoading(true);
    setError(null);
    setResults([]);

    try {
      const response = await axios.post("http://127.0.0.1:5000/api/megablast", {
        sequence: sequence.trim(),
      });
      setResults(response.data.megablast_results);
    } catch (err) {
      console.error("Error running MegaBLAST:", err);
      setError(err.response?.data?.error || "An unknown error occurred.");
    } finally {
      setLoading(false);
    }
  };

  return (
    <div>
      <h2>MegaBLAST Tool</h2>
      <textarea
        rows="5"
        placeholder="Enter your sequence here"
        value={sequence}
        onChange={(e) => setSequence(e.target.value)}
      />
      <br />
      <button onClick={handleMegaBLAST}>Run MegaBLAST</button>

      {loading && <p style={{ color: "blue", marginTop: "10px" }}>Running MegaBLAST. Please wait...</p>}
      {error && <p style={{ color: "red", marginTop: "10px" }}>{error}</p>}

      {results.length > 0 && (
        <div className="datatable-container">
          <table className="datatable">
            <thead>
              <tr>
                <th style={{ width: columnWidths.title }}>Title</th>
                <th style={{ width: columnWidths.length }}>Length</th>
                <th style={{ width: columnWidths.score }}>Score</th>
                <th style={{ width: columnWidths.eValue }}>E-Value</th>
                <th style={{ width: columnWidths.identity }}>Identity</th>
              </tr>
            </thead>
            <tbody>
              {results.map((result, index) => (
                <tr key={index}>
                  <td>{result.title}</td>
                  <td>{result.length}</td>
                  <td>{result.score}</td>
                  <td>{result.e_value}</td>
                  <td>{result.identity}</td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </div>
  );
};

export default MegaBLAST;