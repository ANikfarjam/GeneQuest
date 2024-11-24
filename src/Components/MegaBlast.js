import React, { useState } from "react";
import axios from "axios";
import "./DataTable";

const MegaBLAST = () => {
  const [sequence, setSequence] = useState("");
  const [file, setFile] = useState(null); // State to store uploaded file
  const [results, setResults] = useState([]);
  const [error, setError] = useState(null);
  const [loading, setLoading] = useState(false);

  const handleFileChange = (e) => {
    setFile(e.target.files[0]);
    setSequence(""); // Clear the sequence input if a file is selected
  };

  const handleMegaBLAST = async () => {
    setLoading(true);
    setError(null);
    setResults([]);

    const formData = new FormData();

    if (file) {
      formData.append("file", file); // Add file to form data
    } else if (sequence.trim()) {
      formData.append("sequence", sequence.trim()); // Add sequence to form data
    } else {
      setError("Please provide a sequence or upload a FASTA file.");
      setLoading(false);
      return;
    }

    try {
      const response = await axios.post("http://127.0.0.1:5000/api/megablast", formData, {
        headers: {
          "Content-Type": "multipart/form-data",
        },
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

      {/* Sequence Textarea */}
      <textarea
        rows="5"
        placeholder="Enter your sequence here"
        value={sequence}
        onChange={(e) => setSequence(e.target.value)}
        disabled={!!file} // Disable if a file is uploaded
      />
      <br />

      {/* Hidden File Input */}
      <input
        id="fileInput"
        type="file"
        accept=".fasta"
        onChange={handleFileChange}
        style={{ display: "none" }} // Hide the default file input
      />

      {/* Custom File Button */}
      <button
        onClick={() => document.getElementById("fileInput").click()} // Trigger the file input click
      >
        {file ? file.name : "Choose File"}
      </button>
      <br />

      <button onClick={handleMegaBLAST}>Run MegaBLAST</button>

      {/* Loading Message */}
      {loading && <p style={{ color: "blue", marginTop: "10px" }}>Running MegaBLAST. Please wait...</p>}

      {/* Error Message */}
      {error && <p style={{ color: "red", marginTop: "10px" }}>{error}</p>}

      {/* Results Table */}
      {results.length > 0 && (
        <div className="datatable-container">
          <table className="datatable">
            <thead>
              <tr>
                <th>Title</th>
                <th>Accession</th>
                <th>Length</th>
                <th>Score</th>
                <th>E-Value</th>
                <th>Identity</th>
                <th>Hit Sequence</th>
              </tr>
            </thead>
            <tbody>
              {results.map((result, index) => (
                <tr key={index}>
                  <td>{result.title}</td>
                  <td>{result.accession}</td>
                  <td>{result.length}</td>
                  <td>{result.score}</td>
                  <td>{result.e_value}</td>
                  <td>{result.identity}</td>
                  <td>{result.hit_sequence}</td>
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


