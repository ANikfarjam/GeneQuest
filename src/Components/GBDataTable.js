import React, { useState } from "react";
import GearButton from "./GearBTN";
import axios from "axios";
import BarChartConserved from "./BarChartConserved"; // Import Bar chart
import SequenceLogo from "react-sequence-logos"; // Import Sequence Logo library
import Modal from "@mui/material/Modal";
import Box from "@mui/material/Box";
import "./styling/DataTable.css";

const GBDataTable = ({ data }) => {
  const [selectedRows, setSelectedRows] = useState([]);
  const [conservedData, setConservedData] = useState(null); // Store conservation results
  const [sequenceLogoData, setSequenceLogoData] = useState(null); // Store sequence logo data
  const [showModal, setShowModal] = useState(false); // Modal visibility state

  const handleCheckboxChange = (id) => {
    const isSelected = selectedRows.includes(id);
    const updatedSelection = isSelected
      ? selectedRows.filter((row) => row !== id)
      : [...selectedRows, id];
    setSelectedRows(updatedSelection);
  };

  const downloadSelectedRows = () => {
    if (selectedRows.length === 0) {
      alert("No rows selected to download!");
      return;
    }

    const selectedData = data.filter((item) => selectedRows.includes(item.id));
    const fastaContent = selectedData
      .map((item) => `>${item.id} ${item.description}\n${item.sequence}`)
      .join("\n\n");

    const blob = new Blob([fastaContent], { type: "text/plain" });
    const link = document.createElement("a");
    link.href = URL.createObjectURL(blob);
    link.download = "selected_sequences.fasta";
    link.click();
  };

  const sendSelectedDataForAnalysis = async () => {
    if (selectedRows.length === 0) {
      alert("No rows selected for analysis!");
      return;
    }

    const selectedData = data.filter((item) => selectedRows.includes(item.id));
    try {
      const response = await axios.post(
        "http://127.0.0.1:5000/alignment_conservation/seq_logo",
        selectedData
      );

      // Process the response data for visualization
      setConservedData(response.data.conserved_regions);
      setSequenceLogoData(response.data.alignment); // Alignment for sequence logo
      setShowModal(true); // Show modal with charts
    } catch (error) {
      console.error("Error during analysis:", error);
    }
  };

  const closeModal = () => {
    setShowModal(false);
  };

  if (!data || data.length === 0) {
    return <p className="datatable-no-results">No results to display.</p>;
  }

  return (
    <div className="datatable-container">
      {/* Pass download and data analytics functions to GearButton */}
      <GearButton
        onDownload={downloadSelectedRows}
        onDataAnalytics={sendSelectedDataForAnalysis}
      />
      <table className="datatable">
        <thead>
          <tr>
            <th>Select</th>
            <th>ID</th>
            <th>Description</th>
            <th>Sequence</th>
          </tr>
        </thead>
        <tbody>
          {data.map((item, index) => (
            <tr
              key={index}
              className={selectedRows.includes(item.id) ? "selected" : ""}
            >
              <td>
                <input
                  type="checkbox"
                  checked={selectedRows.includes(item.id)}
                  onChange={() => handleCheckboxChange(item.id)}
                />
              </td>
              <td>{item.id}</td>
              <td>{item.description}</td>
              <td className="datatable-ellipsis">{item.sequence}</td>
            </tr>
          ))}
        </tbody>
      </table>

      {/* Modal for Charts */}
      <Modal open={showModal} onClose={closeModal}>
        <Box
          sx={{
            position: "absolute",
            top: "50%",
            left: "50%",
            transform: "translate(-50%, -50%)",
            width: 600,
            bgcolor: "background.paper",
            boxShadow: 24,
            p: 4,
            borderRadius: 2,
          }}
        >
          <h2>Conservation Analysis Results</h2>
          {conservedData && (
            <BarChartConserved data={conservedData} />
          )}
          {sequenceLogoData && (
            <div>
              <h3>Sequence Logo</h3>
              <SequenceLogo data={sequenceLogoData} />
            </div>
          )}
        </Box>
      </Modal>
    </div>
  );
};

export default GBDataTable;
