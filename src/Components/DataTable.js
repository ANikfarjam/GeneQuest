import React, { useState } from 'react';
import GearButton from "./BLTGearBTN copy";
import './styling/DataTable.css'; // Ensure this path is correct

const DataTable = ({ data, onSelect }) => {
  const [selectedRows, setSelectedRows] = useState([]);

  // Update to store only IDs in selectedRows
  const handleCheckboxChange = (hit) => {
    const isSelected = selectedRows.includes(hit.id);
    const updatedSelection = isSelected
      ? selectedRows.filter((rowId) => rowId !== hit.id)
      : [...selectedRows, hit.id];

    setSelectedRows(updatedSelection);

    if (onSelect) {
      onSelect(updatedSelection);
    }
  };

  // Check if there's any data to display
  if (!data || data.length === 0) {
    return <p className="datatable-no-results">No results to display.</p>;
  }

  // Download selected rows as FASTA
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

  return (
    <div className="datatable-container">
      {/* Pass updated selected data and download logic */}
      <GearButton
        onDownload={downloadSelectedRows}
        selectedData={data.filter((item) => selectedRows.includes(item.id))} // Ensure selectedData is correctly filtered
      />
      <table className="datatable">
        <thead>
          <tr>
            <th>Select</th>
            <th>Num</th>
            <th>ID</th>
            <th>Definition</th>
            <th>Accession</th>
            <th>Length</th>
            <th>Bit Score</th>
            <th>E-value</th>
            <th>Identity</th>
            <th>Hit Sequence</th>
          </tr>
        </thead>
        <tbody>
          {data.map((hit, index) => (
            <tr key={index} className={selectedRows.includes(hit.id) ? 'selected' : ''}>
              <td>
                <input
                  type="checkbox"
                  checked={selectedRows.includes(hit.id)}
                  onChange={() => handleCheckboxChange(hit)}
                />
              </td>
              <td>{hit.num}</td>
              <td>{hit.id}</td>
              <td>{hit.definition}</td>
              <td>{hit.accession}</td>
              <td>{hit.length}</td>
              <td>{hit.bitScore}</td>
              <td>{hit.evalue}</td>
              <td>{hit.identity}</td>
              <td className="datatable-ellipsis">{hit.hitSequence}</td>
            </tr>
          ))}
        </tbody>
      </table>
    </div>
  );
};

export default DataTable;
