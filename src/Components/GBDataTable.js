import React, { useState } from "react";
import GearButton from "./GearBTN";
import "./styling/DataTable.css";

const GBDataTable = ({ data }) => {
  const [selectedRows, setSelectedRows] = useState([]);

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

  return (
    <div className="datatable-container">
      <GearButton
        onDownload={downloadSelectedRows}
        selectedData={data.filter((item) => selectedRows.includes(item.id))} // Pass selected data
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
    </div>
  );
};

export default GBDataTable;
