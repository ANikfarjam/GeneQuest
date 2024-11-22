import React, { useState } from 'react';
import './styling/DataTable.css'; // Make sure this path is correct

const DataTable = ({ data, onSelect }) => {
  const [selectedRows, setSelectedRows] = useState([]);

  const handleCheckboxChange = (hit) => {
    const isSelected = selectedRows.includes(hit);
    const updatedSelection = isSelected
      ? selectedRows.filter((row) => row !== hit)
      : [...selectedRows, hit];

    setSelectedRows(updatedSelection);

    if (onSelect) {
      onSelect(updatedSelection);
    }
  };

  if (!data || data.length === 0) {
    return <p className="datatable-no-results">No results to display.</p>;
  }

  return (
    <div className="datatable-container">
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
            <tr key={index} className={selectedRows.includes(hit) ? 'selected' : ''}>
              <td>
                <input
                  type="checkbox"
                  checked={selectedRows.includes(hit)}
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
