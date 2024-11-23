import React, { useState } from "react";
import { Menu, MenuItem, IconButton, Modal, Box } from "@mui/material";
import SettingsIcon from "@mui/icons-material/Settings";

const GearButton = ({ onDownload, selectedData }) => {
  const [anchorEl, setAnchorEl] = useState(null);
  const [showModal, setShowModal] = useState(false); // Modal visibility state
  const [alignmentData, setAlignmentData] = useState(null); // Alignment data for the modal

  const open = Boolean(anchorEl);

  const handleClick = (event) => {
    setAnchorEl(event.currentTarget);
  };

  const handleClose = () => {
    setAnchorEl(null);
  };

  const handleShowAlignment = (data) => {
    if (!data || data.length === 0) {
      alert("No data selected for alignment!");
      return;
    }

    const alignment = data.map((item) => ({
      querySequence: item.querySequence,
      hitSequence: item.hitSequence,
      midline: item.midline,
    }));

    setAlignmentData(alignment);
    setShowModal(true);
  };

  const closeModal = () => {
    setShowModal(false);
    setAlignmentData(null);
  };

  return (
    <div>
      <IconButton
        aria-label="settings"
        aria-controls="gear-menu"
        aria-haspopup="true"
        onClick={handleClick}
      >
        <SettingsIcon />
      </IconButton>

      <Menu
        id="gear-menu"
        anchorEl={anchorEl}
        open={open}
        onClose={handleClose}
        MenuListProps={{
          "aria-labelledby": "basic-button",
        }}
      >
        <MenuItem onClick={() => handleShowAlignment(selectedData)}>
          Show Alignment
        </MenuItem>
        <MenuItem onClick={onDownload}>Save Selected Data</MenuItem>
      </Menu>

      {/* Modal for Alignment Visualization */}
      <Modal open={showModal} onClose={closeModal}>
        <Box
          sx={{
            position: "absolute",
            top: "50%",
            left: "50%",
            transform: "translate(-50%, -50%)",
            width: 800,
            maxHeight: "90vh", // Ensure it doesn't exceed viewport height
            bgcolor: "background.paper",
            boxShadow: 24,
            p: 4,
            borderRadius: 2,
            overflow: "auto", // Enable vertical scrolling if needed
          }}
        >
          <h2>Alignment Visualization</h2>
          {alignmentData &&
            alignmentData.map((alignment, index) => (
              <div
                key={index}
                style={{
                  marginBottom: "20px",
                  overflowX: "auto", // Enable horizontal scrolling
                  whiteSpace: "pre", // Preserve spacing for sequence alignments
                  border: "1px solid #ccc", // Optional: add a border for clarity
                  padding: "10px", // Add padding for better readability
                  borderRadius: "4px", // Rounded edges
                }}
              >
                <pre>
                  Query: {alignment.querySequence}
                  <br />
                  Match: {alignment.midline}
                  <br />
                  Hit:   {alignment.hitSequence}
                </pre>
              </div>
              ))}
            {!alignmentData && <p>No alignment data available.</p>}
          </Box>
      </Modal>

    </div>
  );
};

export default GearButton;
