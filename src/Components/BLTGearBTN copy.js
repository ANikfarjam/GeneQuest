import React, { useState } from "react";
import { Menu, MenuItem, IconButton, Modal, Box } from "@mui/material";
import SettingsIcon from "@mui/icons-material/Settings";

const GearButton = ({ onDownload, selectedData }) => {
  const [anchorEl, setAnchorEl] = useState(null);
  const [showModal, setShowModal] = useState(false); // Alignment modal state
  const [alignmentData, setAlignmentData] = useState(null); // Alignment data
  const [showTreeModal, setShowTreeModal] = useState(false); // Tree modal state
  const [treeImgUrl, setTreeImgUrl] = useState(null); // Tree image URL

  const open = Boolean(anchorEl);

  const handleClick = (event) => {
    setAnchorEl(event.currentTarget);
  };

  const handleClose = () => {
    setAnchorEl(null);
  };

  // Alignment Visualization
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

  // Generate Phylogenetic Tree
  const handleShowTree = async (data) => {
    if (!data || data.length === 0) {
      alert("No data selected for tree generation!");
      return;
    }

    // Use the first accession as the "requestedQuery"
    const requestedQuery = data[0]?.accession || "";

    try {
      // Send selected data to the backend
      const response = await fetch("http://127.0.0.1:5000/gb_dat_anl", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          requestedQuery,
          sequences: data.map((item) => ({
            requestedQuerySequence: item.querySequence,
            accession: item.accession,
            hitSequence: item.hitSequence,
          })),
        }),
      });

      if (!response.ok) {
        throw new Error(`Error fetching tree image: ${response.statusText}`);
      }

      // Convert response to a blob and create an object URL
      const blob = await response.blob();
      const imgUrl = URL.createObjectURL(blob);
      setTreeImgUrl(imgUrl);
      setShowTreeModal(true); // Open the tree modal
    } catch (error) {
      console.error("Error generating phylogenetic tree:", error);
      alert("Failed to generate the phylogenetic tree.");
    }
  };

  const closeAlignmentModal = () => {
    setShowModal(false);
    setAlignmentData(null);
  };

  const closeTreeModal = () => {
    setShowTreeModal(false);
    setTreeImgUrl(null);
  };

  return (
    <div>
      {/* Gear Icon Button */}
      <IconButton
        aria-label="settings"
        aria-controls="gear-menu"
        aria-haspopup="true"
        onClick={handleClick}
      >
        <SettingsIcon />
      </IconButton>

      {/* Dropdown Menu */}
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
        <MenuItem onClick={() => handleShowTree(selectedData)}>
          Generate Phylogenetic Tree
        </MenuItem>
        <MenuItem onClick={onDownload}>Save Selected Data</MenuItem>
      </Menu>

      {/* Alignment Modal */}
      <Modal open={showModal} onClose={closeAlignmentModal}>
        <Box
          sx={{
            position: "absolute",
            top: "50%",
            left: "50%",
            transform: "translate(-50%, -50%)",
            width: 800,
            bgcolor: "background.paper",
            boxShadow: 24,
            p: 4,
            borderRadius: 2,
            overflow: "auto",
          }}
        >
          <h2>Alignment Visualization</h2>
          {alignmentData &&
            alignmentData.map((alignment, index) => (
              <div
                key={index}
                style={{
                  marginBottom: "20px",
                  overflowX: "auto",
                  whiteSpace: "pre",
                  border: "1px solid #ccc",
                  padding: "10px",
                  borderRadius: "4px",
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

      {/* Phylogenetic Tree Modal */}
      <Modal open={showTreeModal} onClose={closeTreeModal}>
        <Box
          sx={{
            position: "absolute",
            top: "50%",
            left: "50%",
            transform: "translate(-50%, -50%)",
            width: 800,
            bgcolor: "background.paper",
            boxShadow: 24,
            p: 4,
            borderRadius: 2,
            overflow: "auto",
            textAlign: "center",
          }}
        >
          <h2>Phylogenetic Tree</h2>
          {treeImgUrl ? (
            <img
              src={treeImgUrl}
              alt="Phylogenetic Tree"
              style={{ maxWidth: "100%", height: "auto" }}
            />
          ) : (
            <p>Loading...</p>
          )}
        </Box>
      </Modal>
    </div>
  );
};

export default GearButton;
