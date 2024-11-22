import React, { useState } from "react";
import { Menu, MenuItem, IconButton, Modal, Box } from "@mui/material";
import SettingsIcon from "@mui/icons-material/Settings";
import axios from "axios";
import BarChartConserved from "./BarChartConserved"; // Conserved regions chart
import { LogoViewer } from "logojs-react"; // Sequence Logo visualization

const GearButton = ({ onDownload, selectedData }) => {
  const [anchorEl, setAnchorEl] = useState(null);
  const [showModal, setShowModal] = useState(false); // Modal visibility state
  const [conservedData, setConservedData] = useState(null); // Conserved regions
  const [alignmentData, setAlignmentData] = useState(null); // Aligned sequence data
  const [loading, setLoading] = useState(false); // Loading state

  const open = Boolean(anchorEl);

  const handleClick = (event) => {
    setAnchorEl(event.currentTarget);
  };

  const handleClose = () => {
    setAnchorEl(null);
  };

  const handleDownload = () => {
    onDownload();
    handleClose();
  };

  const handleDataAnalytics = async () => {
    if (!selectedData || selectedData.length === 0) {
      alert("No data selected for analysis!");
      return;
    }

    // Open modal and set loading state
    setShowModal(true);
    setLoading(true);
    setConservedData(null);
    setAlignmentData(null);

    try {
      // Send selected data to backend for analysis
      const response = await axios.post(
        "http://127.0.0.1:5000/alignment_conservation/seq_logo",
        selectedData
      );

      console.log("Response from backend:", response.data); // Log the response data

      // Process response data
      setConservedData(response.data.conserved_regions);
      setAlignmentData(response.data.alignment);
    } catch (error) {
      console.error("Error during data analysis:", error);
    } finally {
      // Stop loading
      setLoading(false);
    }
  };

  const closeModal = () => {
    setShowModal(false);
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
        <MenuItem onClick={handleDataAnalytics}>Data Analytics</MenuItem>
        <MenuItem onClick={handleDownload}>Save Selected Data</MenuItem>
        <MenuItem onClick={handleClose}>Blast</MenuItem>
      </Menu>

      {/* Modal for Visualization */}
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
          {loading && <p>Loading... Please wait while we process your data.</p>}
          {!loading && conservedData && (
            <>
              <h3>Conserved Regions</h3>
              <BarChartConserved data={conservedData} />
            </>
          )}
          {!loading && alignmentData && (
            <>
              <h3>Aligned Sequences</h3>
              <LogoViewer
                data={{
                  sequences: alignmentData.map((d) => ({ sequence: d.sequence })),
                }}
              />
            </>
          )}
        </Box>
      </Modal>
    </div>
  );
};

export default GearButton;
