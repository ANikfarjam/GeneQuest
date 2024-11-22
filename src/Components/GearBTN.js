import React, { useState } from "react";
import { Menu, MenuItem, IconButton } from "@mui/material";
import SettingsIcon from "@mui/icons-material/Settings";
// import axios from 'axios';

const GearButton = ({ onDownload }) => {
  const [anchorEl, setAnchorEl] = useState(null);
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
  // Add inside GearButton component
  const handleDataAnalytics = async () => {
    try {
      if (onDataAnalytics) {
        await onDataAnalytics(); // Trigger the data analytics logic in parent
      }
      handleClose();
    } catch (error) {
      console.error("Error sending data to server:", error);
    }
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

      {/* Menu */}
      <Menu
        id="gear-menu"
        anchorEl={anchorEl}
        open={open}
        onClose={handleClose}
        MenuListProps={{
          "aria-labelledby": "basic-button",
        }}
      >
        {/* Menu Options */}
        <MenuItem onClick={handleDataAnalytics}>DataAnalytics</MenuItem>;
        <MenuItem onClick={handleDownload}>Save Selected Data</MenuItem>
        <MenuItem onClick={handleClose}>Blast</MenuItem>
      </Menu>
    </div>
  );
};

export default GearButton;
