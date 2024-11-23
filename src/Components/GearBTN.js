import React, { useState } from "react";
import { Menu, MenuItem, IconButton } from "@mui/material";
import SettingsIcon from "@mui/icons-material/Settings";

const GearButton = ({ onDownload }) => {
  const [anchorEl, setAnchorEl] = useState(null);
  const open = Boolean(anchorEl);

  const handleClick = (event) => {
    setAnchorEl(event.currentTarget);
  };

  const handleClose = () => {
    setAnchorEl(null);
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
        <MenuItem onClick={() => { onDownload(); handleClose(); }}>Save Selected Data</MenuItem>
      </Menu>
    </div>
  );
};

export default GearButton;
