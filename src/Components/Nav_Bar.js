import React from 'react';
import { Link } from 'react-router-dom';
import logo from './media/logo.jpg'; // Adjust this path if the logo is located elsewhere
import './styling/NavBar.css'

const Navbar = () => {
    return (
        <div className="navbar">
            <div className="navbar-left">
                <Link to="/" className="navbar-logo">
                    <img src={logo} alt="Logo" className="logo-image" />
                </Link>
                <Link to="/" className="nimbus-link">
                    <h3>GeneQuest</h3>
                </Link>
            </div>
            <div className="navbar-center">
                <div className="navbar-options">
                    <Link to="/">Home</Link>
                    <Link to="/GenBank">GenBank Search</Link>
                    <Link to="/Blast">Blast a Sequence</Link>
                    <Link to="/about">About Us</Link>
                </div>
            </div>
        </div>
    );
};

export default Navbar;
