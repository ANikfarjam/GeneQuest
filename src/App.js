import React from 'react';
import './Components/styling/App.css';
import Navbar from './Components/Nav_Bar';
import GenomeSearch from './Components/NCBI';
import GenBankSearch from './Components/GenBankSearch';
import PhylogeneticTree from './Components/PhylogeneticTree';
import MegaBLAST from './Components/MegaBlast';
import { BrowserRouter as Router, Route, Routes } from 'react-router-dom';
import logo from './Components/media/logo.jpg';
import TwoCardsComponent from './Components/TwoCardsComponent'

function App() {
  return (
    <Router>
      <div className="App">
        <Navbar />
        <Routes>
          <Route path="/" element={
            <div style={{ display: 'flex', alignItems: 'center', padding: '20px' }}>
              <img src={logo} alt="Logo" style={{ width: '80px', marginRight: '20px' }} />
              <p>
                Welcome to our Genomic Exploration Hub! This web app empowers researchers and enthusiasts to dive deep into the world of genetic information. With easy access to GenBank searches, users can retrieve specific DNA sequences, analyze them through BLAST, and generate phylogenetic trees for comparative insights. Whether you're exploring evolutionary relationships or seeking precise genetic alignments, this platform provides a streamlined, interactive experience for genomic research and discovery. Start your journey to uncover the hidden patterns and connections in DNA sequences!
              </p>
            </div>
          } />
          <Route path="/GenBank" element={<GenBankSearch />} />
          <Route path="/Blast" element={<GenomeSearch />} />
          <Route path="/PhylogeneticTree" element={<PhylogeneticTree />} />
          <Route path="/MegaBLAST" element={<MegaBLAST />} />
          <Route path="/about" element={<TwoCardsComponent/>} />
        </Routes>
      </div>
    </Router>
  );
}

export default App;

