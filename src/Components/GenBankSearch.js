import React, { useState } from 'react';
import axios from 'axios';
import GBDataTable from './GBDataTable';

const GenBankSearch = () => {
  const [searchTerm, setSearchTerm] = useState('');
  const [searchResults, setSearchResults] = useState([]);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);

  const handleSearch = async () => {
    setLoading(true);
    setError(null);
    setSearchResults([]); // Reset previous results

    try {
      const response = await axios.post('http://127.0.0.1:5000/api/genbank_search', {
        query: searchTerm,
      });

      if (response.data.sequences) {
        // Split the sequences into structured data
        const sequences = response.data.sequences.split('\n>').map((seq, index) => {
          const lines = seq.trim().split('\n');
          const idAndDesc = lines[0].replace('>', '').split(' ', 2);
          const sequence = lines.slice(1).join('');
          return {
            id: idAndDesc[0] || `Sequence_${index + 1}`,
            description: idAndDesc[1] || 'No description available',
            sequence: sequence || 'No sequence data',
          };
        });
        setSearchResults(sequences);
      } else {
        setError('No sequences found.');
      }
    } catch (error) {
      console.error('Error searching GenBank:', error);
      setError('Failed to retrieve search results. Please try again.');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div>
      <h2>GenBank Search</h2>
      <input
        type="text"
        value={searchTerm}
        onChange={(e) => setSearchTerm(e.target.value)}
        placeholder="Enter search term"
      />
      <button onClick={handleSearch} disabled={loading}>
        {loading ? 'Searching...' : 'Search'}
      </button>

      {loading && <p>Loading... Please wait.</p>}
      {error && <p style={{ color: 'red' }}>{error}</p>}
      
      {searchResults.length > 0 && (
        <div>
          <h3>Search Results:</h3>
          <GBDataTable data={searchResults} />
        </div>
      )}
    </div>
  );
};

export default GenBankSearch;