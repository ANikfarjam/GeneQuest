import React, { useState } from 'react';
import axios from 'axios';
import xml2js from 'xml2js';
import DataTable from './DataTable';

function GenomeSearch() {
  const [inputSequence, setInputSequence] = useState('');
  const [data, setData] = useState([]);
  const [error, setError] = useState(null);
  const [loading, setLoading] = useState(false);

  const parser = new xml2js.Parser({ explicitArray: false });

  const generateCacheKey = (sequence) => `blast_result_${sequence}`;

  const loadFromCache = (key) => {
    const cachedData = localStorage.getItem(key);
    if (cachedData) {
      console.log('Using cached result from localStorage...');
      return JSON.parse(cachedData);
    }
    return null;
  };

  const saveToCache = (key, data) => {
    localStorage.setItem(key, JSON.stringify(data));
    console.log('Results saved to localStorage.');
  };

  const submitBlast = async (sequence) => {
    const response = await axios.post('https://blast.ncbi.nlm.nih.gov/Blast.cgi', null, {
      params: {
        CMD: 'Put',
        DATABASE: 'nt',
        PROGRAM: 'blastn',
        QUERY: sequence,
        FORMAT_TYPE: 'XML',
      },
    });

    const rid = response.data.match(/RID = (\w+)/)[1];
    const rtoe = parseInt(response.data.match(/RTOE = (\d+)/)[1], 10);
    return { rid, rtoe };
  };

  const checkBlastStatus = async (rid, rtoe) => {
    let status = 'WAITING';

    while (status === 'WAITING') {
      await new Promise((resolve) => setTimeout(resolve, rtoe * 1000)); // Wait based on RTOE
      const statusResponse = await axios.get('https://blast.ncbi.nlm.nih.gov/Blast.cgi', {
        params: {
          CMD: 'Get',
          FORMAT_OBJECT: 'SearchInfo',
          RID: rid,
        },
      });

      if (statusResponse.data.includes('Status=WAITING')) {
        console.log('Job is still processing...');
      } else if (statusResponse.data.includes('Status=FAILED')) {
        throw new Error('BLAST search failed.');
      } else {
        status = 'READY';
      }
    }

    const resultResponse = await axios.get('https://blast.ncbi.nlm.nih.gov/Blast.cgi', {
      params: {
        CMD: 'Get',
        FORMAT_TYPE: 'XML',
        RID: rid,
      },
    });

    return resultResponse.data;
  };

  const runBLASTSearch = async () => {
    try {
      setLoading(true);
      setError(null);
      setData([]);

      const cacheKey = generateCacheKey(inputSequence);

      // Check localStorage for cached results
      const cachedData = loadFromCache(cacheKey);
      if (cachedData) {
        setData(cachedData);
        return;
      }

      // Perform BLAST search
      const { rid, rtoe } = await submitBlast(inputSequence);
      const xmlResults = await checkBlastStatus(rid, rtoe);

      parser.parseString(xmlResults, (err, result) => {
        if (err) {
          setError('Error parsing XML results.');
          console.error('XML Parsing Error:', err);
          return;
        }

        const hits = result?.BlastOutput?.BlastOutput_iterations?.Iteration?.Iteration_hits?.Hit || [];
        const formattedHits = hits.map((hit) => ({
          num: hit.Hit_num,
          id: hit.Hit_id,
          definition: hit.Hit_def,
          accession: hit.Hit_accession,
          length: hit.Hit_len,
          bitScore: hit.Hit_hsps?.Hsp?.Hsp_bit_score || 'N/A',
          evalue: hit.Hit_hsps?.Hsp?.Hsp_evalue || 'N/A',
          identity: hit.Hit_hsps?.Hsp?.Hsp_identity || 'N/A',
          querySequence: hit.Hit_hsps?.Hsp?.Hsp_qseq || 'N/A',
          hitSequence: hit.Hit_hsps?.Hsp?.Hsp_hseq || 'N/A',
          midline: hit.Hit_hsps?.Hsp?.Hsp_midline || 'N/A',
        }));

        setData(formattedHits);

        // Save results to localStorage for future use
        saveToCache(cacheKey, formattedHits);
      });
    } catch (error) {
      setError(error.message || 'Error performing BLAST search.');
      console.error('Error performing BLAST search:', error);
    } finally {
      setLoading(false);
    }
  };

  const handleSearch = () => {
    if (inputSequence.trim() !== '') {
      runBLASTSearch();
    } else {
      setError('Please enter a valid sequence.');
    }
  };

  return (
    <div style={{ margin: '20px', padding: '20px', border: '1px solid #ccc', borderRadius: '8px', backgroundColor: '#f9f9f9' }}>
      <h2>Genome Search</h2>
      <input
        type="text"
        placeholder="Please enter your sequence"
        value={inputSequence}
        onChange={(e) => setInputSequence(e.target.value)}
        style={{ width: '80%', padding: '10px', marginBottom: '10px' }}
      />
      <button onClick={handleSearch} style={{ padding: '10px 20px', backgroundColor: '#0077b6', color: 'white', border: 'none', borderRadius: '4px' }}>
        Search Genome
      </button>

      {loading && <p>Loading... This may take several minutes. Please wait.</p>}
      {error && <p style={{ color: 'red' }}>Error: {error}</p>}

      {/* Render the DataTable instead of JSON */}
      <DataTable data={data} />
    </div>
  );

  
}

export default GenomeSearch;
