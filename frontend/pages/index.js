import React, { useState } from 'react';
import Head from 'next/head';
import dynamic from 'next/dynamic';

const Chart = dynamic(() => import('../components/Chart'), { ssr: false });

export default function Home() {
  const [file, setFile] = useState(null);
  const [loading, setLoading] = useState(false);
  const [results, setResults] = useState(null);
  const [error, setError] = useState(null);

  const API_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

  const handleFileChange = (e) => {
    setFile(e.target.files[0]);
    setError(null);
  };

  const handleUpload = async () => {
    if (!file) {
      setError('Please select a file');
      return;
    }

    setLoading(true);
    setError(null);
    
    const formData = new FormData();
    formData.append('file', file);

    try {
      const response = await fetch(`${API_URL}/api/upload`, {
        method: 'POST',
        body: formData,
      });
      
      if (!response.ok) {
        throw new Error('Analysis failed');
      }
      
      const data = await response.json();
      setResults(data.results);
    } catch (err) {
      setError(err.message || 'Analysis failed');
    } finally {
      setLoading(false);
    }
  };

  return (
    <div className="container">
      <Head>
        <title>Single-Cell RNA-seq Analysis</title>
        <meta name="description" content="Analyze single-cell RNA sequencing data" />
      </Head>

      <main>
        <h1>🧬 Single-Cell RNA-seq Analysis</h1>
        
        <div className="upload-section">
          <h2>Upload Your Data</h2>
          <p>Supported formats: .h5ad, .csv (max 50MB)</p>
          
          <input 
            type="file" 
            onChange={handleFileChange}
            accept=".h5ad,.csv"
          />
          
          <button 
            onClick={handleUpload} 
            disabled={loading || !file}
          >
            {loading ? 'Analyzing...' : 'Start Analysis'}
          </button>
          
          {error && <div className="error">Error: {error}</div>}
        </div>

        {loading && (
          <div className="loading">
            <div className="spinner"></div>
            <p>Processing your data...</p>
          </div>
        )}

        {results && (
          <div className="results">
            <h2>Results</h2>
            <div className="stats">
              <div className="stat">
                <h3>Cells</h3>
                <p>{results.n_cells}</p>
              </div>
              <div className="stat">
                <h3>Genes</h3>
                <p>{results.n_genes}</p>
              </div>
              <div className="stat">
                <h3>Clusters</h3>
                <p>{results.n_clusters}</p>
              </div>
            </div>
            {results.umap && <Chart data={results} />}
          </div>
        )}
      </main>

      <style jsx>{`
        .container {
          min-height: 100vh;
          padding: 2rem;
        }
        main {
          max-width: 1200px;
          margin: 0 auto;
        }
        h1 {
          text-align: center;
          color: #333;
        }
        .upload-section {
          background: #f5f5f5;
          padding: 2rem;
          border-radius: 10px;
          margin: 2rem 0;
        }
        input[type="file"] {
          display: block;
          margin: 1rem 0;
        }
        button {
          background: #0070f3;
          color: white;
          border: none;
          padding: 12px 24px;
          border-radius: 5px;
          cursor: pointer;
        }
        button:disabled {
          background: #ccc;
        }
        .error {
          color: red;
          margin-top: 1rem;
        }
        .loading {
          text-align: center;
          margin: 2rem 0;
        }
        .spinner {
          border: 4px solid #f3f3f3;
          border-top: 4px solid #0070f3;
          border-radius: 50%;
          width: 40px;
          height: 40px;
          animation: spin 1s linear infinite;
          margin: 0 auto;
        }
        @keyframes spin {
          0% { transform: rotate(0deg); }
          100% { transform: rotate(360deg); }
        }
        .stats {
          display: grid;
          grid-template-columns: repeat(3, 1fr);
          gap: 1rem;
          margin: 2rem 0;
        }
        .stat {
          background: #f0f0f0;
          padding: 1rem;
          border-radius: 8px;
          text-align: center;
        }
        .stat h3 {
          margin: 0;
          font-size: 0.9rem;
          color: #666;
        }
        .stat p {
          margin: 0.5rem 0 0;
          font-size: 2rem;
          font-weight: bold;
          color: #0070f3;
        }
      `}</style>
    </div>
  );
}
