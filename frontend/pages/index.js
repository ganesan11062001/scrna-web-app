import React, { useState } from 'react';
import axios from 'axios';
import { Line, Scatter } from 'react-chartjs-2';
import {
  Chart as ChartJS,
  CategoryScale,
  LinearScale,
  PointElement,
  LineElement,
  Title,
  Tooltip,
  Legend
} from 'chart.js';

ChartJS.register(
  CategoryScale,
  LinearScale,
  PointElement,
  LineElement,
  Title,
  Tooltip,
  Legend
);

const API_URL = process.env.NEXT_PUBLIC_API_URL || 'https://your-app.railway.app';

export default function Home() {
  const [file, setFile] = useState(null);
  const [loading, setLoading] = useState(false);
  const [results, setResults] = useState(null);
  const [error, setError] = useState(null);

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
      const response = await axios.post(`${API_URL}/api/upload`, formData, {
        headers: {
          'Content-Type': 'multipart/form-data',
        },
      });
      
      setResults(response.data.results);
    } catch (err) {
      setError(err.response?.data?.detail || 'Analysis failed');
    } finally {
      setLoading(false);
    }
  };

  const renderResults = () => {
    if (!results) return null;

    // Prepare UMAP data
    const scatterData = {
      datasets: [{
        label: 'Cells',
        data: results.umap.x.map((x, i) => ({
          x: x,
          y: results.umap.y[i]
        })),
        backgroundColor: results.umap.clusters.map(c => 
          `hsl(${parseInt(c) * 40}, 70%, 50%)`
        ),
      }]
    };

    return (
      <div className="results-container">
        <h2>Analysis Results</h2>
        
        <div className="stats">
          <div className="stat-card">
            <h3>Total Cells</h3>
            <p>{results.n_cells}</p>
          </div>
          <div className="stat-card">
            <h3>Total Genes</h3>
            <p>{results.n_genes}</p>
          </div>
          <div className="stat-card">
            <h3>Clusters Found</h3>
            <p>{results.n_clusters}</p>
          </div>
        </div>

        <div className="qc-metrics">
          <h3>QC Metrics</h3>
          <ul>
            <li>Median genes per cell: {results.qc.median_genes.toFixed(0)}</li>
            <li>Median counts per cell: {results.qc.median_counts.toFixed(0)}</li>
            <li>Median MT%: {results.qc.median_mt.toFixed(2)}%</li>
          </ul>
        </div>

        <div className="umap-plot">
          <h3>UMAP Visualization</h3>
          <div style={{ width: '600px', height: '600px' }}>
            <Scatter 
              data={scatterData}
              options={{
                responsive: true,
                plugins: {
                  title: {
                    display: true,
                    text: 'UMAP - Cell Clusters'
                  }
                },
                scales: {
                  x: {
                    title: {
                      display: true,
                      text: 'UMAP 1'
                    }
                  },
                  y: {
                    title: {
                      display: true,
                      text: 'UMAP 2'
                    }
                  }
                }
              }}
            />
          </div>
        </div>

        <div className="markers">
          <h3>Top Marker Genes per Cluster</h3>
          {Object.entries(results.markers).map(([cluster, genes]) => (
            <div key={cluster} className="cluster-markers">
              <h4>Cluster {cluster}</h4>
              <ul>
                {genes.map(gene => (
                  <li key={gene}>{gene}</li>
                ))}
              </ul>
            </div>
          ))}
        </div>
      </div>
    );
  };

  return (
    <div className="container">
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
        
        {error && (
          <div className="error">
            Error: {error}
          </div>
        )}
      </div>

      {loading && (
        <div className="loading">
          <div className="spinner"></div>
          <p>Processing your data... This may take a few minutes.</p>
        </div>
      )}

      {results && renderResults()}

      <style jsx>{`
        .container {
          max-width: 1200px;
          margin: 0 auto;
          padding: 20px;
          font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
        }

        h1 {
          color: #333;
          text-align: center;
        }

        .upload-section {
          background: #f5f5f5;
          padding: 30px;
          border-radius: 10px;
          margin: 30px 0;
        }

        input[type="file"] {
          margin: 20px 0;
          display: block;
        }

        button {
          background: #0070f3;
          color: white;
          border: none;
          padding: 12px 24px;
          border-radius: 5px;
          cursor: pointer;
          font-size: 16px;
        }

        button:disabled {
          background: #ccc;
          cursor: not-allowed;
        }

        .error {
          color: red;
          margin-top: 10px;
        }

        .loading {
          text-align: center;
          margin: 50px 0;
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

        .results-container {
          margin-top: 40px;
        }

        .stats {
          display: grid;
          grid-template-columns: repeat(3, 1fr);
          gap: 20px;
          margin: 30px 0;
        }

        .stat-card {
          background: #f0f0f0;
          padding: 20px;
          border-radius: 8px;
          text-align: center;
        }

        .stat-card h3 {
          margin: 0;
          color: #555;
          font-size: 14px;
        }

        .stat-card p {
          margin: 10px 0 0 0;
          font-size: 28px;
          font-weight: bold;
          color: #0070f3;
        }

        .qc-metrics {
          background: #fafafa;
          padding: 20px;
          border-radius: 8px;
          margin: 20px 0;
        }

        .umap-plot {
          margin: 30px 0;
        }

        .markers {
          margin: 30px 0;
        }

        .cluster-markers {
          background: #f9f9f9;
          padding: 15px;
          margin: 10px 0;
          border-radius: 5px;
        }

        .cluster-markers h4 {
          margin: 0 0 10px 0;
          color: #0070f3;
        }

        .cluster-markers ul {
          margin: 0;
          padding-left: 20px;
        }
      `}</style>
    </div>
  );
}
