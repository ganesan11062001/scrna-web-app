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

const API_URL = process.env.NEXT_PUBLIC_API_URL || 'http://localhost:8000';

export default function Home() {
  const [file, setFile] = useState(null);
  const [loading, setLoading] = useState(false);
  const [results, setResults] = useState(null);
  const [error, setError] = useState(null);
  const [sampleDatasets, setSampleDatasets] = useState([]);
  const [showSamples, setShowSamples] = useState(false);
  const [currentAnalysisId, setCurrentAnalysisId] = useState(null);
  const [testStatus, setTestStatus] = useState(null);

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
      setCurrentAnalysisId(response.data.analysis_id);
    } catch (err) {
      console.error('Upload error:', err);
      const errorMessage = err.response?.data?.detail || err.message || 'Analysis failed';
      setError(errorMessage);
    } finally {
      setLoading(false);
    }
  };

  const loadSampleDatasets = async () => {
    try {
      const response = await axios.get(`${API_URL}/api/sample-datasets`);
      setSampleDatasets(response.data.datasets);
      setShowSamples(true);
    } catch (err) {
      setError('Failed to load sample datasets');
    }
  };

  const loadSampleDataset = async (datasetId) => {
    setLoading(true);
    setError(null);
    
    try {
      const response = await axios.post(`${API_URL}/api/load-sample/${datasetId}`);
      setResults(response.data.results);
      setCurrentAnalysisId(response.data.analysis_id);
      setShowSamples(false);
    } catch (err) {
      setError(err.response?.data?.detail || 'Failed to load sample dataset');
    } finally {
      setLoading(false);
    }
  };

  const exportResults = async (format) => {
    if (!currentAnalysisId) {
      setError('No analysis to export');
      return;
    }

    try {
      const response = await axios.get(`${API_URL}/api/export/${currentAnalysisId}?format=${format}`, {
        responseType: 'blob'
      });
      
      // Create download link
      const url = window.URL.createObjectURL(new Blob([response.data]));
      const link = document.createElement('a');
      link.href = url;
      link.setAttribute('download', `analysis_${currentAnalysisId}.${format}`);
      document.body.appendChild(link);
      link.click();
      link.remove();
      window.URL.revokeObjectURL(url);
    } catch (err) {
      setError('Export failed');
    }
  };

  const testAnalysisPipeline = async () => {
    try {
      // First test basic connectivity
      console.log(`Testing connection to: ${API_URL}`);
      const healthResponse = await axios.get(`${API_URL}/health`);
      console.log('Health check:', healthResponse.data);
      
      // Then test the analysis pipeline
      const response = await axios.get(`${API_URL}/api/test-upload`);
      setTestStatus({
        ...response.data,
        connection: "Connected successfully"
      });
    } catch (err) {
      console.error('Test error:', err);
      let errorMessage = "Test failed";
      
      if (err.code === 'ECONNREFUSED' || err.message.includes('Network Error')) {
        errorMessage = `Cannot connect to backend at ${API_URL}. Make sure the backend is running on port 8000.`;
      } else if (err.response?.status === 404) {
        errorMessage = "Backend endpoint not found. Check if the backend is running the latest version.";
      } else {
        errorMessage = err.response?.data?.message || err.message || "Test failed";
      }
      
      setTestStatus({
        status: "error",
        message: errorMessage,
        connection: "Connection failed"
      });
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
          <div className="results-header">
            <h2>Analysis Results</h2>
            <div className="export-buttons">
              <button onClick={() => exportResults('json')} className="export-btn">
                Export JSON
              </button>
              <button onClick={() => exportResults('csv')} className="export-btn">
                Export CSV
              </button>
            </div>
          </div>
        
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
          <h3>Quality Control Metrics</h3>
          <div className="qc-grid">
            <div className="qc-section">
              <h4>Gene Expression</h4>
              <ul>
                <li>Median genes per cell: {results.qc.median_genes.toFixed(0)}</li>
                <li>Median counts per cell: {results.qc.median_counts.toFixed(0)}</li>
              </ul>
            </div>
            <div className="qc-section">
              <h4>Gene Type Percentages</h4>
              <ul>
                <li>Mitochondrial genes: {results.qc.median_mt.toFixed(2)}%</li>
                <li>Ribosomal genes: {results.qc.median_ribo?.toFixed(2) || 'N/A'}%</li>
                <li>Hemoglobin genes: {results.qc.median_hb?.toFixed(2) || 'N/A'}%</li>
              </ul>
            </div>
            <div className="qc-section">
              <h4>Filtering Results</h4>
              <ul>
                <li>Cells before filter: {results.qc.cells_before_filter}</li>
                <li>Cells after filter: {results.qc.cells_after_filter}</li>
                <li>Cells removed: {results.qc.cells_removed}</li>
                <li>Genes removed: {results.qc.genes_removed}</li>
              </ul>
            </div>
          </div>
        </div>

        <div className="visualizations">
          <h3>Dimensionality Reduction</h3>
          <div className="plot-grid">
            <div className="plot-container">
              <h4>PCA Plot</h4>
              {results.pca && (
                <div style={{ width: '500px', height: '500px' }}>
                  <Scatter 
                    data={{
                      datasets: [{
                        label: 'Cells',
                        data: results.pca.x.map((x, i) => ({
                          x: x,
                          y: results.pca.y[i]
                        })),
                        backgroundColor: results.pca.clusters.map(c => 
                          `hsl(${parseInt(c) * 40}, 70%, 50%)`
                        ),
                      }]
                    }}
                    options={{
                      responsive: true,
                      plugins: {
                        title: {
                          display: true,
                          text: 'PCA - Cell Clusters'
                        }
                      },
                      scales: {
                        x: {
                          title: {
                            display: true,
                            text: `PC1 (${(results.pca_variance_explained?.pc1 * 100).toFixed(1)}%)`
                          }
                        },
                        y: {
                          title: {
                            display: true,
                            text: `PC2 (${(results.pca_variance_explained?.pc2 * 100).toFixed(1)}%)`
                          }
                        }
                      }
                    }}
                  />
                </div>
              )}
            </div>
            
            <div className="plot-container">
              <h4>UMAP Plot</h4>
              <div style={{ width: '500px', height: '500px' }}>
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
          </div>
        </div>

        <div className="cluster-analysis">
          <h3>Cluster Analysis</h3>
          
          {results.cluster_stats && (
            <div className="cluster-stats">
              <h4>Cluster Statistics</h4>
              <div className="cluster-stats-grid">
                {Object.entries(results.cluster_stats).map(([cluster, stats]) => (
                  <div key={cluster} className="cluster-stat-card">
                    <h5>Cluster {cluster}</h5>
                    <p>Cells: {stats.n_cells}</p>
                    <p>Percentage: {stats.percentage.toFixed(1)}%</p>
                  </div>
                ))}
              </div>
            </div>
          )}

          <div className="markers">
            <h4>Top Marker Genes per Cluster</h4>
            {Object.entries(results.markers).map(([cluster, genes]) => (
              <div key={cluster} className="cluster-markers">
                <h5>Cluster {cluster}</h5>
                <div className="marker-genes">
                  {genes.map((gene, idx) => (
                    <div key={gene} className="marker-gene">
                      <span className="gene-name">{gene}</span>
                      {results.marker_scores && (
                        <span className="gene-score">
                          Score: {results.marker_scores[cluster][idx]?.toFixed(2) || 'N/A'}
                        </span>
                      )}
                    </div>
                  ))}
                </div>
              </div>
            ))}
          </div>
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
        
        <div className="divider">
          <span>OR</span>
        </div>
        
        <div className="sample-section">
          <h3>Try Sample Datasets</h3>
          <p>Explore the app with pre-loaded datasets</p>
          
          <div className="test-section">
            <button onClick={testAnalysisPipeline} className="test-button">
              Test Analysis Pipeline
            </button>
            {testStatus && (
              <div className={`test-status ${testStatus.status}`}>
                <strong>Test Result:</strong> {testStatus.message}
                {testStatus.connection && (
                  <div className="connection-info">
                    <strong>Connection:</strong> {testStatus.connection}
                  </div>
                )}
                {testStatus.test_results && (
                  <div className="test-details">
                    <p>Cells: {testStatus.test_results.n_cells}</p>
                    <p>Genes: {testStatus.test_results.n_genes}</p>
                    <p>Clusters: {testStatus.test_results.n_clusters}</p>
                  </div>
                )}
              </div>
            )}
          </div>
          
          {!showSamples ? (
            <button onClick={loadSampleDatasets} className="sample-button">
              View Sample Datasets
            </button>
          ) : (
            <div className="sample-datasets">
              {sampleDatasets.map((dataset) => (
                <div key={dataset.id} className="sample-dataset-card">
                  <h4>{dataset.name}</h4>
                  <p>{dataset.description}</p>
                  <div className="dataset-stats">
                    <span>Cells: {dataset.cells.toLocaleString()}</span>
                    <span>Genes: {dataset.genes.toLocaleString()}</span>
                    <span>Source: {dataset.source}</span>
                  </div>
                  <button 
                    onClick={() => loadSampleDataset(dataset.id)}
                    disabled={loading}
                    className="load-sample-btn"
                  >
                    {loading ? 'Loading...' : 'Load Dataset'}
                  </button>
                </div>
              ))}
              <button 
                onClick={() => setShowSamples(false)}
                className="back-button"
              >
                Back to Upload
              </button>
            </div>
          )}
        </div>
        
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

        .divider {
          text-align: center;
          margin: 20px 0;
          position: relative;
        }

        .divider::before {
          content: '';
          position: absolute;
          top: 50%;
          left: 0;
          right: 0;
          height: 1px;
          background: #ddd;
        }

        .divider span {
          background: #f5f5f5;
          padding: 0 15px;
          color: #666;
          font-weight: bold;
        }

        .sample-section {
          margin-top: 20px;
          padding-top: 20px;
          border-top: 1px solid #ddd;
        }

        .sample-section h3 {
          margin: 0 0 10px 0;
          color: #333;
        }

        .sample-section p {
          margin: 0 0 15px 0;
          color: #666;
        }

        .sample-button {
          background: #28a745;
          color: white;
          border: none;
          padding: 12px 24px;
          border-radius: 5px;
          cursor: pointer;
          font-size: 16px;
        }

        .sample-button:hover {
          background: #218838;
        }

        .sample-datasets {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
          gap: 20px;
          margin-top: 20px;
        }

        .sample-dataset-card {
          background: white;
          padding: 20px;
          border-radius: 8px;
          border: 1px solid #ddd;
          box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }

        .sample-dataset-card h4 {
          margin: 0 0 10px 0;
          color: #0070f3;
        }

        .sample-dataset-card p {
          margin: 0 0 15px 0;
          color: #666;
          font-size: 14px;
        }

        .dataset-stats {
          display: flex;
          flex-direction: column;
          gap: 5px;
          margin-bottom: 15px;
        }

        .dataset-stats span {
          font-size: 12px;
          color: #888;
        }

        .load-sample-btn {
          background: #0070f3;
          color: white;
          border: none;
          padding: 8px 16px;
          border-radius: 4px;
          cursor: pointer;
          font-size: 14px;
          width: 100%;
        }

        .load-sample-btn:disabled {
          background: #ccc;
          cursor: not-allowed;
        }

        .back-button {
          background: #6c757d;
          color: white;
          border: none;
          padding: 10px 20px;
          border-radius: 5px;
          cursor: pointer;
          font-size: 14px;
          margin-top: 20px;
        }

        .test-section {
          margin-bottom: 20px;
          padding: 15px;
          background: #f8f9fa;
          border-radius: 5px;
          border: 1px solid #dee2e6;
        }

        .test-button {
          background: #17a2b8;
          color: white;
          border: none;
          padding: 8px 16px;
          border-radius: 4px;
          cursor: pointer;
          font-size: 14px;
          margin-bottom: 10px;
        }

        .test-button:hover {
          background: #138496;
        }

        .test-status {
          padding: 10px;
          border-radius: 4px;
          margin-top: 10px;
        }

        .test-status.success {
          background: #d4edda;
          color: #155724;
          border: 1px solid #c3e6cb;
        }

        .test-status.error {
          background: #f8d7da;
          color: #721c24;
          border: 1px solid #f5c6cb;
        }

        .test-details {
          margin-top: 8px;
          font-size: 12px;
        }

        .test-details p {
          margin: 2px 0;
        }

        .connection-info {
          margin-top: 8px;
          font-size: 12px;
          font-weight: normal;
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

        .results-header {
          display: flex;
          justify-content: space-between;
          align-items: center;
          margin-bottom: 30px;
          flex-wrap: wrap;
          gap: 15px;
        }

        .export-buttons {
          display: flex;
          gap: 10px;
        }

        .export-btn {
          background: #28a745;
          color: white;
          border: none;
          padding: 8px 16px;
          border-radius: 4px;
          cursor: pointer;
          font-size: 14px;
        }

        .export-btn:hover {
          background: #218838;
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

        .qc-grid {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
          gap: 20px;
          margin-top: 15px;
        }

        .qc-section {
          background: white;
          padding: 15px;
          border-radius: 6px;
          border-left: 4px solid #0070f3;
        }

        .qc-section h4 {
          margin: 0 0 10px 0;
          color: #333;
          font-size: 14px;
        }

        .visualizations {
          margin: 30px 0;
        }

        .plot-grid {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(500px, 1fr));
          gap: 30px;
          margin-top: 20px;
        }

        .plot-container {
          background: #fafafa;
          padding: 20px;
          border-radius: 8px;
          text-align: center;
        }

        .plot-container h4 {
          margin: 0 0 15px 0;
          color: #333;
        }

        .cluster-analysis {
          margin: 30px 0;
        }

        .cluster-stats {
          margin-bottom: 30px;
        }

        .cluster-stats-grid {
          display: grid;
          grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
          gap: 15px;
          margin-top: 15px;
        }

        .cluster-stat-card {
          background: #f0f8ff;
          padding: 15px;
          border-radius: 6px;
          text-align: center;
          border: 2px solid #0070f3;
        }

        .cluster-stat-card h5 {
          margin: 0 0 8px 0;
          color: #0070f3;
          font-size: 16px;
        }

        .cluster-stat-card p {
          margin: 4px 0;
          font-size: 14px;
          color: #666;
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

        .cluster-markers h5 {
          margin: 0 0 10px 0;
          color: #0070f3;
        }

        .marker-genes {
          display: flex;
          flex-wrap: wrap;
          gap: 10px;
        }

        .marker-gene {
          background: white;
          padding: 8px 12px;
          border-radius: 4px;
          border: 1px solid #ddd;
          display: flex;
          flex-direction: column;
          align-items: center;
          min-width: 100px;
        }

        .gene-name {
          font-weight: bold;
          color: #333;
          margin-bottom: 4px;
        }

        .gene-score {
          font-size: 12px;
          color: #666;
        }
      `}</style>
    </div>
  );
}
