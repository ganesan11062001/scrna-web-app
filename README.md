# Single-Cell RNA-seq Analysis Web App

A comprehensive web application for analyzing single-cell RNA sequencing data with advanced quality control, clustering, and visualization capabilities.

## 🚀 Features

### Core Analysis
- 📊 **File Upload**: Support for CSV and H5AD files (max 50MB)
- 🔬 **Quality Control**: Comprehensive QC metrics including mitochondrial, ribosomal, and hemoglobin gene percentages
- 📈 **Dimensionality Reduction**: PCA and UMAP visualization with variance explained
- 🎯 **Clustering**: Leiden clustering algorithm with configurable resolution
- 🧬 **Marker Gene Analysis**: Top marker genes per cluster with statistical scores

### Enhanced Features
- 📋 **Sample Datasets**: Pre-loaded datasets (PBMC 3K, PBMC 68K, T Cell 4K) for testing
- 📊 **Interactive Visualizations**: Side-by-side PCA and UMAP plots with cluster information
- 📈 **Detailed Statistics**: Cell and gene filtering statistics, cluster composition
- 💾 **Data Export**: Export results as JSON or marker genes as CSV
- 🎨 **Modern UI**: Responsive design with improved user experience

### Quality Control Metrics
- Gene expression statistics (median genes/counts per cell)
- Gene type percentages (mitochondrial, ribosomal, hemoglobin)
- Pre and post-filtering cell/gene counts
- Filtering efficiency metrics

## 🛠 Tech Stack

- **Frontend**: Next.js 14, React 18, Chart.js, Axios
- **Backend**: FastAPI, Scanpy, SQLAlchemy, PostgreSQL
- **Analysis**: Scanpy, NumPy, Pandas, SciPy
- **Deployment**: Vercel (frontend), Railway (backend)

## 🌐 Live Demo

- Frontend: https://scrna-web-app.vercel.app
- API: https://scrna-web-app-production.up.railway.app

## 🚀 Local Development

### Backend Setup
```bash
cd backend
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
pip install -r requirements.txt
uvicorn main:app --reload --port 8000
```

### Frontend Setup
```bash
cd frontend
npm install
npm run dev
```

## 📊 API Endpoints

### Analysis
- `POST /api/upload` - Upload and analyze data file
- `POST /api/load-sample/{dataset_id}` - Load sample dataset
- `GET /api/analysis/{analysis_id}` - Get analysis results
- `GET /api/analyses` - List all analyses

### Sample Data
- `GET /api/sample-datasets` - Get available sample datasets

### Export
- `GET /api/export/{analysis_id}?format=json` - Export full results as JSON
- `GET /api/export/{analysis_id}?format=csv` - Export marker genes as CSV

## 📈 Analysis Pipeline

1. **Data Loading**: Support for CSV and H5AD formats
2. **Quality Control**: 
   - Calculate QC metrics (genes, counts, gene type percentages)
   - Filter cells (min 200 genes) and genes (min 3 cells)
   - Remove cells with >5% mitochondrial genes
3. **Normalization**: Total count normalization and log transformation
4. **Feature Selection**: Identify top 1000 highly variable genes
5. **Dimensionality Reduction**: 
   - PCA with variance explained calculation
   - UMAP embedding
6. **Clustering**: Leiden clustering with resolution 0.5
7. **Marker Gene Analysis**: Wilcoxon rank-sum test for each cluster

## 🎯 Sample Datasets

- **PBMC 3K**: 2,700 peripheral blood mononuclear cells
- **PBMC 68K**: 68,579 peripheral blood mononuclear cells  
- **T Cell 4K**: 4,000 T cells from peripheral blood

## 📋 Output Files

### JSON Export
Complete analysis results including:
- QC metrics and filtering statistics
- PCA and UMAP coordinates
- Cluster assignments and statistics
- Marker genes with scores
- PCA variance explained

### CSV Export
Marker genes table with:
- Cluster ID
- Gene name
- Statistical score

## 🔧 Configuration

### Environment Variables
- `DATABASE_URL`: PostgreSQL connection string
- `NEXT_PUBLIC_API_URL`: Backend API URL for frontend

### Analysis Parameters
- Maximum cells: 5,000 (for free hosting)
- Minimum genes per cell: 200
- Minimum cells per gene: 3
- Maximum mitochondrial percentage: 5%
- Highly variable genes: 1,000
- Clustering resolution: 0.5

## 📚 Dependencies

### Backend
- fastapi==0.104.0
- scanpy==1.9.5
- pandas==2.1.0
- numpy==1.24.0
- sqlalchemy==2.0.0
- psycopg2-binary==2.9.0
- leidenalg==0.10.0

### Frontend
- next==14.0.0
- react==18.2.0
- chart.js==4.4.0
- react-chartjs-2==5.2.0
- axios==1.5.0

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Test thoroughly
5. Submit a pull request

## 📄 License

This project is licensed under the MIT License.

## 🙏 Acknowledgments

- Scanpy team for the excellent single-cell analysis library
- 10X Genomics for sample datasets
- FastAPI and Next.js communities
