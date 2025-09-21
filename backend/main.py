from fastapi import FastAPI, File, UploadFile, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse, FileResponse
import pandas as pd
import scanpy as sc
import numpy as np
import json
import uuid
import io
from typing import Dict, Any
from sqlalchemy import create_engine, Column, String, JSON, DateTime, Integer
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
from datetime import datetime
import os
import tempfile
import csv

# FastAPI app
app = FastAPI(title="Single Cell RNA-seq Analysis")

# CORS for Vercel frontend
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Will update with actual Vercel URL later
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Database setup (Railway provides DATABASE_URL)
DATABASE_URL = os.getenv("DATABASE_URL", "sqlite:///./test.db")
if DATABASE_URL and DATABASE_URL.startswith("postgres://"):
    DATABASE_URL = DATABASE_URL.replace("postgres://", "postgresql://", 1)
    
engine = create_engine(DATABASE_URL)
SessionLocal = sessionmaker(bind=engine)
Base = declarative_base()

# Database Models
class Analysis(Base):
    __tablename__ = "analyses"
    
    id = Column(String, primary_key=True)
    created_at = Column(DateTime, default=datetime.utcnow)
    filename = Column(String)
    n_cells = Column(Integer)
    n_genes = Column(Integer)
    results = Column(JSON)
    status = Column(String, default="pending")

Base.metadata.create_all(bind=engine)

# Robust Analysis Pipeline with Comprehensive Error Handling
def run_analysis(file_content: bytes, filename: str) -> Dict[str, Any]:
    """
    Robust single-cell analysis pipeline with comprehensive error handling
    """
    try:
        print(f"Starting analysis for file: {filename}")
        print(f"File size: {len(file_content)} bytes")
        
        # Load data with comprehensive error handling
        if filename.endswith('.csv'):
            print("Loading CSV file...")
            try:
                df = pd.read_csv(io.BytesIO(file_content))
                print(f"CSV loaded: {df.shape[0]} rows, {df.shape[1]} columns")
                
                # Check if data is empty
                if df.empty:
                    raise ValueError("CSV file is empty")
                
                # Check if we have enough data
                if df.shape[0] < 10 or df.shape[1] < 10:
                    raise ValueError(f"Dataset too small: {df.shape[0]} cells, {df.shape[1]} genes. Need at least 10 cells and 10 genes.")
                
                # Check for non-numeric data
                if not df.select_dtypes(include=[np.number]).shape[1] > 0:
                    raise ValueError("CSV file must contain numeric data (gene expression values)")
                
                # Create AnnData object
                adata = sc.AnnData(df)
                print(f"AnnData created: {adata.n_obs} cells, {adata.n_vars} genes")
                
            except Exception as e:
                raise ValueError(f"Error loading CSV file: {str(e)}")
            
        elif filename.endswith('.h5ad'):
            print("Loading H5AD file...")
            try:
                adata = sc.read_h5ad(io.BytesIO(file_content))
                print(f"H5AD loaded: {adata.n_obs} cells, {adata.n_vars} genes")
            except Exception as e:
                raise ValueError(f"Error loading H5AD file: {str(e)}")
            
        else:
            print("Unsupported format, using demo data...")
            adata = sc.datasets.pbmc3k()
            print(f"Demo data loaded: {adata.n_obs} cells, {adata.n_vars} genes")
        
        # Limit size for free hosting
        if adata.n_obs > 5000:
            sc.pp.subsample(adata, n_obs=5000)
        
        results = {
            'n_cells': int(adata.n_obs),
            'n_genes': int(adata.n_vars)
        }
        
        # Enhanced QC metrics
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        adata.var['ribo'] = adata.var_names.str.match('^RP[SL]')
        adata.var['hb'] = adata.var_names.str.contains('^HB[^(P)]')
        
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'], inplace=True)
        
        # Calculate additional QC metrics
        n_cells_before = adata.n_obs
        n_genes_before = adata.n_vars
        
        results['qc'] = {
            'median_genes': float(np.median(adata.obs['n_genes_by_counts'])),
            'median_counts': float(np.median(adata.obs['total_counts'])),
            'median_mt': float(np.median(adata.obs['pct_counts_mt'])),
            'median_ribo': float(np.median(adata.obs['pct_counts_ribo'])),
            'median_hb': float(np.median(adata.obs['pct_counts_hb'])),
            'cells_before_filter': int(n_cells_before),
            'genes_before_filter': int(n_genes_before)
        }
        
        # Filter cells and genes
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        adata = adata[adata.obs['pct_counts_mt'] < 5, :]
        
        # Update results with post-filtering metrics
        results['qc']['cells_after_filter'] = int(adata.n_obs)
        results['qc']['genes_after_filter'] = int(adata.n_vars)
        results['qc']['cells_removed'] = int(n_cells_before - adata.n_obs)
        results['qc']['genes_removed'] = int(n_genes_before - adata.n_vars)
        
        # Normalize and log transform
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # Find highly variable genes
        sc.pp.highly_variable_genes(adata, n_top_genes=1000)
        adata = adata[:, adata.var.highly_variable]
        
        # PCA
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver='arpack')
        
        # Store PCA variance explained
        pca_variance = adata.uns['pca']['variance_ratio']
        results['pca_variance_explained'] = {
            'pc1': float(pca_variance[0]),
            'pc2': float(pca_variance[1]),
            'pc3': float(pca_variance[2]),
            'total_top10': float(np.sum(pca_variance[:10]))
        }
        
        # Compute neighborhood graph
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        
        # UMAP and clustering
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=0.5)
        
        results['n_clusters'] = len(adata.obs['leiden'].unique())
        
        # Get PCA coordinates for visualization
        pca_coords = adata.obsm['X_pca']
        if len(pca_coords) > 1000:
            idx = np.random.choice(len(pca_coords), 1000, replace=False)
            pca_coords = pca_coords[idx]
            pca_clusters = adata.obs['leiden'].iloc[idx].tolist()
        else:
            pca_clusters = adata.obs['leiden'].tolist()
        
        results['pca'] = {
            'x': pca_coords[:, 0].tolist(),
            'y': pca_coords[:, 1].tolist(),
            'clusters': pca_clusters
        }
        
        # Get UMAP coordinates for visualization
        umap_coords = adata.obsm['X_umap']
        # Sample if too many cells
        if len(umap_coords) > 1000:
            idx = np.random.choice(len(umap_coords), 1000, replace=False)
            umap_coords = umap_coords[idx]
            clusters = adata.obs['leiden'].iloc[idx].tolist()
        else:
            clusters = adata.obs['leiden'].tolist()
        
        results['umap'] = {
            'x': umap_coords[:, 0].tolist(),
            'y': umap_coords[:, 1].tolist(),
            'clusters': clusters
        }
        
        # Find marker genes
        sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
        
        # Get top 5 genes per cluster with scores
        markers = {}
        marker_scores = {}
        for cluster in adata.obs['leiden'].unique():
            genes = adata.uns['rank_genes_groups']['names'][cluster][:5]
            scores = adata.uns['rank_genes_groups']['scores'][cluster][:5]
            markers[str(cluster)] = genes.tolist()
            marker_scores[str(cluster)] = scores.tolist()
        
        results['markers'] = markers
        results['marker_scores'] = marker_scores
        
        # Calculate cluster statistics
        cluster_stats = {}
        for cluster in adata.obs['leiden'].unique():
            cluster_cells = adata.obs['leiden'] == cluster
            cluster_stats[str(cluster)] = {
                'n_cells': int(np.sum(cluster_cells)),
                'percentage': float(np.sum(cluster_cells) / len(cluster_cells) * 100)
            }
        
        results['cluster_stats'] = cluster_stats
        
        return results
        
    except Exception as e:
        print(f"Analysis error: {str(e)}")
        raise

# API Endpoints
@app.get("/")
async def root():
    return {"message": "Single Cell RNA-seq Analysis API"}

@app.post("/api/upload")
async def upload_and_analyze(file: UploadFile = File(...)):
    """
    Upload file and run analysis
    """
    # Check file size (limit to 50MB for free hosting)
    contents = await file.read()
    if len(contents) > 50 * 1024 * 1024:
        raise HTTPException(400, "File too large. Maximum 50MB.")
    
    # Create analysis record
    analysis_id = str(uuid.uuid4())
    db = SessionLocal()
    
    try:
        # Run analysis
        results = run_analysis(contents, file.filename)
        
        # Save to database
        analysis = Analysis(
            id=analysis_id,
            filename=file.filename,
            n_cells=results['n_cells'],
            n_genes=results['n_genes'],
            results=results,
            status="completed"
        )
        db.add(analysis)
        db.commit()
        
        return {
            "analysis_id": analysis_id,
            "status": "completed",
            "results": results
        }
        
    except Exception as e:
        # Save error with detailed information
        error_message = str(e)
        print(f"Analysis failed for {file.filename}: {error_message}")
        
        analysis = Analysis(
            id=analysis_id,
            filename=file.filename,
            status="failed",
            results={"error": error_message, "error_type": type(e).__name__}
        )
        db.add(analysis)
        db.commit()
        
        # Provide more helpful error messages
        if "empty" in error_message.lower():
            raise HTTPException(400, "The uploaded file appears to be empty. Please check your file and try again.")
        elif "too small" in error_message.lower():
            raise HTTPException(400, "The dataset is too small for analysis. Please ensure you have at least 10 cells and 10 genes.")
        elif "numeric data" in error_message.lower():
            raise HTTPException(400, "The CSV file must contain numeric gene expression data. Please check your file format.")
        elif "format" in error_message.lower():
            raise HTTPException(400, "Unsupported file format. Please upload a CSV or H5AD file.")
        else:
            raise HTTPException(500, f"Analysis failed: {error_message}")
    
    finally:
        db.close()

@app.get("/api/analysis/{analysis_id}")
async def get_analysis(analysis_id: str):
    """
    Get analysis results
    """
    db = SessionLocal()
    analysis = db.query(Analysis).filter(Analysis.id == analysis_id).first()
    db.close()
    
    if not analysis:
        raise HTTPException(404, "Analysis not found")
    
    return {
        "id": analysis.id,
        "status": analysis.status,
        "filename": analysis.filename,
        "n_cells": analysis.n_cells,
        "n_genes": analysis.n_genes,
        "results": analysis.results,
        "created_at": analysis.created_at
    }

@app.get("/api/analyses")
async def list_analyses():
    """
    List all analyses
    """
    db = SessionLocal()
    analyses = db.query(Analysis).order_by(Analysis.created_at.desc()).limit(20).all()
    db.close()
    
    return [
        {
            "id": a.id,
            "filename": a.filename,
            "status": a.status,
            "created_at": a.created_at,
            "n_cells": a.n_cells,
            "n_genes": a.n_genes
        }
        for a in analyses
    ]

@app.get("/api/sample-datasets")
async def get_sample_datasets():
    """
    Get list of available sample datasets
    """
    sample_datasets = [
        {
            "id": "pbmc3k",
            "name": "PBMC 3K Dataset",
            "description": "3,000 peripheral blood mononuclear cells from a healthy donor",
            "cells": 2700,
            "genes": 32738,
            "source": "10X Genomics"
        },
        {
            "id": "pbmc68k", 
            "name": "PBMC 68K Dataset",
            "description": "68,000 peripheral blood mononuclear cells from a healthy donor",
            "cells": 68579,
            "genes": 32738,
            "source": "10X Genomics"
        },
        {
            "id": "t_4k",
            "name": "T Cell 4K Dataset", 
            "description": "4,000 T cells from peripheral blood",
            "cells": 4000,
            "genes": 2000,
            "source": "10X Genomics"
        }
    ]
    return {"datasets": sample_datasets}

@app.post("/api/load-sample/{dataset_id}")
async def load_sample_dataset(dataset_id: str):
    """
    Load a sample dataset for analysis
    """
    try:
        # Create analysis record
        analysis_id = str(uuid.uuid4())
        db = SessionLocal()
        
        # Load sample data based on ID
        if dataset_id == "pbmc3k":
            adata = sc.datasets.pbmc3k()
        elif dataset_id == "pbmc68k":
            adata = sc.datasets.pbmc68k_reduced()
        elif dataset_id == "t_4k":
            adata = sc.datasets.t_4k()
        else:
            raise HTTPException(404, "Sample dataset not found")
        
        # Run analysis on sample data
        results = run_analysis_on_adata(adata, f"sample_{dataset_id}")
        
        # Save to database
        analysis = Analysis(
            id=analysis_id,
            filename=f"sample_{dataset_id}.h5ad",
            n_cells=results['n_cells'],
            n_genes=results['n_genes'],
            results=results,
            status="completed"
        )
        db.add(analysis)
        db.commit()
        db.close()
        
        return {
            "analysis_id": analysis_id,
            "status": "completed",
            "results": results
        }
        
    except Exception as e:
        raise HTTPException(500, f"Failed to load sample dataset: {str(e)}")

def run_analysis_on_adata(adata, dataset_name: str) -> Dict[str, Any]:
    """
    Run analysis pipeline on an AnnData object
    """
    try:
        # Limit size for free hosting
        if adata.n_obs > 5000:
            sc.pp.subsample(adata, n_obs=5000)
        
        results = {
            'n_cells': int(adata.n_obs),
            'n_genes': int(adata.n_vars)
        }
        
        # Enhanced QC metrics
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        adata.var['ribo'] = adata.var_names.str.match('^RP[SL]')
        adata.var['hb'] = adata.var_names.str.contains('^HB[^(P)]')
        
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribo', 'hb'], inplace=True)
        
        # Calculate additional QC metrics
        n_cells_before = adata.n_obs
        n_genes_before = adata.n_vars
        
        results['qc'] = {
            'median_genes': float(np.median(adata.obs['n_genes_by_counts'])),
            'median_counts': float(np.median(adata.obs['total_counts'])),
            'median_mt': float(np.median(adata.obs['pct_counts_mt'])),
            'median_ribo': float(np.median(adata.obs['pct_counts_ribo'])),
            'median_hb': float(np.median(adata.obs['pct_counts_hb'])),
            'cells_before_filter': int(n_cells_before),
            'genes_before_filter': int(n_genes_before)
        }
        
        # Filter cells and genes
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        adata = adata[adata.obs['pct_counts_mt'] < 5, :]
        
        # Update results with post-filtering metrics
        results['qc']['cells_after_filter'] = int(adata.n_obs)
        results['qc']['genes_after_filter'] = int(adata.n_vars)
        results['qc']['cells_removed'] = int(n_cells_before - adata.n_obs)
        results['qc']['genes_removed'] = int(n_genes_before - adata.n_vars)
        
        # Normalize and log transform
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # Find highly variable genes
        sc.pp.highly_variable_genes(adata, n_top_genes=1000)
        adata = adata[:, adata.var.highly_variable]
        
        # PCA
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver='arpack')
        
        # Store PCA variance explained
        pca_variance = adata.uns['pca']['variance_ratio']
        results['pca_variance_explained'] = {
            'pc1': float(pca_variance[0]),
            'pc2': float(pca_variance[1]),
            'pc3': float(pca_variance[2]),
            'total_top10': float(np.sum(pca_variance[:10]))
        }
        
        # Compute neighborhood graph
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        
        # UMAP and clustering
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=0.5)
        
        results['n_clusters'] = len(adata.obs['leiden'].unique())
        
        # Get PCA coordinates for visualization
        pca_coords = adata.obsm['X_pca']
        if len(pca_coords) > 1000:
            idx = np.random.choice(len(pca_coords), 1000, replace=False)
            pca_coords = pca_coords[idx]
            pca_clusters = adata.obs['leiden'].iloc[idx].tolist()
        else:
            pca_clusters = adata.obs['leiden'].tolist()
        
        results['pca'] = {
            'x': pca_coords[:, 0].tolist(),
            'y': pca_coords[:, 1].tolist(),
            'clusters': pca_clusters
        }
        
        # Get UMAP coordinates for visualization
        umap_coords = adata.obsm['X_umap']
        # Sample if too many cells
        if len(umap_coords) > 1000:
            idx = np.random.choice(len(umap_coords), 1000, replace=False)
            umap_coords = umap_coords[idx]
            clusters = adata.obs['leiden'].iloc[idx].tolist()
        else:
            clusters = adata.obs['leiden'].tolist()
        
        results['umap'] = {
            'x': umap_coords[:, 0].tolist(),
            'y': umap_coords[:, 1].tolist(),
            'clusters': clusters
        }
        
        # Find marker genes
        sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
        
        # Get top 5 genes per cluster with scores
        markers = {}
        marker_scores = {}
        for cluster in adata.obs['leiden'].unique():
            genes = adata.uns['rank_genes_groups']['names'][cluster][:5]
            scores = adata.uns['rank_genes_groups']['scores'][cluster][:5]
            markers[str(cluster)] = genes.tolist()
            marker_scores[str(cluster)] = scores.tolist()
        
        results['markers'] = markers
        results['marker_scores'] = marker_scores
        
        # Calculate cluster statistics
        cluster_stats = {}
        for cluster in adata.obs['leiden'].unique():
            cluster_cells = adata.obs['leiden'] == cluster
            cluster_stats[str(cluster)] = {
                'n_cells': int(np.sum(cluster_cells)),
                'percentage': float(np.sum(cluster_cells) / len(cluster_cells) * 100)
            }
        
        results['cluster_stats'] = cluster_stats
        
        return results
        
    except Exception as e:
        print(f"Analysis error: {str(e)}")
        raise

@app.get("/api/export/{analysis_id}")
async def export_analysis(analysis_id: str, format: str = "json"):
    """
    Export analysis results in various formats
    """
    db = SessionLocal()
    analysis = db.query(Analysis).filter(Analysis.id == analysis_id).first()
    db.close()
    
    if not analysis:
        raise HTTPException(404, "Analysis not found")
    
    if analysis.status != "completed":
        raise HTTPException(400, "Analysis not completed")
    
    try:
        if format == "json":
            # Export as JSON
            export_data = {
                "analysis_id": analysis.id,
                "filename": analysis.filename,
                "created_at": analysis.created_at.isoformat(),
                "n_cells": analysis.n_cells,
                "n_genes": analysis.n_genes,
                "results": analysis.results
            }
            
            # Create temporary file
            with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
                json.dump(export_data, f, indent=2)
                temp_path = f.name
            
            return FileResponse(
                temp_path,
                media_type='application/json',
                filename=f"analysis_{analysis_id}.json"
            )
            
        elif format == "csv":
            # Export marker genes as CSV
            results = analysis.results
            if 'markers' not in results:
                raise HTTPException(400, "No marker genes found")
            
            with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
                writer = csv.writer(f)
                writer.writerow(['Cluster', 'Gene', 'Score'])
                
                for cluster, genes in results['markers'].items():
                    scores = results.get('marker_scores', {}).get(cluster, [])
                    for i, gene in enumerate(genes):
                        score = scores[i] if i < len(scores) else 'N/A'
                        writer.writerow([cluster, gene, score])
                
                temp_path = f.name
            
            return FileResponse(
                temp_path,
                media_type='text/csv',
                filename=f"marker_genes_{analysis_id}.csv"
            )
            
        else:
            raise HTTPException(400, "Unsupported export format")
            
    except Exception as e:
        raise HTTPException(500, f"Export failed: {str(e)}")

@app.get("/api/test-upload")
async def test_upload():
    """
    Test endpoint to check if the analysis pipeline works with sample data
    """
    try:
        # Test with sample data
        adata = sc.datasets.pbmc3k()
        results = run_analysis_on_adata(adata, "test_sample")
        return {
            "status": "success",
            "message": "Analysis pipeline is working correctly",
            "test_results": {
                "n_cells": results['n_cells'],
                "n_genes": results['n_genes'],
                "n_clusters": results['n_clusters']
            }
        }
    except Exception as e:
        return {
            "status": "error",
            "message": f"Analysis pipeline test failed: {str(e)}"
        }

@app.get("/health")
async def health():
    return {"status": "healthy"}
    