from fastapi import FastAPI, File, UploadFile, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
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

# Simple Analysis Pipeline
def run_analysis(file_content: bytes, filename: str) -> Dict[str, Any]:
    """
    Simple single-cell analysis pipeline
    """
    try:
        # Load data
        if filename.endswith('.csv'):
            df = pd.read_csv(io.BytesIO(file_content))
            adata = sc.AnnData(df)
        elif filename.endswith('.h5ad'):
            adata = sc.read_h5ad(io.BytesIO(file_content))
        else:
            # Create demo data if format not supported
            adata = sc.datasets.pbmc3k()
        
        # Limit size for free hosting
        if adata.n_obs > 5000:
            sc.pp.subsample(adata, n_obs=5000)
        
        results = {
            'n_cells': int(adata.n_obs),
            'n_genes': int(adata.n_vars)
        }
        
        # Basic QC
        adata.var['mt'] = adata.var_names.str.startswith('MT-')
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
        
        results['qc'] = {
            'median_genes': float(np.median(adata.obs['n_genes_by_counts'])),
            'median_counts': float(np.median(adata.obs['total_counts'])),
            'median_mt': float(np.median(adata.obs.pct_counts_mt']))
        }
        
        # Filter cells and genes
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        adata = adata[adata.obs.pct_counts_mt < 5, :]
        
        # Normalize and log transform
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        
        # Find highly variable genes
        sc.pp.highly_variable_genes(adata, n_top_genes=1000)
        adata = adata[:, adata.var.highly_variable]
        
        # PCA
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver='arpack')
        
        # Compute neighborhood graph
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        
        # UMAP and clustering
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=0.5)
        
        results['n_clusters'] = len(adata.obs['leiden'].unique())
        
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
        
        # Get top 5 genes per cluster
        markers = {}
        for cluster in adata.obs['leiden'].unique():
            genes = adata.uns['rank_genes_groups']['names'][cluster][:5]
            markers[str(cluster)] = genes.tolist()
        
        results['markers'] = markers
        
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
        # Save error
        analysis = Analysis(
            id=analysis_id,
            filename=file.filename,
            status="failed",
            results={"error": str(e)}
        )
        db.add(analysis)
        db.commit()
        
        raise HTTPException(500, f"Analysis failed: {str(e)}")
    
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

@app.get("/health")
async def health():
    return {"status": "healthy"}
    