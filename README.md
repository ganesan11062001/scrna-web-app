# Single-Cell RNA-seq Analysis Web App

Web application for analyzing single-cell RNA sequencing data.

## Tech Stack
- Frontend: Next.js (Deployed on Vercel)
- Backend: FastAPI (Deployed on Railway)
- Database: PostgreSQL (Railway)

## Features
- Upload CSV or H5AD files (max 50MB)
- Quality control and filtering
- UMAP visualization
- Clustering analysis (Leiden algorithm)
- Marker gene identification
- Interactive scatter plots

## Project Structure

scrna-app/
├── frontend/          # Next.js frontend
│   ├── pages/        # React pages
│   ├── components/   # React components
│   └── package.json
├── backend/          # FastAPI backend
│   ├── main.py      # API endpoints
│   ├── requirements.txt
│   └── railway.json
└── README.md

## Local Development

### Backend
cd backend
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
uvicorn main:app --reload

### Frontend
cd frontend
npm install
npm run dev

## Deployment
- Frontend: Vercel
- Backend: Railway

## License
MIT
