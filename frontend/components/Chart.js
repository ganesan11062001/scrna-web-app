import React from 'react';
import { Scatter } from 'react-chartjs-2';
import {
  Chart as ChartJS,
  LinearScale,
  PointElement,
  Tooltip,
  Legend
} from 'chart.js';

ChartJS.register(LinearScale, PointElement, Tooltip, Legend);

export default function Chart({ data }) {
  if (!data.umap) return null;

  const scatterData = {
    datasets: [{
      label: 'Cells',
      data: data.umap.x.map((x, i) => ({
        x: x,
        y: data.umap.y[i]
      })),
      backgroundColor: data.umap.clusters.map(c => 
        `hsl(${parseInt(c) * 40}, 70%, 50%)`
      ),
      pointRadius: 3,
    }]
  };

  return (
    <div>
      <h3>UMAP Visualization</h3>
      <div style={{ height: '400px' }}>
        <Scatter 
          data={scatterData}
          options={{
            responsive: true,
            maintainAspectRatio: false,
            scales: {
              x: { title: { display: true, text: 'UMAP 1' }},
              y: { title: { display: true, text: 'UMAP 2' }}
            }
          }}
        />
      </div>
    </div>
  );
}
