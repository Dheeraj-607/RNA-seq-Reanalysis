// Initialize UMAP Visualization
let umapData = [];
let geneExpressionData = {};
let allGenes = [];

// Load UMAP data from JSON file
fetch('data/umap_clusters.json')
    .then(response => response.json())
    .then(data => {
        umapData = data;
        initializeUMAP();
        processGeneExpressionData();
    })
    .catch(error => console.error('Error loading UMAP data:', error));

// Process gene expression data from the loaded JSON
function processGeneExpressionData() {
    if (umapData.length > 0 && umapData[0].gene_expression) {
        const geneNameMapping = {};
        const rawGeneNames = Object.keys(umapData[0].gene_expression);

        // Process all gene names - trim whitespace and standardize case
        allGenes = rawGeneNames.map(raw => {
            const clean = raw.trim().toUpperCase();
            geneNameMapping[clean] = raw;  // Map clean name back to original
            return clean;
        });

        // Create gene expression matrix using cleaned names
        allGenes.forEach(clean => {
            const original = geneNameMapping[clean];
            geneExpressionData[clean] = umapData.map(cell =>
                (cell.gene_expression && original in cell.gene_expression)
                    ? cell.gene_expression[original]
                    : 0
            );
        });

        updateGeneSuggestions();
    }
}

// Update gene suggestions based on available genes
function updateGeneSuggestions() {
    const container = document.querySelector('.search-suggestions');
    if (container && allGenes.length > 0) {
        const suggestions = allGenes.slice(0, 4);
        container.innerHTML = `
            <span>Try:</span>
            ${suggestions.map(g => `<a href="#" class="gene-suggestion">${g}</a>`).join('')}
        `;

        document.querySelectorAll('.gene-suggestion').forEach(link => {
            link.addEventListener('click', function (e) {
                e.preventDefault();
                document.getElementById('geneSearch').value = this.textContent;
                searchGene();
            });
        });
    }
}

// Initialize UMAP plot with real data
function initializeUMAP() {
    const trace = {
        x: umapData.map(d => d.x),
        y: umapData.map(d => d.y),
        mode: 'markers',
        type: 'scatter',
        marker: {
            size: 4,
            color: umapData.map(d => d.cluster),
            colorscale: 'Portland',
            opacity: 0.8,
            line: {
                width: 0.5,
                color: 'rgba(0,0,0,0.1)'
            }
        },
        text: umapData.map(d => `Cell ID: ${d.cell_id}<br>Cluster: ${d.cluster}`),
        hoverinfo: 'text'
    };

    const layout = {
        title: '',
        xaxis: { 
            title: 'UMAP 1',
            gridcolor: 'rgba(0,0,0,0.05)'
        },
        yaxis: { 
            title: 'UMAP 2',
            gridcolor: 'rgba(0,0,0,0.05)'
        },
        hovermode: 'closest',
        margin: { t: 30, r: 20, b: 50, l: 50 },
        height: 500,
        plot_bgcolor: 'rgba(0,0,0,0)',
        paper_bgcolor: 'rgba(0,0,0,0)',
        font: {
            family: 'Poppins, sans-serif',
            size: 12,
            color: '#263238'
        }
    };

    Plotly.newPlot('umapPlot', [trace], layout);

    // Set up event listeners for controls
    document.getElementById('colorScheme').addEventListener('change', updateUMAP);
    document.getElementById('pointSize').addEventListener('input', updateUMAP);
    document.getElementById('opacity').addEventListener('input', updateUMAP);
}

// Update UMAP visualization based on controls
function updateUMAP() {
    const colorBy = document.getElementById('colorScheme').value;
    const pointSize = parseInt(document.getElementById('pointSize').value);
    const opacity = parseFloat(document.getElementById('opacity').value);
    
    let colorValues;
    let colorscale = 'Portland';
    
    switch(colorBy) {
        case 'clusters':
            colorValues = umapData.map(d => d.cluster);
            break;
        case 'cellType':
            colorValues = umapData.map(d => d.cell_type || d.cluster);
            colorscale = 'Viridis';
            break;
        case 'nGenes':
            colorValues = umapData.map(d => d.n_genes || 0);
            colorscale = 'Plasma';
            break;
    }
    
    const update = {
        'marker.size': pointSize,
        'marker.color': [colorValues],
        'marker.opacity': opacity,
        'marker.colorscale': colorscale
    };
    
    Plotly.update('umapPlot', update);
}

// Gene search functionality
document.getElementById('searchBtn').addEventListener('click', searchGene);
document.getElementById('geneSearch').addEventListener('keypress', function(e) {
    if (e.key === 'Enter') {
        searchGene();
    }
});

function searchGene() {
    const rawInput = document.getElementById('geneSearch').value;
    const gene = rawInput.trim().toUpperCase();
    const resultsDiv = document.getElementById('geneResults');

    if (!gene) {
        resultsDiv.innerHTML = `
            <div class="error-message">
                <i class="fas fa-exclamation-circle"></i>
                <p>Please enter a gene symbol</p>
            </div>
        `;
        return;
    }

    if (gene in geneExpressionData) {
        const expr = geneExpressionData[gene];
        const avg = (expr.reduce((a, b) => a + b, 0) / expr.length).toFixed(2);
        const max = Math.max(...expr).toFixed(2);
        const min = Math.min(...expr).toFixed(2);

        resultsDiv.innerHTML = `
            <div class="gene-header">
                <h3>${gene} Expression</h3>
                <div class="expression-summary">
                    <div class="summary-item"><span class="summary-value">${avg}</span><span class="summary-label">Average</span></div>
                    <div class="summary-item"><span class="summary-value">${max}</span><span class="summary-label">Max</span></div>
                    <div class="summary-item"><span class="summary-value">${min}</span><span class="summary-label">Min</span></div>
                </div>
            </div>
            <div class="expression-distribution">
                <h4>Expression Distribution</h4>
                <div class="distribution-chart" id="geneDistributionChart"></div>
            </div>
            <div class="cluster-expression">
                <h4>Cluster Expression</h4>
                <div class="cluster-bars"></div>
            </div>
        `;

        createDistributionPlot(gene);
        createClusterBars(gene);
        highlightGeneExpression(gene);
    } else {
        const suggestions = findSimilarGenes(gene).slice(0, 3).join(', ');
        resultsDiv.innerHTML = `
            <div class="error-message">
                <i class="fas fa-times-circle"></i>
                <p>Gene "${gene}" not found in dataset.</p>
                ${suggestions ? `<p class="similar-genes">Did you mean: ${suggestions}?</p>` : ''}
            </div>
        `;
    }
}

// Suggest similar genes if search fails
function findSimilarGenes(input) {
    const inputLower = input.toLowerCase();
    return allGenes
        .map(g => ({
            gene: g,
            score: (g.toLowerCase().startsWith(inputLower) ? 3 : 0) +
                   (g.toLowerCase().includes(inputLower) ? 2 : 0) -
                   Math.abs(g.length - input.length) * 0.1
        }))
        .sort((a, b) => b.score - a.score)
        .map(g => g.gene)
        .filter(g => g !== input);
}

function createDistributionPlot(geneName) {
    const expressionValues = geneExpressionData[geneName];
    const maxValue = Math.max(...expressionValues);
    const binSize = maxValue / 20;
    
    const bins = Array.from({length: 20}, (_, i) => i * binSize);
    const counts = Array(20).fill(0);
    
    expressionValues.forEach(value => {
        const binIndex = Math.min(Math.floor(value / binSize), 19);
        counts[binIndex]++;
    });
    
    const trace = {
        x: bins,
        y: counts,
        type: 'bar',
        marker: {
            color: 'rgba(94, 53, 177, 0.7)'
        }
    };
    
    const layout = {
        title: '',
        xaxis: { title: 'Expression Level' },
        yaxis: { title: 'Number of Cells' },
        margin: { t: 10, r: 20, b: 40, l: 50 },
        height: 200
    };
    
    Plotly.newPlot('geneDistributionChart', [trace], layout);
}

function createClusterBars(geneName) {
    const clusterExpression = {};
    const clusterCounts = {};
    
    umapData.forEach((cell, index) => {
        const cluster = cell.cluster;
        const expression = geneExpressionData[geneName][index];
        
        if (!clusterExpression[cluster]) {
            clusterExpression[cluster] = 0;
            clusterCounts[cluster] = 0;
        }
        
        clusterExpression[cluster] += expression;
        clusterCounts[cluster]++;
    });
    
    const clusters = Object.keys(clusterExpression).sort();
    const avgExpression = clusters.map(cluster => 
        (clusterExpression[cluster] / clusterCounts[cluster]).toFixed(2)
    );
    
    const maxExpression = Math.max(...avgExpression.map(Number));
    
    const container = document.querySelector('.cluster-bars');
    container.innerHTML = '';
    
    clusters.forEach((cluster, i) => {
        const value = avgExpression[i];
        const percentage = (value / maxExpression) * 100;
        
        const barContainer = document.createElement('div');
        barContainer.className = 'cluster-bar-container';
        
        const clusterLabel = document.createElement('span');
        clusterLabel.className = 'cluster-label';
        clusterLabel.textContent = `Cluster ${cluster}`;
        
        const barWrapper = document.createElement('div');
        barWrapper.className = 'bar-wrapper';
        
        const bar = document.createElement('div');
        bar.className = 'cluster-bar';
        bar.style.width = `${percentage}%`;
        
        const valueLabel = document.createElement('span');
        valueLabel.className = 'value-label';
        valueLabel.textContent = value;
        
        barWrapper.appendChild(bar);
        barContainer.appendChild(clusterLabel);
        barContainer.appendChild(barWrapper);
        barContainer.appendChild(valueLabel);
        container.appendChild(barContainer);
    });
}

function highlightGeneExpression(geneName) {
    const expressionValues = geneExpressionData[geneName];
    const maxExpression = Math.max(...expressionValues);
    
    const expressionTrace = {
        x: umapData.map(d => d.x),
        y: umapData.map(d => d.y),
        mode: 'markers',
        type: 'scatter',
        marker: {
            size: 6,
            color: expressionValues,
            colorscale: 'Viridis',
            cmin: 0,
            cmax: maxExpression,
            opacity: 0.8,
            line: {
                width: 0.5,
                color: 'rgba(0,0,0,0.1)'
            },
            colorbar: {
                title: 'Expression',
                thickness: 10
            }
        },
        text: umapData.map((d, i) => 
            `Cell ID: ${d.cell_id}<br>Cluster: ${d.cluster}<br>${geneName}: ${expressionValues[i].toFixed(2)}`
        ),
        hoverinfo: 'text',
        name: geneName
    };
    
    Plotly.newPlot('umapPlot', [expressionTrace], {
        title: '',
        xaxis: { title: 'UMAP 1' },
        yaxis: { title: 'UMAP 2' },
        hovermode: 'closest',
        margin: { t: 30, r: 20, b: 50, l: 50 },
        height: 500,
        plot_bgcolor: 'rgba(0,0,0,0)',
        paper_bgcolor: 'rgba(0,0,0,0)',
        font: {
            family: 'Poppins, sans-serif',
            size: 12,
            color: '#263238'
        }
    });
}

// Smooth scrolling for navigation
document.querySelectorAll('a[href^="#"]').forEach(anchor => {
    anchor.addEventListener('click', function(e) {
        e.preventDefault();
        
        const targetId = this.getAttribute('href');
        const targetElement = document.querySelector(targetId);
        
        if (targetElement) {
            window.scrollTo({
                top: targetElement.offsetTop - 80,
                behavior: 'smooth'
            });
            
            document.querySelectorAll('.nav-link').forEach(link => {
                link.classList.remove('active');
            });
            this.classList.add('active');
        }
    });
});

// Mobile navigation toggle
document.querySelector('.nav-toggle').addEventListener('click', function() {
    document.querySelector('.nav-links').classList.toggle('show');
});

// Highlight active section in navigation
window.addEventListener('scroll', function() {
    const scrollPosition = window.scrollY;
    
    document.querySelectorAll('section').forEach(section => {
        const sectionTop = section.offsetTop - 100;
        const sectionBottom = sectionTop + section.offsetHeight;
        
        if (scrollPosition >= sectionTop && scrollPosition < sectionBottom) {
            const sectionId = section.getAttribute('id');
            const correspondingNavLink = document.querySelector(`.nav-link[href="#${sectionId}"]`);
            
            if (correspondingNavLink) {
                document.querySelectorAll('.nav-link').forEach(link => {
                    link.classList.remove('active');
                });
                correspondingNavLink.classList.add('active');
            }
        }
    });
});

// Initialize when DOM is loaded
document.addEventListener('DOMContentLoaded', function() {
    // Initialization code if needed
});