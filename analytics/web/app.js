/**
 * NucleiTaxa Analytics Dashboard - Client
 * 
 * Real-time visualization of amplicon analysis pipeline
 * Connects to C++ analytics server via WebSocket
 */

class AnalyticsDashboard {
    constructor() {
        this.serverUrl = `ws://${window.location.hostname}:8888`;
        this.socket = null;
        this.metrics = {
            asv_count: 0,
            chimera_count: 0,
            mean_quality: 0,
            gc_percent: 0,
            taxonomy_distribution: {},
            elapsed_sec: 0,
            active_stage: 0
        };
        this.logEntries = [];
        this.initCharts();
        this.connect();
    }

    connect() {
        console.log('Connecting to analytics server...');
        
        try {
            this.socket = new WebSocket(this.serverUrl);
            
            this.socket.onopen = () => {
                this.setStatus('connected');
                this.addLog('Connected to analytics server', 'success');
            };
            
            this.socket.onmessage = (event) => {
                this.handleMetrics(JSON.parse(event.data));
            };
            
            this.socket.onerror = (error) => {
                this.addLog(`WebSocket error: ${error.message}`, 'error');
                this.setStatus('disconnected');
            };
            
            this.socket.onclose = () => {
                this.setStatus('disconnected');
                this.addLog('Disconnected from server. Reconnecting in 3s...', 'warning');
                setTimeout(() => this.connect(), 3000);
            };
        } catch (error) {
            console.error('Connection failed:', error);
            this.addLog(`Connection failed: ${error.message}`, 'error');
        }
    }

    handleMetrics(data) {
        this.metrics = data;
        this.updateUI();
    }

    updateUI() {
        // Update metric cards
        document.getElementById('asv-count').textContent = 
            this.metrics.asv_count.toLocaleString();
        document.getElementById('mean-quality').textContent = 
            this.metrics.mean_quality.toFixed(1);
        document.getElementById('gc-percent').textContent = 
            (this.metrics.gc_percent * 100).toFixed(1) + '%';
        document.getElementById('chimera-count').textContent = 
            this.metrics.chimera_count.toLocaleString();

        // Update elapsed time
        const mins = Math.floor(this.metrics.elapsed_sec / 60);
        const secs = Math.floor(this.metrics.elapsed_sec % 60);
        document.getElementById('elapsed-time').textContent = 
            `${String(mins).padStart(2, '0')}:${String(secs).padStart(2, '0')}`;

        // Update pipeline stage
        const stages = ['Waiting', 'QC & Filter', 'DADA2 Denoise', 'Chimera Remove', 
                       'Taxonomy RDP', 'Phylo FastTree', 'Krona Viz'];
        document.getElementById('pipeline-stage').textContent = 
            stages[this.metrics.active_stage] || 'Complete';

        // Update charts
        this.updateTaxaPie();
        this.updateQualityHistogram();
        this.updateStageProgress();
    }

    updateTaxaPie() {
        const taxa = Object.entries(this.metrics.taxonomy_distribution);
        
        const labels = taxa.map(([name]) => {
            // Truncate long names
            return name.length > 20 ? name.substring(0, 17) + '...' : name;
        });
        
        const values = taxa.map(([_, count]) => count);

        // Filter out zero values
        const filtered = taxa.filter(([_, count]) => count > 0);
        const labels_filtered = filtered.map(([name]) => {
            return name.length > 20 ? name.substring(0, 17) + '...' : name;
        });
        const values_filtered = filtered.map(([_, count]) => count);

        if (values_filtered.length === 0) return;

        const data = [{
            labels: labels_filtered,
            values: values_filtered,
            type: 'pie',
            marker: {
                colors: this.generateColors(labels_filtered.length)
            },
            hovertemplate: '<b>%{label}</b><br>Count: %{value}<extra></extra>'
        }];

        const layout = {
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            font: { color: '#e2e8f0', family: 'Arial' },
            margin: { l: 0, r: 0, t: 0, b: 0 },
            showlegend: true,
            legend: {
                x: 1.1,
                y: 1,
                bgcolor: 'rgba(30, 41, 59, 0.8)',
                bordercolor: '#334155',
                borderwidth: 1
            }
        };

        Plotly.newPlot('taxa-pie', data, layout, { responsive: true });
    }

    updateQualityHistogram() {
        // Simulated quality score distribution
        // In production, read from actual QC metrics file
        const qualityScores = this.generateQualityDistribution();
        
        const data = [{
            x: qualityScores,
            type: 'histogram',
            nbinsx: 30,
            marker: { color: '#1e40af' },
            hovertemplate: 'Q-score: %{x}<br>Frequency: %{y}<extra></extra>'
        }];

        const layout = {
            paper_bgcolor: 'rgba(0,0,0,0)',
            plot_bgcolor: 'rgba(0,0,0,0)',
            font: { color: '#e2e8f0', family: 'Arial' },
            xaxis: { title: 'Phred Score (Q)' },
            yaxis: { title: 'Frequency' },
            margin: { l: 50, r: 20, t: 20, b: 40 }
        };

        Plotly.newPlot('quality-hist', data, layout, { responsive: true });
    }

    updateStageProgress() {
        const stageItems = document.querySelectorAll('.stage-item');
        stageItems.forEach((item, idx) => {
            item.classList.remove('active', 'completed');
            if (idx < this.metrics.active_stage) {
                item.classList.add('completed');
            } else if (idx === this.metrics.active_stage) {
                item.classList.add('active');
            }
        });
    }

    initCharts() {
        // Initialize empty charts
        Plotly.newPlot('taxa-pie', [], {}, { responsive: true });
        Plotly.newPlot('quality-hist', [], {}, { responsive: true });
    }

    setStatus(status) {
        const dot = document.getElementById('connection-status');
        const text = document.getElementById('connection-text');
        
        if (status === 'connected') {
            dot.classList.remove('disconnected');
            dot.classList.add('connected');
            text.textContent = 'Connected';
        } else {
            dot.classList.remove('connected');
            dot.classList.add('disconnected');
            text.textContent = 'Disconnected';
        }
    }

    addLog(message, level = 'info') {
        const logContainer = document.getElementById('live-log');
        const entry = document.createElement('p');
        entry.className = `log-entry ${level}`;
        
        const timestamp = new Date().toLocaleTimeString();
        entry.textContent = `[${timestamp}] ${message}`;
        
        logContainer.appendChild(entry);
        logContainer.scrollTop = logContainer.scrollHeight;
        
        // Keep only last 50 entries
        while (logContainer.children.length > 50) {
            logContainer.removeChild(logContainer.firstChild);
        }
        
        this.logEntries.push({ timestamp, message, level });
    }

    generateColors(count) {
        const colors = [
            '#1e40af', '#16a34a', '#ea580c', '#dc2626', '#7c3aed',
            '#0891b2', '#d97706', '#059669', '#2563eb', '#9333ea'
        ];
        
        return Array(count).fill().map((_, i) => colors[i % colors.length]);
    }

    generateQualityDistribution() {
        // Generate realistic Q-score distribution
        // Most reads have high Q, tail of low Q
        const scores = [];
        for (let i = 0; i < 1000; i++) {
            // Beta distribution biased towards high quality
            const u = Math.random();
            const q = Math.floor(
                (Math.pow(u, 0.3) * 40)  // Skew towards high quality
            );
            scores.push(Math.max(0, Math.min(40, q)));
        }
        return scores;
    }
}

// Initialize dashboard on page load
document.addEventListener('DOMContentLoaded', () => {
    window.dashboard = new AnalyticsDashboard();
    
    // Log initial message
    window.dashboard.addLog('NucleiTaxa Analytics Dashboard v1.0', 'info');
    window.dashboard.addLog('Waiting for connection to analytics server...', 'info');
});

// Auto-refresh metrics every 500ms (if manually polling)
setInterval(() => {
    // Optional: Implement polling fallback if WebSocket unavailable
}, 500);
