/**
 * Heesch Polyform Explorer
 *
 * JavaScript application for viewing polyform renderings via Modal endpoints.
 */

// Configuration - Update these URLs after deploying to Modal
const CONFIG = {
    // Base URL for Modal endpoints - update after deployment
    baseUrl: 'https://YOUR_MODAL_USERNAME--heesch-renderings',

    // Endpoint paths
    endpoints: {
        render: '/render',
        renderSync: '/render_sync',
        list: '/list_polyforms',
        gridTypes: '/grid_types'
    },

    // Polling settings for async rendering
    pollInterval: 1000,  // ms
    maxPollAttempts: 30
};

// Sample polyforms for quick testing
const SAMPLE_POLYFORMS = [
    { name: 'Polyhex L', gridType: 'H', coords: '0,0_1,0_2,0_2,1' },
    { name: 'Polyhex Bar', gridType: 'H', coords: '-2,2_-1,1_0,0_1,0_2,0_2,1' },
    { name: 'Polyomino T', gridType: 'O', coords: '0,0_1,0_2,0_1,1' },
    { name: 'Polyomino L', gridType: 'O', coords: '0,0_0,1_0,2_1,2' },
    { name: 'Polyomino Square', gridType: 'O', coords: '0,0_0,1_1,0_1,1' },
    { name: 'Polyiamond Triangle', gridType: 'I', coords: '0,0_1,0_3,0' },
    { name: 'Polykite Trio', gridType: 'K', coords: '0,0_1,0_2,0' },
];

// Grid type info cache
let gridTypesCache = null;

/**
 * Initialize the application
 */
async function init() {
    console.log('Initializing Heesch Polyform Explorer...');

    // Set up event listeners
    document.getElementById('render-btn').addEventListener('click', handleRenderClick);
    document.getElementById('refresh-list-btn').addEventListener('click', loadPolyformList);
    document.getElementById('coords-input').addEventListener('keypress', (e) => {
        if (e.key === 'Enter') handleRenderClick();
    });

    // Load grid types
    await loadGridTypes();

    // Load available polyforms
    await loadPolyformList();

    // Populate sample polyforms
    populateSamples();
}

/**
 * Load grid types from API
 */
async function loadGridTypes() {
    const select = document.getElementById('grid-type');

    try {
        const response = await fetch(`${CONFIG.baseUrl}${CONFIG.endpoints.gridTypes}`);

        if (!response.ok) {
            throw new Error(`HTTP ${response.status}`);
        }

        const data = await response.json();
        gridTypesCache = data.grid_types;

        select.innerHTML = '';
        for (const gt of data.grid_types) {
            const option = document.createElement('option');
            option.value = gt.abbrev;
            option.textContent = `${gt.abbrev} - ${gt.full_name}`;
            select.appendChild(option);
        }
    } catch (error) {
        console.warn('Failed to load grid types from API, using defaults:', error);

        // Fallback to hardcoded grid types
        const defaultTypes = [
            { abbrev: 'O', name: 'Polyomino' },
            { abbrev: 'H', name: 'Polyhex' },
            { abbrev: 'I', name: 'Polyiamond' },
            { abbrev: 'o', name: 'Poly-[4.8.8]' },
            { abbrev: 'T', name: 'Poly-[3.6.3.6]' },
            { abbrev: 'A', name: 'Polyabolo' },
            { abbrev: 'D', name: 'Polydrafter' },
            { abbrev: 'K', name: 'Polykite' },
            { abbrev: 'h', name: 'Polyhalfcairo' },
            { abbrev: 'B', name: 'Polybevelhex' },
        ];

        select.innerHTML = '';
        for (const gt of defaultTypes) {
            const option = document.createElement('option');
            option.value = gt.abbrev;
            option.textContent = `${gt.abbrev} - ${gt.name}`;
            select.appendChild(option);
        }
    }
}

/**
 * Load list of available polyforms
 */
async function loadPolyformList() {
    const listContainer = document.getElementById('polyform-list');
    listContainer.innerHTML = '<p class="placeholder">Loading...</p>';

    try {
        const gridType = document.getElementById('grid-type').value;
        const url = gridType
            ? `${CONFIG.baseUrl}${CONFIG.endpoints.list}?grid_type=${gridType}`
            : `${CONFIG.baseUrl}${CONFIG.endpoints.list}`;

        const response = await fetch(url);

        if (!response.ok) {
            throw new Error(`HTTP ${response.status}`);
        }

        const data = await response.json();

        if (data.polyforms && data.polyforms.length > 0) {
            listContainer.innerHTML = '';
            for (const pf of data.polyforms) {
                const item = document.createElement('div');
                item.className = 'polyform-item';
                item.innerHTML = `
                    <span class="grid-type">${pf.grid_type}</span>
                    <span class="coords">${pf.coords}</span>
                `;
                item.addEventListener('click', () => {
                    document.getElementById('grid-type').value = pf.grid_type;
                    document.getElementById('coords-input').value = pf.coords;
                    handleRenderClick();
                });
                listContainer.appendChild(item);
            }
        } else {
            listContainer.innerHTML = '<p class="placeholder">No polyforms available yet. Render some to add them!</p>';
        }
    } catch (error) {
        console.warn('Failed to load polyform list:', error);
        listContainer.innerHTML = '<p class="placeholder">Could not load polyforms. API may be unavailable.</p>';
    }
}

/**
 * Populate sample polyforms
 */
function populateSamples() {
    const container = document.getElementById('samples');
    container.innerHTML = '';

    for (const sample of SAMPLE_POLYFORMS) {
        const item = document.createElement('div');
        item.className = 'sample-item';
        item.innerHTML = `
            <div class="sample-name">${sample.name}</div>
            <div class="sample-coords">${sample.gridType}: ${sample.coords}</div>
        `;
        item.addEventListener('click', () => {
            document.getElementById('grid-type').value = sample.gridType;
            document.getElementById('coords-input').value = sample.coords;
            handleRenderClick();
        });
        container.appendChild(item);
    }
}

/**
 * Handle render button click
 */
async function handleRenderClick() {
    const gridType = document.getElementById('grid-type').value;
    const coords = document.getElementById('coords-input').value.trim();

    if (!coords) {
        showStatus('Please enter coordinates', 'error');
        return;
    }

    await renderPolyform(gridType, coords);
}

/**
 * Render a polyform
 */
async function renderPolyform(gridType, coords) {
    const outputContainer = document.getElementById('render-output');
    const statusContainer = document.getElementById('render-status');

    outputContainer.innerHTML = '<div class="loading"></div>';
    showStatus('Requesting rendering...', 'computing');

    try {
        // First try the async endpoint
        const response = await fetch(
            `${CONFIG.baseUrl}${CONFIG.endpoints.render}?grid_type=${encodeURIComponent(gridType)}&coords=${encodeURIComponent(coords)}`
        );

        if (!response.ok) {
            throw new Error(`HTTP ${response.status}`);
        }

        const data = await response.json();

        if (data.status === 'available') {
            displayRendering(data);
        } else if (data.status === 'computing') {
            showStatus('Computing rendering... This may take a moment.', 'computing');
            // Poll for result
            await pollForRendering(gridType, data.coords);
        } else if (data.status === 'error') {
            showStatus(data.message, 'error');
            outputContainer.innerHTML = `<p class="placeholder">${data.message}</p>`;
        }
    } catch (error) {
        console.error('Render error:', error);

        // Fallback: Try to render locally using client-side SVG generation
        try {
            const svg = renderLocalSVG(gridType, coords);
            outputContainer.innerHTML = svg;
            showStatus('Rendered locally (API unavailable)', 'available');
        } catch (localError) {
            showStatus(`Failed to render: ${error.message}`, 'error');
            outputContainer.innerHTML = '<p class="placeholder">Rendering failed. Please check the API connection.</p>';
        }
    }
}

/**
 * Poll for async rendering result
 */
async function pollForRendering(gridType, coords, attempts = 0) {
    if (attempts >= CONFIG.maxPollAttempts) {
        showStatus('Rendering is taking longer than expected. Please try again later.', 'error');
        return;
    }

    await new Promise(resolve => setTimeout(resolve, CONFIG.pollInterval));

    try {
        const response = await fetch(
            `${CONFIG.baseUrl}${CONFIG.endpoints.render}?grid_type=${encodeURIComponent(gridType)}&coords=${encodeURIComponent(coords)}`
        );

        const data = await response.json();

        if (data.status === 'available') {
            displayRendering(data);
        } else if (data.status === 'computing') {
            showStatus(`Computing rendering... (attempt ${attempts + 1}/${CONFIG.maxPollAttempts})`, 'computing');
            await pollForRendering(gridType, coords, attempts + 1);
        } else {
            showStatus(data.message || 'Unknown error', 'error');
        }
    } catch (error) {
        showStatus(`Poll failed: ${error.message}`, 'error');
    }
}

/**
 * Display a rendering result
 */
function displayRendering(data) {
    const outputContainer = document.getElementById('render-output');
    outputContainer.innerHTML = data.svg;
    showStatus(`Rendered ${data.grid_name} (${data.coords})`, 'available');

    // Refresh the list to show newly computed polyforms
    loadPolyformList();
}

/**
 * Show a status message
 */
function showStatus(message, type) {
    const statusContainer = document.getElementById('render-status');
    statusContainer.textContent = message;
    statusContainer.className = `status-message ${type}`;
}

/**
 * Client-side SVG rendering fallback
 * This provides basic rendering when the API is unavailable
 */
function renderLocalSVG(gridType, coordsStr) {
    const coords = coordsStr.split('_').map(pair => {
        const [x, y] = pair.split(',').map(Number);
        return { x, y };
    });

    if (gridType === 'O') {
        return renderOminoSVG(coords);
    } else if (gridType === 'H') {
        return renderHexSVG(coords);
    } else {
        // Fallback to square rendering
        return renderOminoSVG(coords);
    }
}

/**
 * Render polyomino SVG locally
 */
function renderOminoSVG(coords) {
    const cellSize = 30;
    const padding = 20;

    const xs = coords.map(c => c.x);
    const ys = coords.map(c => c.y);
    const minX = Math.min(...xs);
    const minY = Math.min(...ys);
    const maxX = Math.max(...xs);
    const maxY = Math.max(...ys);

    const width = (maxX - minX + 1) * cellSize + 2 * padding;
    const height = (maxY - minY + 1) * cellSize + 2 * padding;

    let cells = '';
    for (const {x, y} of coords) {
        const px = (x - minX) * cellSize + padding;
        const py = (y - minY) * cellSize + padding;
        cells += `<rect x="${px}" y="${py}" width="${cellSize}" height="${cellSize}" fill="#FFD700" stroke="#000" stroke-width="2"/>`;
    }

    return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" width="${width}" height="${height}">${cells}</svg>`;
}

/**
 * Render polyhex SVG locally
 */
function renderHexSVG(coords) {
    const hexSize = 20;
    const padding = 40;

    function hexToPixel(x, y) {
        return {
            px: hexSize * 1.5 * x,
            py: hexSize * Math.sqrt(3) * (y + x * 0.5)
        };
    }

    function hexVertices(cx, cy) {
        const vertices = [];
        for (let i = 0; i < 6; i++) {
            const angle = Math.PI / 3 * i;
            vertices.push({
                x: cx + hexSize * Math.cos(angle),
                y: cy + hexSize * Math.sin(angle)
            });
        }
        return vertices;
    }

    // Calculate bounds
    let allVertices = [];
    for (const {x, y} of coords) {
        const {px, py} = hexToPixel(x, y);
        allVertices = allVertices.concat(hexVertices(px, py));
    }

    const vxs = allVertices.map(v => v.x);
    const vys = allVertices.map(v => v.y);
    const minX = Math.min(...vxs);
    const minY = Math.min(...vys);
    const maxX = Math.max(...vxs);
    const maxY = Math.max(...vys);

    const width = maxX - minX + 2 * padding;
    const height = maxY - minY + 2 * padding;
    const offsetX = padding - minX;
    const offsetY = padding - minY;

    let hexagons = '';
    for (const {x, y} of coords) {
        const {px, py} = hexToPixel(x, y);
        const vertices = hexVertices(px + offsetX, py + offsetY);
        const points = vertices.map(v => `${v.x},${v.y}`).join(' ');
        hexagons += `<polygon points="${points}" fill="#FFD700" stroke="#000" stroke-width="2"/>`;
    }

    return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" width="${width}" height="${height}">${hexagons}</svg>`;
}

// Initialize when DOM is ready
document.addEventListener('DOMContentLoaded', init);
