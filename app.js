/**
 * WebLucy - Cheminformatics Teaching Tool
 * Bond insertion algorithm adapted from LUCY (Steinbeck, 1995)
 */

// Element data: valence, typical max bond order
const ELEMENT_DATA = {
    'C': { valence: 4, maxBond: 3 },
    'N': { valence: 3, maxBond: 3 },
    'O': { valence: 2, maxBond: 2 },
    'S': { valence: 2, maxBond: 2 },
    'P': { valence: 3, maxBond: 3 },
    'F': { valence: 1, maxBond: 1 },
    'Cl': { valence: 1, maxBond: 1 },
    'Br': { valence: 1, maxBond: 1 },
    'I': { valence: 1, maxBond: 1 },
    'Si': { valence: 4, maxBond: 2 },
    'B': { valence: 3, maxBond: 3 }
};

// Global state
let atoms = [];           // Array of atom objects {symbol, valence, hCount, maxBond, freeValence}
let N = 0;                // Number of heavy atoms
let matrix = [];          // Connection matrix (N x N)
let bondSum = [];         // Current bond sum for each atom
let states = [];          // History of states for backward navigation
let currentStateIndex = -1;
let autoStepInterval = null;

// API endpoint for 2D depiction
const DEPICT_API_URL = 'https://api.naturalproducts.net/latest/depict/2D';

/**
 * Parse molecular formula into atom counts
 * e.g., "C5H12" -> {C: 5, H: 12}
 */
function parseFormula(formula) {
    const counts = {};
    // Match element symbols (1-2 chars) followed by optional number
    const regex = /([A-Z][a-z]?)(\d*)/g;
    let match;

    while ((match = regex.exec(formula)) !== null) {
        if (match[1]) {
            const element = match[1];
            const count = match[2] ? parseInt(match[2]) : 1;
            counts[element] = (counts[element] || 0) + count;
        }
    }

    return counts;
}

/**
 * Calculate and distribute hydrogens among heavy atoms
 * Each atom gets a random valid H count, sum must equal totalH
 * Max H per atom = valence - 1 (to ensure at least 1 bond to other heavy atoms)
 */
function distributeHydrogens(heavyAtoms, totalH) {
    // Calculate max possible H (each atom keeps at least 1 valence for heavy-atom bonds)
    let maxPossibleH = 0;
    for (const atom of heavyAtoms) {
        atom.maxH = atom.valence - 1; // Reserve 1 valence for bonds to other heavy atoms
        maxPossibleH += atom.maxH;
    }

    // Check if formula is valid
    if (totalH > maxPossibleH) {
        throw new Error(`Too many hydrogens (${totalH}) for the given heavy atoms (max: ${maxPossibleH})`);
    }

    // Calculate minimum H needed (if all atoms use max bonds to each other)
    // This is a simplification - we allow 0 H per atom as minimum
    const minPossibleH = 0;

    if (totalH < minPossibleH) {
        throw new Error(`Too few hydrogens (${totalH}) for a valid structure`);
    }

    // Distribute H randomly but validly
    let remainingH = totalH;
    const atomIndices = heavyAtoms.map((_, i) => i);

    // First pass: assign random H to each atom (0 to maxH)
    for (const atom of heavyAtoms) {
        atom.hCount = 0;
    }

    // Randomly distribute hydrogens
    while (remainingH > 0) {
        // Shuffle to randomize distribution
        const availableAtoms = heavyAtoms.filter(a => a.hCount < a.maxH);
        if (availableAtoms.length === 0) break;

        const randomAtom = availableAtoms[Math.floor(Math.random() * availableAtoms.length)];
        randomAtom.hCount++;
        remainingH--;
    }

    // Calculate free valence for each atom
    for (const atom of heavyAtoms) {
        atom.freeValence = atom.valence - atom.hCount;
    }

    // Calculate bonds needed between heavy atoms
    let totalFreeValence = 0;
    for (const atom of heavyAtoms) {
        totalFreeValence += atom.freeValence;
    }
    const bondsNeeded = totalFreeValence / 2;

    return bondsNeeded;
}

/**
 * Initialize the molecule from a formula
 */
function initializeMolecule() {
    const formulaInput = document.getElementById('formula').value.trim();

    if (!formulaInput) {
        showStatus('Please enter a molecular formula', 'error');
        return;
    }

    try {
        const counts = parseFormula(formulaInput);

        // Extract hydrogens
        const totalH = counts['H'] || 0;
        delete counts['H'];

        // Build atoms array from heavy atoms
        atoms = [];
        let atomIndex = 1;

        for (const [element, count] of Object.entries(counts)) {
            if (!ELEMENT_DATA[element]) {
                throw new Error(`Unknown element: ${element}`);
            }

            for (let i = 0; i < count; i++) {
                atoms.push({
                    index: atomIndex++,
                    symbol: element,
                    valence: ELEMENT_DATA[element].valence,
                    maxBond: ELEMENT_DATA[element].maxBond,
                    hCount: 0,
                    freeValence: ELEMENT_DATA[element].valence
                });
            }
        }

        N = atoms.length;

        if (N < 2) {
            throw new Error('Need at least 2 heavy atoms to create a connection matrix');
        }

        // Store total H for validation
        atoms.totalH = totalH;

        // Distribute hydrogens and calculate bonds needed
        const bondsNeeded = distributeHydrogens(atoms, totalH);
        atoms.bondsNeeded = bondsNeeded;

        // Initialize matrix and bond sums
        matrix = [];
        bondSum = [];
        for (let i = 0; i <= N; i++) {
            matrix[i] = [];
            bondSum[i] = 0;
            for (let j = 0; j <= N; j++) {
                matrix[i][j] = 0;
            }
        }

        // Reset state
        states = [];
        currentStateIndex = -1;

        // Save initial state
        saveState(0, 0, 'Initial empty matrix');

        // Update UI
        renderMatrix();
        renderAtomInfo();
        updateControls();
        showStatus(`Created ${N}x${N} matrix for formula ${formulaInput}. Need ${bondsNeeded} bonds between heavy atoms.`, 'success');
        document.getElementById('step-info').textContent = `Step 0 / ? - Click Forward to begin bond insertion`;

    } catch (error) {
        showStatus(error.message, 'error');
    }
}

/**
 * Save current state for history
 */
function saveState(i, j, description) {
    // Deep copy matrix
    const matrixCopy = matrix.map(row => [...row]);
    const bondSumCopy = [...bondSum];

    // Remove any states after current index (if we went backward and now forward)
    states = states.slice(0, currentStateIndex + 1);

    states.push({
        matrix: matrixCopy,
        bondSum: bondSumCopy,
        i: i,
        j: j,
        description: description
    });

    currentStateIndex = states.length - 1;
}

/**
 * Restore a state from history
 */
function restoreState(index) {
    if (index < 0 || index >= states.length) return;

    const state = states[index];
    matrix = state.matrix.map(row => [...row]);
    bondSum = [...state.bondSum];
    currentStateIndex = index;

    renderMatrix(state.i, state.j);
    renderAtomInfo();
    updateControls();

    document.getElementById('step-info').textContent =
        `Step ${currentStateIndex} / ${states.length - 1} - ${state.description}`;
}

/**
 * Get next position in matrix traversal (column by column, lower triangle)
 * Traversal order: (2,1), (3,1), (4,1)... (N,1), (3,2), (4,2)... (N,2), etc.
 */
function getNextPosition(i, j) {
    if (i < N) {
        return { i: i + 1, j: j };
    } else {
        // Move to next column
        if (j + 1 < N) {
            return { i: j + 2, j: j + 1 };
        } else {
            return null; // End of matrix
        }
    }
}

/**
 * Get atom data (1-indexed like in original C code)
 */
function getAtom(index) {
    return atoms[index - 1];
}

/**
 * Try to set a bond at position (i,j) with given value
 * Returns true if successful, false if violates constraints
 */
function trySetBond(i, j, bondValue) {
    const atomI = getAtom(i);
    const atomJ = getAtom(j);

    // Check max bond order constraint
    if (bondValue > Math.min(atomI.maxBond, atomJ.maxBond)) {
        return false;
    }

    // Check free valence constraints
    const currentBond = matrix[i][j];
    const newBondSumI = bondSum[i] - currentBond + bondValue;
    const newBondSumJ = bondSum[j] - currentBond + bondValue;

    if (newBondSumI > atomI.freeValence || newBondSumJ > atomJ.freeValence) {
        return false;
    }

    // Apply the bond
    bondSum[i] = newBondSumI;
    bondSum[j] = newBondSumJ;
    matrix[i][j] = bondValue;
    matrix[j][i] = bondValue;

    return true;
}

/**
 * Check if we're at the last position in the matrix
 */
function isAtEnd(i, j) {
    return i === N && j === N - 1;
}

/**
 * Step forward: depth-first traversal like mol.c
 * 1. If not at end, move to next cell with bond=0
 * 2. If at end, backtrack and increment
 */
function stepForward() {
    if (states.length === 0) return;

    const currentState = states[currentStateIndex];
    let i = currentState.i;
    let j = currentState.j;

    // If we're at the beginning, start at position (2,1) with bond=1
    if (i === 0 && j === 0) {
        i = 2;
        j = 1;
        trySetBond(i, j, 1);
        saveState(i, j, `Set bond (${i},${j}) = 1`);
        renderMatrix(i, j);
        renderAtomInfo();
        updateControls();
        document.getElementById('step-info').textContent =
            `Step ${currentStateIndex} - Set bond (${i},${j}) = 1`;
        return;
    }

    // Try to move to next position (depth-first: go deeper first)
    const next = getNextPosition(i, j);
    if (next !== null) {
        // Move to next cell, start with bond=1
        trySetBond(next.i, next.j, 1);
        saveState(next.i, next.j, `Set bond (${next.i},${next.j}) = 1`);
        renderMatrix(next.i, next.j);
        renderAtomInfo();
        updateControls();
        document.getElementById('step-info').textContent =
            `Step ${currentStateIndex} - Set bond (${next.i},${next.j}) = 1`;
        checkCompletion();
        return;
    }

    // At end of matrix - need to backtrack and try incrementing
    backtrackAndIncrement();
}

/**
 * Get next bond value to try in the sequence: 1, 2, 3, 0
 * Returns null if all values exhausted
 */
function getNextBondValue(currentBond) {
    if (currentBond === 0) {
        return null; // Exhausted all options (we tried 1, 2, 3, 0)
    } else if (currentBond < 3) {
        return currentBond + 1; // Try next higher: 1→2, 2→3
    } else {
        return 0; // After 3, try 0 (no bond)
    }
}

/**
 * Backtrack through states and find a position where we can try the next bond value
 */
function backtrackAndIncrement() {
    // We need to find the most recent position where we can try a different bond
    // Start from current state and work backwards

    for (let stateIdx = currentStateIndex; stateIdx >= 1; stateIdx--) {
        const state = states[stateIdx];
        const i = state.i;
        const j = state.j;
        const currentBondAtState = state.matrix[i][j];

        // Get next bond value to try (sequence: 1, 2, 3, 0)
        let nextBond = getNextBondValue(currentBondAtState);

        while (nextBond !== null) {
            // Restore to the state BEFORE this bond was set
            const prevState = states[stateIdx - 1];
            matrix = prevState.matrix.map(row => [...row]);
            bondSum = [...prevState.bondSum];

            if (trySetBond(i, j, nextBond)) {
                // Successfully set new bond! Save this new state
                // Truncate history to this point
                states = states.slice(0, stateIdx);
                currentStateIndex = stateIdx - 1;

                saveState(i, j, `Set bond (${i},${j}) = ${nextBond}`);
                renderMatrix(i, j);
                renderAtomInfo();
                updateControls();
                document.getElementById('step-info').textContent =
                    `Step ${currentStateIndex} - Backtrack: Set bond (${i},${j}) = ${nextBond}`;
                checkCompletion();
                return;
            }

            // This bond value didn't work, try next in sequence
            nextBond = getNextBondValue(nextBond);
        }

        // Exhausted all bond values at this position, continue backtracking
    }

    // Exhausted all possibilities
    // Restore current state
    restoreState(currentStateIndex);
    showStatus('Explored all valid bond configurations for this H-distribution.', 'info');
}

/**
 * Step backward: go to previous state
 */
function stepBackward() {
    if (currentStateIndex > 0) {
        restoreState(currentStateIndex - 1);
    }
}

/**
 * Reset matrix to initial state
 */
function resetMatrix() {
    if (states.length > 0) {
        restoreState(0);
        hideMolecule();
        showStatus('Matrix reset to initial state', 'info');
    }
}

/**
 * Toggle auto-stepping
 */
function toggleAutoStep() {
    const btn = document.getElementById('btn-auto');

    if (autoStepInterval) {
        clearInterval(autoStepInterval);
        autoStepInterval = null;
        btn.textContent = 'Auto Step';
        btn.classList.remove('btn-primary');
        btn.classList.add('btn-secondary');
    } else {
        autoStepInterval = setInterval(() => {
            stepForward();
            // Stop if we can't go forward anymore
            const currentState = states[currentStateIndex];
            if (currentState.i === N && currentState.j === N - 1) {
                toggleAutoStep();
            }
        }, 500);
        btn.textContent = 'Stop';
        btn.classList.remove('btn-secondary');
        btn.classList.add('btn-primary');
    }
}

/**
 * Check if the molecule is connected (all atoms reachable via bonds)
 * Uses depth-first search starting from atom 1
 */
function isConnected() {
    if (N <= 1) return true;

    const visited = new Array(N + 1).fill(false);
    let visitedCount = 0;

    // Depth-first search
    function dfs(atomIndex) {
        if (visited[atomIndex]) return;
        visited[atomIndex] = true;
        visitedCount++;

        // Visit all bonded neighbors
        for (let j = 1; j <= N; j++) {
            if (matrix[atomIndex][j] > 0 || matrix[j][atomIndex] > 0) {
                if (!visited[j]) {
                    dfs(j);
                }
            }
        }
    }

    // Start DFS from atom 1
    dfs(1);

    // Molecule is connected if all atoms were visited
    return visitedCount === N;
}

/**
 * Check if we've reached a valid complete molecule
 */
function checkCompletion() {
    // Check if all atoms have their free valences exactly satisfied
    let totalBonds = 0;
    let allSatisfied = true;

    for (let k = 1; k <= N; k++) {
        const atom = getAtom(k);
        // Each atom's bondSum should equal its freeValence (valence - hCount)
        if (bondSum[k] !== atom.freeValence) {
            allSatisfied = false;
        }
        totalBonds += bondSum[k];
    }

    totalBonds /= 2; // Each bond counted twice

    if (allSatisfied && totalBonds === atoms.bondsNeeded) {
        // Check if molecule is connected
        if (isConnected()) {
            showStatus(`Valid connected structure found! ${Math.round(totalBonds)} bonds between heavy atoms.`, 'success');
            // Display the molecule
            displayMolecule();
        } else {
            showStatus(`Structure has correct valences but is not connected (disconnected fragments).`, 'warning');
        }
    } else {
        // Just update status, keep the last valid molecule displayed
        showStatus('Exploring bond configurations...', 'info');
    }
}

/**
 * Render the connection matrix
 */
function renderMatrix(highlightI = 0, highlightJ = 0) {
    if (N === 0) return;

    let html = '<table class="matrix">';

    // Header row
    html += '<tr><th></th>';
    for (let j = 1; j <= N; j++) {
        const atom = getAtom(j);
        html += `<th>${atom.symbol}${j}</th>`;
    }
    html += '</tr>';

    // Matrix rows
    for (let i = 1; i <= N; i++) {
        const atomI = getAtom(i);
        html += `<tr><th>${atomI.symbol}${i}</th>`;

        for (let j = 1; j <= N; j++) {
            let cellClass = '';
            let value = '';

            if (i === j) {
                cellClass = 'diagonal';
                value = '-';
            } else {
                const bond = matrix[i][j];
                cellClass = `bond-${bond}`;
                value = bond;

                if ((i === highlightI && j === highlightJ) ||
                    (i === highlightJ && j === highlightI)) {
                    cellClass += ' current';
                }
            }

            html += `<td class="${cellClass}">${value}</td>`;
        }
        html += '</tr>';
    }

    html += '</table>';
    document.getElementById('matrix-container').innerHTML = html;
}

/**
 * Handle H-Count change from input field
 */
function updateHCount(atomIndex, newValue) {
    const atom = getAtom(atomIndex);
    const oldValue = atom.hCount;
    const parsed = parseInt(newValue);

    // Validate
    if (isNaN(parsed) || parsed < 0 || parsed > atom.maxH) {
        renderAtomInfo(); // Reset to valid value
        return;
    }

    // Calculate what total H would be with this change
    let newTotalH = 0;
    for (let i = 1; i <= N; i++) {
        if (i === atomIndex) {
            newTotalH += parsed;
        } else {
            newTotalH += getAtom(i).hCount;
        }
    }

    // Check if total H matches required
    if (newTotalH !== atoms.totalH) {
        showStatus(`H-Count sum must equal ${atoms.totalH}. Current sum would be ${newTotalH}.`, 'warning');
        renderAtomInfo(); // Reset
        return;
    }

    // Update the atom
    atom.hCount = parsed;
    atom.freeValence = atom.valence - atom.hCount;

    // Recalculate bonds needed
    let totalFreeValence = 0;
    for (let i = 1; i <= N; i++) {
        totalFreeValence += getAtom(i).freeValence;
    }
    atoms.bondsNeeded = totalFreeValence / 2;

    // Reset matrix since H distribution changed
    resetMatrixState();

    showStatus(`Updated H-Count. Now need ${atoms.bondsNeeded} bonds between heavy atoms.`, 'success');
}

/**
 * Reset matrix state after H distribution change
 */
function resetMatrixState() {
    // Reset matrix and bond sums
    for (let i = 0; i <= N; i++) {
        bondSum[i] = 0;
        for (let j = 0; j <= N; j++) {
            matrix[i][j] = 0;
        }
    }

    // Reset state history
    states = [];
    currentStateIndex = -1;
    saveState(0, 0, 'Initial empty matrix');

    // Hide molecule panel
    hideMolecule();

    renderMatrix();
    renderAtomInfo();
    updateControls();
    document.getElementById('step-info').textContent = `Step 0 / ? - Click Forward to begin bond insertion`;
}

/**
 * Render atom information table
 */
function renderAtomInfo() {
    if (atoms.length === 0) return;

    let html = '<table>';
    html += '<tr><th>Atom</th><th>Symbol</th><th>Valence</th><th>H-Count</th><th>Bond Sum</th><th>Free Valence</th></tr>';

    for (let i = 1; i <= N; i++) {
        const atom = getAtom(i);
        const freeValence = atom.freeValence - bondSum[i];
        const status = freeValence < 0 ? ' style="color: red;"' : '';

        html += `<tr${status}>`;
        html += `<td>${i}</td>`;
        html += `<td>${atom.symbol}</td>`;
        html += `<td>${atom.valence}</td>`;
        html += `<td><input type="number" min="0" max="${atom.maxH}" value="${atom.hCount}"
                    onchange="updateHCount(${i}, this.value)"
                    style="width: 50px; text-align: center;"></td>`;
        html += `<td>${bondSum[i]}</td>`;
        html += `<td>${freeValence}</td>`;
        html += '</tr>';
    }

    html += '</table>';

    // Add summary
    let currentHSum = 0;
    for (let i = 1; i <= N; i++) {
        currentHSum += getAtom(i).hCount;
    }

    const hSumClass = currentHSum === atoms.totalH ? 'color: green;' : 'color: red;';
    html += `<p style="margin-top: 10px;">`;
    html += `<strong>H-Count sum:</strong> <span style="${hSumClass}">${currentHSum}</span> (required: ${atoms.totalH || 0})<br>`;
    html += `<strong>Bonds needed:</strong> ${atoms.bondsNeeded || 0}`;
    html += `</p>`;

    document.getElementById('atom-info').innerHTML = html;
}

/**
 * Update control button states
 */
function updateControls() {
    document.getElementById('btn-reset').disabled = states.length === 0;
    document.getElementById('btn-back').disabled = currentStateIndex <= 0;
    document.getElementById('btn-forward').disabled = states.length === 0;
    document.getElementById('btn-auto').disabled = states.length === 0;
    document.getElementById('btn-silent').disabled = N === 0;
}

/**
 * Show status message
 */
function showStatus(message, type = 'info') {
    const statusEl = document.getElementById('status');
    statusEl.textContent = message;
    statusEl.className = `status ${type}`;
}

/**
 * Create a Kekule molecule from the current connection matrix (with circular layout)
 */
function createMoleculeWithCircularLayout() {
    const mol = new Kekule.Molecule();

    // Create atoms and position them in a circle
    const kekuleAtoms = [];
    const radius = 1.5;
    const centerX = 0;
    const centerY = 0;

    for (let i = 1; i <= N; i++) {
        const atom = getAtom(i);
        const angle = (2 * Math.PI * (i - 1)) / N - Math.PI / 2; // Start from top
        const x = centerX + radius * Math.cos(angle);
        const y = centerY + radius * Math.sin(angle);

        const kAtom = new Kekule.Atom();
        kAtom.setSymbol(atom.symbol);
        kAtom.setCoord2D({ x: x, y: y });
        mol.appendNode(kAtom);
        kekuleAtoms.push(kAtom);
    }

    // Create bonds from connection matrix
    for (let i = 1; i <= N; i++) {
        for (let j = i + 1; j <= N; j++) {
            const bondOrder = matrix[i][j];
            if (bondOrder > 0) {
                const bond = new Kekule.Bond();
                bond.setBondOrder(bondOrder);
                bond.setConnectedObjs([kekuleAtoms[i - 1], kekuleAtoms[j - 1]]);
                mol.appendConnector(bond);
            }
        }
    }

    return mol;
}

/**
 * Create a Kekule molecule without coordinates (for 2D generation)
 */
function createMoleculeWithoutCoords() {
    const mol = new Kekule.Molecule();
    const kekuleAtoms = [];

    // Create atoms without coordinates
    for (let i = 1; i <= N; i++) {
        const atom = getAtom(i);
        const kAtom = new Kekule.Atom();
        kAtom.setSymbol(atom.symbol);
        mol.appendNode(kAtom);
        kekuleAtoms.push(kAtom);
    }

    // Create bonds from connection matrix
    for (let i = 1; i <= N; i++) {
        for (let j = i + 1; j <= N; j++) {
            const bondOrder = matrix[i][j];
            if (bondOrder > 0) {
                const bond = new Kekule.Bond();
                bond.setBondOrder(bondOrder);
                bond.setConnectedObjs([kekuleAtoms[i - 1], kekuleAtoms[j - 1]]);
                mol.appendConnector(bond);
            }
        }
    }

    return mol;
}

/**
 * Generate 2D coordinates using OpenBabel
 */
function generate2DCoordinates(mol, callback) {
    if (!openBabelReady || !Kekule.Calculator || !Kekule.Calculator.ObStructure2DGenerator) {
        console.log('OpenBabel 2D generator not available');
        callback(null);
        return;
    }

    try {
        const generator = new Kekule.Calculator.ObStructure2DGenerator();
        generator.setSourceMol(mol);

        // Use async execution
        generator.execute(function(generatedMol) {
            if (generatedMol) {
                console.log('2D coordinates generated successfully');
                callback(generator.getGeneratedMol());
            } else {
                console.warn('2D generation returned no molecule');
                callback(null);
            }
        }, function(error) {
            console.warn('2D generation failed:', error);
            callback(null);
        });
    } catch (e) {
        console.error('Error in 2D generation:', e);
        callback(null);
    }
}

/**
 * Create a Kekule molecule from the current connection matrix
 */
function createMoleculeFromMatrix() {
    try {
        // Always create with circular layout as fallback
        return createMoleculeWithCircularLayout();
    } catch (e) {
        console.error('Error creating molecule:', e);
        return null;
    }
}

/**
 * Convert molecule to SMILES using Kekule.js
 */
function moleculeToSmiles() {
    try {
        const mol = createMoleculeWithoutCoords();
        if (!mol) {
            console.error('Failed to create molecule for SMILES conversion');
            return null;
        }

        // Use Kekule.IO to save as SMILES
        const smiles = Kekule.IO.saveFormatData(mol, 'smi');
        console.log('Generated SMILES:', smiles);
        return smiles;
    } catch (e) {
        console.error('Error converting to SMILES:', e);
        return null;
    }
}

/**
 * Fetch 2D depiction SVG from Natural Products API
 */
async function fetch2DDepiction(smiles) {
    const params = new URLSearchParams({
        smiles: smiles,
        toolkit: 'rdkit',
        width: 300,
        height: 300
    });

    const url = `${DEPICT_API_URL}?${params.toString()}`;
    console.log('Fetching 2D depiction from:', url);

    try {
        const response = await fetch(url);
        if (!response.ok) {
            throw new Error(`API returned ${response.status}: ${response.statusText}`);
        }
        const svgText = await response.text();
        return svgText;
    } catch (e) {
        console.error('Error fetching 2D depiction:', e);
        return null;
    }
}

/**
 * Display the molecule in the viewer
 */
async function displayMolecule() {
    console.log('displayMolecule called');

    const panel = document.getElementById('molecule-panel');
    const viewerElement = document.getElementById('molecule-viewer');
    const statusElement = document.getElementById('molecule-status');

    if (!viewerElement) {
        console.error('molecule-viewer element not found');
        return;
    }

    // Check if Kekule is available for SMILES generation
    if (typeof Kekule === 'undefined') {
        console.error('Kekule.js is not loaded');
        viewerElement.innerHTML = '<p style="color: red; text-align: center; padding-top: 130px;">Kekule.js not loaded</p>';
        return;
    }

    // Show loading state
    viewerElement.innerHTML = '<p style="text-align: center; color: #666; padding-top: 130px;">Generating 2D structure...</p>';

    // Convert molecule to SMILES
    const smiles = moleculeToSmiles();
    if (!smiles) {
        viewerElement.innerHTML = '<p style="color: red; text-align: center; padding-top: 130px;">Failed to generate SMILES</p>';
        return;
    }

    // Fetch 2D depiction from API
    const svg = await fetch2DDepiction(smiles);

    if (svg) {
        // Display the SVG
        viewerElement.innerHTML = svg;

        // Make the SVG responsive
        const svgElement = viewerElement.querySelector('svg');
        if (svgElement) {
            svgElement.setAttribute('width', '100%');
            svgElement.setAttribute('height', '100%');
            svgElement.style.maxWidth = '100%';
            svgElement.style.maxHeight = '300px';
        }

        if (panel) {
            panel.classList.add('has-molecule');
        }

        if (statusElement) {
            statusElement.textContent = `Valid structure (SMILES: ${smiles})`;
            statusElement.style.color = '#28a745';
        }

        console.log('2D depiction displayed successfully');
    } else {
        viewerElement.innerHTML = '<p style="color: orange; text-align: center; padding-top: 130px;">Could not fetch 2D depiction</p>';

        if (statusElement) {
            statusElement.textContent = `Valid structure but depiction failed (SMILES: ${smiles})`;
            statusElement.style.color = '#ffc107';
        }
    }
}

/**
 * Clear the molecule from the viewer (only called on full reset)
 */
function hideMolecule() {
    const panel = document.getElementById('molecule-panel');
    const viewerElement = document.getElementById('molecule-viewer');
    const statusElement = document.getElementById('molecule-status');

    if (panel) {
        panel.classList.remove('has-molecule');
    }

    if (viewerElement) {
        viewerElement.innerHTML = '<p style="text-align: center; color: #999; padding-top: 130px;">Molecule will appear here</p>';
    }

    if (statusElement) {
        statusElement.textContent = 'Find a valid structure to see the 2D depiction.';
        statusElement.style.color = '#666';
    }
}

/**
 * Silent generation - enumerate all valid structures and download as zip
 * This enumerates ALL possible H-Count distributions, then for each
 * distribution enumerates all valid bond configurations.
 */
async function silentGeneration() {
    if (N === 0) {
        showStatus('Please create a matrix first', 'error');
        return;
    }

    const btn = document.getElementById('btn-silent');
    btn.disabled = true;
    btn.textContent = 'Generating...';
    showStatus('Generating all valid structures (enumerating all H distributions)...', 'info');

    // Use setTimeout to allow UI to update
    await new Promise(resolve => setTimeout(resolve, 100));

    const validSmiles = [];

    // Get atom info
    const atomInfo = [];
    for (let i = 1; i <= N; i++) {
        const atom = getAtom(i);
        atomInfo.push({
            index: i,
            symbol: atom.symbol,
            valence: atom.valence,
            maxBond: atom.maxBond,
            maxH: atom.valence - 1  // Max H per atom (leave at least 1 for bonds)
        });
    }

    const totalH = atoms.totalH;

    // Get all matrix positions to fill (lower triangle)
    const positions = [];
    for (let j = 1; j < N; j++) {
        for (let i = j + 1; i <= N; i++) {
            positions.push({ i, j });
        }
    }

    console.log('Silent generation starting...');
    console.log('N =', N);
    console.log('Total H =', totalH);
    console.log('Positions to fill:', positions.length);

    let hDistributionsChecked = 0;
    let bondConfigurationsChecked = 0;

    /**
     * Generate unique (sorted, non-increasing) H distributions that sum to totalH.
     * Since all atoms of the same type are initially equivalent, we only need
     * to generate canonical distributions where H counts are in non-increasing order.
     * This dramatically reduces duplicates from symmetry.
     */
    function* generateUniqueHDistributions(atomIndex, remainingH, currentDistribution, maxAllowed) {
        if (atomIndex >= N) {
            if (remainingH === 0) {
                yield [...currentDistribution];
            }
            return;
        }

        const maxH = atomInfo[atomIndex].maxH;
        // maxAllowed ensures non-increasing order (each H count <= previous)
        const effectiveMax = Math.min(maxH, remainingH, maxAllowed);
        const remainingAtoms = N - atomIndex - 1;

        // Try from high to low to get non-increasing order
        for (let h = effectiveMax; h >= 0; h--) {
            // Check if remaining atoms can absorb remaining hydrogens
            // (each can take at most min(h, maxH) to maintain non-increasing order)
            const maxRemainingH = remainingAtoms * Math.min(h, maxH);
            if (maxRemainingH >= remainingH - h) {
                currentDistribution[atomIndex] = h;
                yield* generateUniqueHDistributions(atomIndex + 1, remainingH - h, currentDistribution, h);
            }
        }
    }

    /**
     * Check if current bond configuration is valid and connected
     */
    function isValidComplete(freeValences, workMatrix, workBondSum) {
        // Check all valences are satisfied
        for (let k = 1; k <= N; k++) {
            if (workBondSum[k] !== freeValences[k]) {
                return false;
            }
        }

        // Check connectivity using DFS
        const visited = new Array(N + 1).fill(false);
        let visitedCount = 0;

        function dfs(atomIndex) {
            if (visited[atomIndex]) return;
            visited[atomIndex] = true;
            visitedCount++;
            for (let j = 1; j <= N; j++) {
                if ((workMatrix[atomIndex][j] > 0 || workMatrix[j][atomIndex] > 0) && !visited[j]) {
                    dfs(j);
                }
            }
        }

        dfs(1);
        return visitedCount === N;
    }

    /**
     * Generate SMILES from working matrix
     */
    function generateSmiles(workMatrix) {
        try {
            const mol = new Kekule.Molecule();
            const kekuleAtoms = [];

            for (let i = 1; i <= N; i++) {
                const kAtom = new Kekule.Atom();
                kAtom.setSymbol(atomInfo[i - 1].symbol);
                mol.appendNode(kAtom);
                kekuleAtoms.push(kAtom);
            }

            for (let i = 1; i <= N; i++) {
                for (let j = i + 1; j <= N; j++) {
                    const bondOrder = workMatrix[i][j];
                    if (bondOrder > 0) {
                        const bond = new Kekule.Bond();
                        bond.setBondOrder(bondOrder);
                        bond.setConnectedObjs([kekuleAtoms[i - 1], kekuleAtoms[j - 1]]);
                        mol.appendConnector(bond);
                    }
                }
            }

            return Kekule.IO.saveFormatData(mol, 'smi');
        } catch (e) {
            return null;
        }
    }

    /**
     * Enumerate all bond configurations for a given H distribution
     */
    function enumerateBonds(freeValences) {
        // Create working matrix and bondSum arrays
        const workMatrix = [];
        const workBondSum = [];
        for (let i = 0; i <= N; i++) {
            workMatrix[i] = [];
            workBondSum[i] = 0;
            for (let j = 0; j <= N; j++) {
                workMatrix[i][j] = 0;
            }
        }

        function enumerate(posIndex) {
            if (posIndex >= positions.length) {
                // Reached end - check if valid
                bondConfigurationsChecked++;
                if (isValidComplete(freeValences, workMatrix, workBondSum)) {
                    const smiles = generateSmiles(workMatrix);
                    if (smiles) {
                        validSmiles.push(smiles);
                    }
                }
                return;
            }

            const { i, j } = positions[posIndex];
            const maxBondI = atomInfo[i - 1].maxBond;
            const maxBondJ = atomInfo[j - 1].maxBond;
            const maxBond = Math.min(maxBondI, maxBondJ);

            // Try bond orders: 0, 1, 2, 3
            for (let bondOrder = 0; bondOrder <= maxBond; bondOrder++) {
                // Check if this bond order is feasible
                const newBondSumI = workBondSum[i] + bondOrder;
                const newBondSumJ = workBondSum[j] + bondOrder;

                if (newBondSumI <= freeValences[i] && newBondSumJ <= freeValences[j]) {
                    // Apply bond
                    workMatrix[i][j] = bondOrder;
                    workMatrix[j][i] = bondOrder;
                    workBondSum[i] = newBondSumI;
                    workBondSum[j] = newBondSumJ;

                    // Recurse
                    enumerate(posIndex + 1);

                    // Undo bond
                    workMatrix[i][j] = 0;
                    workMatrix[j][i] = 0;
                    workBondSum[i] -= bondOrder;
                    workBondSum[j] -= bondOrder;
                }
            }
        }

        enumerate(0);
    }

    // Enumerate all unique H distributions and for each, enumerate bonds
    // Start with maxAllowed = 3 (the maximum H any carbon can have)
    const initialMaxH = Math.max(...atomInfo.map(a => a.maxH));
    try {
        for (const hDistribution of generateUniqueHDistributions(0, totalH, new Array(N), initialMaxH)) {
            hDistributionsChecked++;

            // Calculate freeValences (1-indexed)
            const freeValences = [0]; // Index 0 unused
            for (let i = 0; i < N; i++) {
                freeValences.push(atomInfo[i].valence - hDistribution[i]);
            }

            // Enumerate all bond configurations for this H distribution
            enumerateBonds(freeValences);
        }
    } catch (e) {
        console.error('Error during enumeration:', e);
    }

    console.log('H distributions checked:', hDistributionsChecked);
    console.log('Bond configurations checked:', bondConfigurationsChecked);
    console.log('Valid structures found:', validSmiles.length);
    console.log('SMILES:', validSmiles);

    // Create and download zip file
    if (validSmiles.length > 0) {
        const zip = new JSZip();
        const smilesContent = validSmiles.join('\n');
        zip.file('structures.smi', smilesContent);

        const blob = await zip.generateAsync({ type: 'blob' });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = `weblucy_structures_${validSmiles.length}.zip`;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);

        showStatus(`Generated ${validSmiles.length} structures (including duplicates) from ${hDistributionsChecked} unique H-distributions. Download started.`, 'success');
    } else {
        showStatus('No valid structures found.', 'warning');
    }

    btn.disabled = false;
    btn.textContent = 'Silent Generation';
}

// Initialize on page load
document.addEventListener('DOMContentLoaded', () => {
    updateControls();
    console.log('WebLucy initialized. Using Natural Products API for 2D depiction.');
});
