# WebLucy

A web-based cheminformatics teaching tool for exploring molecular connection matrices and understanding structure elucidation algorithms.

**[Try WebLucy Online](https://steinbeck.github.io/weblucy/)**

## Overview

WebLucy is an interactive webapp designed to teach students how molecular structures can be systematically generated from a molecular formula. It visualizes the process of bond insertion in a connection matrix, demonstrating the depth-first search algorithm used in structure elucidation.

The algorithm is based on the original LUCY (semiautomatische Strukturaufklärung) software by Christoph Steinbeck (1994-1995). The bond generation algorithm is extremely basic and not suited for the generation of larger molecules due to the lack of canonicalization. In the context of NMR-generated constraints, however, it works reasonably well.

For production-level structure generation from molecular formulas, consider using [SURGE](https://github.com/StructureGenerator/surge) (Structure Generator), a highly efficient and comprehensive tool for exhaustive structure enumeration with advanced canonicalization and symmetry handling.

## Features

### Core Functionality
- **Molecular Formula Input**: Enter any molecular formula (e.g., C5H12, C2H6O, C6H6)
- **Connection Matrix Visualization**: See the NxN matrix of heavy atoms with color-coded bond orders
- **Editable Hydrogen Distribution**: Freely adjust how hydrogens are distributed among heavy atoms
- **Step-by-Step Algorithm**: Walk through the bond insertion process with Forward/Back buttons
- **Connectivity Validation**: Only connected molecules are accepted as valid structures
- **2D Structure Depiction**: Valid molecules are displayed using the Natural Products API
- **Structure Counter**: Track how many valid structures have been found (shown in panel heading)

### Auto-Step Mode
- **Continuous Exploration**: Auto-step continues through all bond configurations until exhausted
- **Adjustable Speed**: Choose from Slow (500ms), Medium (250ms), Fast (100ms), or Very Fast (50ms)
- **Stop Anytime**: Click the Stop button to pause auto-stepping at any point

### Bond Order Options
- **Ascending Order (1→2→3)**: Try single bonds first, then double, then triple (default)
- **Descending Order (3→2→1)**: Try triple bonds first, then double, then single
- Useful for demonstrating how bond order affects structure discovery

### Silent Generation
- **Batch Generation**: Generate all valid structures for the current formula
- **Unique H-Distributions**: Automatically enumerates all canonical hydrogen distributions
- **SMILES Export**: Downloads all structures as a zipped SMILES file
- **Includes Duplicates**: Intentionally preserves duplicates to demonstrate the need for canonicalization

## How It Works

1. **Enter a molecular formula** - The app parses the formula and creates a connection matrix for heavy atoms (non-hydrogen atoms)

2. **Adjust hydrogen distribution** - Edit how many hydrogens are attached to each atom. The sum must match the formula (validated when you click Forward or Auto Step).

3. **Select options** - Choose auto-step speed and bond order sequence (ascending or descending)

4. **Explore bond configurations** - Use Forward/Back buttons or Auto Step to traverse the depth-first search:
   - The algorithm traverses the lower triangle of the matrix column by column
   - At each position, it tries bond orders in the selected sequence
   - Valence constraints are checked at each step
   - When the matrix is complete, connectivity is verified

5. **View valid structures** - When a valid, connected structure is found, it's displayed as a 2D diagram with a structure count

6. **Silent Generation** - Or click "Silent Generation" to enumerate all valid structures and download them as SMILES

## Running Locally

Simply open `index.html` in a web browser:

```bash
open index.html
```

Or serve via a local HTTP server:

```bash
python3 -m http.server 8000
```

Then open http://localhost:8000

## Dependencies

- [Kekule.js](https://partridgejiang.github.io/Kekule.js/) - For molecule representation and SMILES generation
- [JSZip](https://stuk.github.io/jszip/) - For creating downloadable zip files
- [Natural Products API](https://api.naturalproducts.net/latest/docs) - For 2D structure depiction

## Algorithm

The bond insertion algorithm performs a depth-first traversal of possible bond configurations:

1. Start at matrix position (2,1)
2. Try bond orders in the selected sequence (1→2→3→0 or 3→2→1→0)
3. If a bond value exceeds available valence, try the next value in sequence
4. Continue until the matrix is complete
5. Check if all valences are satisfied and the molecule is connected
6. If valid, display the structure and increment the counter
7. Backtrack and try the next bond order in sequence
8. Repeat until all configurations are explored

## License

MIT License

## Author

Christoph Steinbeck

## Acknowledgments

- Original LUCY algorithm (1994-1995)
- [Kekule.js](https://partridgejiang.github.io/Kekule.js/) by Partridge Jiang
- [JSZip](https://stuk.github.io/jszip/) by Stuart Knightley
- [Cheminformatics Microservices](https://api.naturalproducts.net) by the Steinbeck Group
