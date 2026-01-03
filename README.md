# WebLucy

A web-based cheminformatics teaching tool for exploring molecular connection matrices and understanding structure elucidation algorithms.

## Overview

WebLucy is an interactive webapp designed to teach students how molecular structures can be systematically generated from a molecular formula. It visualizes the process of bond insertion in a connection matrix, demonstrating the depth-first search algorithm used in structure elucidation.

The algorithm is based on the original LUCY (semiautomatische Strukturaufklärung) software by Christoph Steinbeck (1994-1995). The bond generation algorithm is extremely basic and not suited for the generation of larger molecules due to the lack of canonicalization. In the context of NMR-generated contrains, however, it works reasonably well. 

## Features

- **Molecular Formula Input**: Enter any molecular formula (e.g., C5H12, C2H6O, C6H6)
- **Connection Matrix Visualization**: See the NxN matrix of heavy atoms with color-coded bond orders
- **Editable Hydrogen Distribution**: Adjust how hydrogens are distributed among heavy atoms
- **Step-by-Step Algorithm**: Walk through the bond insertion process with Forward/Back buttons
- **Auto-Step Mode**: Automatically step through configurations
- **Connectivity Validation**: Only connected molecules are accepted as valid structures
- **2D Structure Depiction**: Valid molecules are displayed using the Natural Products API

## How It Works

1. **Enter a molecular formula** - The app parses the formula and creates a connection matrix for heavy atoms (non-hydrogen atoms)

2. **Adjust hydrogen distribution** - Edit how many hydrogens are attached to each atom. The sum must match the formula.

3. **Explore bond configurations** - Use Forward/Back buttons to step through the depth-first search:
   - The algorithm traverses the lower triangle of the matrix column by column
   - At each position, it tries bond orders: 1, 2, 3, then 0
   - Valence constraints are checked at each step
   - When the matrix is complete, connectivity is verified

4. **View valid structures** - When a valid, connected structure is found, it's displayed as a 2D diagram

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
- [Natural Products API](https://api.naturalproducts.net/latest/docs) - For 2D structure depiction

## Algorithm

The bond insertion algorithm performs a depth-first traversal of possible bond configurations:

1. Start at matrix position (2,1)
2. Try bond order 1, then proceed to next position
3. Continue until the matrix is complete
4. Check if all valences are satisfied and the molecule is connected
5. If valid, display the structure
6. Backtrack and try higher bond orders (2, 3) or no bond (0)
7. Repeat until all configurations are explored

## License

MIT License

## Author

Christoph Steinbeck

## Acknowledgments

- Original LUCY algorithm (1994-1995)
- [Kekule.js](https://partridgejiang.github.io/Kekule.js/) by Partridge Jiang
- [Cheminformatics Microservices](https://api.naturalproducts.net) by the Steinbeck Group
