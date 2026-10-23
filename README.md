# Topology-Informed Indices for Neurovascular Data

This repository provides a Python implementation of four topology-informed indices for the quantitative analysis of 3D neurovascular structures:

- **Bifurcation loop index (β):** Characterizes loop structures.
- **Bifurcation angle index (α):** Quantifies potential angle-forming locations at branching junctions.
- **Curvature index (c):** Characterizes multiscale geometric curvature.
- **Entanglement index (e):** Quantifies spatial proximity between vascular and neural structures.

## Requirements

- Python 3
- NumPy
- GUDHI

Install the required packages using:

```bash
pip install numpy gudhi
```

## Usage

Run the following command to compute the indices for the example HOC file:

```bash
python ./code/test.py ex.hoc
```

To analyze your own data, replace `ex.hoc` with the path to your HOC file:

```bash
python ./code/test.py path/to/your_data.hoc
```

The script computes and prints the bifurcation loop, bifurcation angle, curvature, and entanglement indices.

## Input Data

The input should be a `.hoc` file containing 3D filament data, with each filament represented as an ordered sequence of spatial coordinates.
