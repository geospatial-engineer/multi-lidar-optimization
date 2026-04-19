# Symmetric Multi-LiDAR Mount Configuration Optimization

This repository contains MATLAB source code for mixed discrete–continuous optimization of symmetric multi-LiDAR mount configurations for mobile mapping systems. The current release is code-only and focuses on the OS1-64 optimization workflows implemented with Bayesian Optimization, Genetic Algorithm, and Particle Swarm Optimization.

## Repository contents

The main MATLAB entry points are stored in the `code/` directory.

`optimize_OS1_BO.m` implements the Bayesian Optimization workflow.

`optimize_OS1_GA.m` implements the Genetic Algorithm workflow.

`optimize_OS1_PSO.m` implements the Particle Swarm Optimization workflow.

These scripts include trajectory-based voxel simulation, overlap regularization, angular-diversity scoring, symmetry expansion across the YZ plane, and export of optimized configurations.

## Scope of this release

This repository is intended as a lightweight code release that can accompany a manuscript submission or accepted article. At the moment, it does not bundle datasets, figures, or large result files. If desired, those can be added later in a structured release with versioned archives.

## Recommended repository structure

The current structure is ready for GitHub Pages publishing from the repository root.

The landing page is `index.html`.

The stylesheet is `assets/styles.css`.

The MATLAB code is stored under `code/`.

## How to publish with GitHub Pages

Create a new GitHub repository and upload the full contents of this folder.

In the repository settings, open the Pages section and set the source to deploy from the main branch root.

Once Pages is enabled, GitHub will publish `index.html` automatically.

## Citation

If you use this repository in academic work, please cite the associated manuscript and, if available, the repository archive DOI.

Suggested software citation:

Author(s). Symmetric Multi-LiDAR Mount Configuration Optimization: MATLAB Code Release. GitHub repository.

## License

No license file is included yet. Add the license you prefer before making the repository public.
