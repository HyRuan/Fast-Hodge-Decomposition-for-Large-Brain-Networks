# Fast Hodge Decomposition for Large Brain Networks
Based on Hodge decomposition method (https://github.com/laplcebeltrami/hodge) from **Vijay D. Anand** and **M.K. Chung**. This modified implementation of Hodge decomposition is designed to efficiently handle large, fully connected functional brain networks (e.g., with 300+ nodes) where conventional formulations become computationally prohibitive. The goal is to decompose an edge flow on a graph into its gradient, curl, and harmonic components, while keeping both memory usage and runtime within practical limits on standard hardware (e.g., 16 GB RAM).

## Practical Notes
- Use this version when processing dense functional connectivity matrices or any complete graph (e.g., fMRI correlation networks). You can simply replace the previous functions with same names.
- The harmonic component may show values on the order of 1e-12 due to iterative tolerance—these can safely be thresholded to zero if needed.
- PCG settings are set default as 1e-8 relative tolerance and 500 iterations. Modify if needed.
- Requires MATLAB R2017b or later for efficient decomposition() and pcg() functions.

(C) 2025 Oct. Hanyang Ruan. Technical University of Munich
