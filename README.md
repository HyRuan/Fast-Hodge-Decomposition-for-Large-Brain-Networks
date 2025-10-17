# Fast Hodge Decomposition for Large Brain Networks
Based on Hodge decomposition method (https://github.com/laplcebeltrami/hodge) from **Vijay D. Anand** and **M.K. Chung**. This modified implementation of Hodge decomposition is designed to efficiently handle large, fully connected functional brain networks (e.g., with 300+ nodes) where conventional formulations become computationally prohibitive. It maintains mathematical equivalence with the original algorithm described in _Hodge Decomposition of Functional Human Brain Networks_ while optimizing for scalability and numerical stability. The goal is to decompose an edge flow on a graph into its gradient, curl, and harmonic components, while keeping both memory usage and runtime within practical limits on standard hardware (e.g., 16 GB RAM).

## Practical Notes
- Use this version when processing dense functional connectivity matrices or any complete graph (e.g., fMRI correlation networks). You can simply replace the previous functions with same names.
- Oct 17 update: now the codes can also properly handle sparse graph, with tolerable error (max 1e-06).

(C) 2025 Oct. Hanyang Ruan. Technical University of Munich
