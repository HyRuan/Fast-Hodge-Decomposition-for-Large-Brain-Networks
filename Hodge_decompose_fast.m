function [Yg, Yc, Yh] = Hodge_decompose_fast(Yvec, B)
% Decompose a 1-cochain (edge flow) Y into orthogonal components:
% Y = Yg + Yc + Yh
% 
% Gradient: solve L0 * φ = B1 * Y using the precomputed Cholesky
%      decomposition (Requires MATLAB 2017b or newer)
% Curl:     project the remainder Y⊥ = Y − Yg onto Im(B2) by solving
%      (B2*B2') * X = B2*B2' * Y⊥   with PCG,
%      using only operator products Afun(x) = B2 * (B2' * x), i.e., no
%      explicit E×E or T×T matrices are formed. A simple Jacobi preconditioner
%      M^{-1} ≈ diag(B2*B2')^{-1} (edge participation count) is applied.
% 
% INPUT
% Yvec : [E x 1] double. Edge flow stacked in the SAME order as pre.ei/pre.ej
%       (i.e., the order of Skel.EdgeList with i<j).
% B : struct. Incidence matrices containing B1 and B2
% 
% (C) 2025 Oct Hanyang Ruan

    B1 = B.B1;      % [P x E] sparse
    B2 = B.B2;      % [E x T] sparse
    P  = size(B1,1);
    L0 = B1*B1.' + 1e-10*speye(P);
    dec0 = decomposition(L0,'chol'); 

    % Gradient flow Yg
    rhs0 = B1 * Yvec;             % P×1
    phi = dec0 \ rhs0;
    Yg = B1.' * phi;              % E×1

    % Curl flow Yc: edge-space projection (avoid extremely large
    % triangle-space)
    yperp = Yvec - Yg;
    rhs   = B2 * (B2.' * yperp);       % E×1

    % Use only operator products A(x) = B2*(B2'*x) 
    Afun = @(x, varargin) B2 * (B2.' * x);

    % Jacobi precondition：M^{-1} x = x ./ deg_edge
    deg_edge = full(sum(spones(B2), 2)); 
    deg_edge(deg_edge==0) = 1;
    Minv = @(x, varargin) x ./ deg_edge;

    % PCG 
    x0 = zeros(size(B1, 2),1);
    [Yc, ~, ~, ~] = pcg(Afun, rhs, 1e-8, 500, Minv, [], x0);

    % Harmony flow Yh
    Yh = Yvec - Yg - Yc;
end