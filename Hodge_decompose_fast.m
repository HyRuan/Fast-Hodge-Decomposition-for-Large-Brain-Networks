function [Yg, Yc, Yh] = Hodge_decompose_fast(Yvec, B)
% Decompose a 1-cochain (edge flow) Y into orthogonal components:
% Y = Yg + Yc + Yh
% 
% INPUT
% Yvec : [E x 1] double of vectorized connectivity matrix
% B : struct. Incidence matrices containing B1 and B2
% 
% OUTPUT
% Yg [E × 1] double of Gradient flow
% Yc [E × 1] double of Curl flow
% Yh [E × 1] double of Harmonic flow
%
% (C) 2025 Oct Hanyang Ruan

    B1 = B.B1;      % [P x E] sparse
    B2 = B.B2;      % [E x T] sparse
    
    % --- Gradient flow Yg calculation ---
    % Using QR decomposition to achieve faster calculation. 
    % Numerically identical to the lsqr, and has ~1e-15 precision in complete
    % graph while ~1e-6 precision in sparse graph.
    
    A = full(B1.');
    [Q, R, ~] = qr(A, 0);
    tol_qr = max(size(A)) * eps(norm(R, 'inf')); 
    r = sum(abs(diag(R)) > tol_qr);
    Qr = Q(:, 1:r);
    Yg = Qr * (Qr' * Yvec);
    
    % --- Curl flow Yc calculation ---
    % Define operator instead of construct a full L1, which can be extremely
    % large for big FC matrices. This has ~1e-16 precision in complete graph
    % while ~1e-12 precision in sparse graph.
    Afun_T = @(x, varargin) B2.' * (B2 * x);
    rhs_T = B2.' * Yvec;
    [z, ~, ~, ~] = lsqr(Afun_T, rhs_T);
    Yc = B2 * z;

    % Harmony flow Yh
    Yh = Yvec - Yg - Yc;
end
