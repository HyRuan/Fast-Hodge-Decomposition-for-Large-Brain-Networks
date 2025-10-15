% Create a toy 150x150 FC matrix
% The test was performed using an AMD Ryzen 9 16-core CPU and 128 GB of RAM
n = 150;
rm = 2 * rand(n) - 1;
rm_fc = triu(rm) + triu(rm, 1)';
rm_fc(1:n+1:end) = 0;

% ------- Original version -------
tS0 = tic;
Skel0 = Hodge_2Skeleton(rm_fc);
Bmat0 = Hodge_incidence(Skel0);
% A matrix with larger size will cause MATLAB unresponsive problem by 
% creating an array that exceeds maximum array size (128 GB)
Yvec0 = Hodge_vec(rm_fc);
[Yg0, Yc0, Yh0] = Hodge_decompose(Yvec0, Bmat0);
Gmat0 = Hodge_project(Yg0, Skel0);
Cmat0 = Hodge_project(Yc0, Skel0);
toc(tS0);
% Run time: 2925.9982s

% ------- Fast version -------
% Note: the outputs are all arranged as structs for clearer referencing
tS = tic;
Skel  = Hodge_2Skeleton_fast(rm_fc);
Bmat  = Hodge_incidence_sparse(Skel);
Yvec = Hodge_vec(rm_fc); 
[Yg,Yc,Yh] = Hodge_decompose_fast(Yvec, Bmat);
Gmat = Hodge_project_fast(Yg, Skel);
Cmat = Hodge_project_fast(Yc, Skel);
toc(tS);
% Run time: 0.0225s

% ------- Validation -------
% Compare results from original and from fast versions
tolerance = 1e-12;
disp(isequal(size(Gmat), size(Gmat0)) && max(abs(Gmat(:) - Gmat0(:))) <= tolerance); % 1
disp(isequal(size(Cmat), size(Cmat0)) && max(abs(Cmat(:) - Cmat0(:))) <= tolerance); % 1
% Expect Yh to be ~1e-12…1e-13 (not machine-zero) due to iterative tolerance
% and the small λ0 added in L0 for stability.
disp(isequal(size(Yh), size(Yh0)) && max(abs(Yh(:) - Yh0(:))) <= tolerance); % 1