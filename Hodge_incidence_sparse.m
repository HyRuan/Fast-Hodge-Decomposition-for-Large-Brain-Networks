function B = Hodge_incidence_sparse(Skel)
% Takes a pre-computed 2-skeleton (struct) and creates sparse Hodge
% incidence matrices. Create sparse B2 matrix to overcome potential out-of-
% memory or "requested array exceeds maximum array size preference" error 
% when the connectivity matrix is very large. Also runs faster.
% 
% INPUT
% Skel (Struct): Including node, edge, triangle lists.
%
% OUTPUT
% B(Struct):
%   B0: empty list (theoretically doesn't exist)
%   B1: [P x E] sparse node-edge incidence
%   B2: [E x T] sparse edge-triangle incidence
% 
% (C) 2025 Hanyang Ruan
%          Technical University of Munich

P  = size(Skel.NodeList, 1);  % Number of nodes
Elist = Skel.EdgeList;  % [E x 2], i<j
Tlist = Skel.TriList;   % [T x 3], i<j<k

E = size(Elist,1);
T = size(Tlist,1);

% B1: Node - Edge. Direction: i->j
i_ind = double(Elist(:,1));
j_ind = double(Elist(:,2));
col = (1:E).';
B1 = sparse([i_ind; j_ind], [col; col], [-ones(E,1); +ones(E,1)], P, E);

% A look up table EID for edge signs 
% EID(a,b) = +e  if the e-th edge has a<b
%           = -e  if the e-th edge has b<a
EID = sparse(P,P);
eid = (1:E)';
EID = EID + sparse(Elist(:,1), Elist(:,2),  eid, P, P);
EID = EID + sparse(Elist(:,2), Elist(:,1), -eid, P, P);

% B2: Edge - Triangle
i = Tlist(:,1); j = Tlist(:,2); k = Tlist(:,3);

e_ij = full(EID(sub2ind([P,P], i, j)));
e_jk = full(EID(sub2ind([P,P], j, k)));
e_ik = full(EID(sub2ind([P,P], i, k)));

row_ij = abs(e_ij);  s_ij = sign(e_ij);
row_jk = abs(e_jk);  s_jk = sign(e_jk);
row_ik = abs(e_ik);  s_ik = sign(e_ik);

rows = [row_ij; row_jk; row_ik];
cols = [(1:T)'; (1:T)'; (1:T)'];
vals = [+s_ij; +s_jk; -s_ik];

B2 = sparse(rows, cols, vals, E, T);

B = struct('B0', [], 'B1', B1, 'B2', B2);

end
