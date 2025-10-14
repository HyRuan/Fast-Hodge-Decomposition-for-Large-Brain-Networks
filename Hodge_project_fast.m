function M = Hodge_project_fast(Yvec, skel)
% Project [E x 1] gradient or curl vector back to [P x P] matrix
% INPUT
% Yvec: [E x 1] double. Gradient or curl flow vector
% skel: Struct. Node, edge, triangle lists.
% 
% OUTPUT
% M [P x P] full matrix of gradient or curl flow
% 
% (C) 2025 Oct Hanyang Ruan
    ei = skel.EdgeList(:,1);  ej = skel.EdgeList(:,2);  
    P = size(skel.NodeList, 1);
    M_upper = sparse(double(ei), double(ej), Yvec, P, P);
    M = M_upper + M_upper.';
    M(1:P+1:end) = 0;
    M = full(M);

end
