function Skel = Hodge_2Skeleton_fast(A_bin)
% Create 2-skeleton from adjacency matrix. Include a fast track for
% complete graph.
%
% INPUT:
% A_bin: Binary adjacency matrix (accept sparse matrix) for finding edges
% and triangles.
%
% OUTPUT
% Skel (Struct): Node, edge, triangle lists.
%
% (C) 2025 Oct Hanyang Ruan

    if ~issparse(A_bin), A_bin = sparse(A_bin); end
    A_bin = spones(A_bin);
    P = size(A_bin,1);
    U = triu(A_bin,1);
    E = nnz(U);
    is_complete = (E == P*(P-1)/2);

    if is_complete
        % ------- Fast track for a complete graph -------
        % Edge: All pairwise combinations with i<j
        EdgeList = double(nchoosek(double(1:P), 2));  % [E x 2]
        % Triangle: All triplets with i<j<k
        TriList  = double(nchoosek(double(1:P), 3));  % [T x 3]
    else
        % ------- Sparse graph track -------
        % Edge list 
        [i_ind, j_ind] = find(U);                         % i<j
        EdgeList = double([i_ind j_ind]);
        EdgeList = sparse(sortrows(EdgeList, [1 2]));
        % Triangle list 
        TriList = zeros(0,3,'double');
        chunk = 5e5; buf = zeros(chunk,3,'double'); ptr = 0;
        N = cell(P,1);
        for u = 1:P
            N{u} = double(find(U(u,u+1:end)~=0) + u);
        end

        for i = 1:P-2
            Ni = N{i};
            if numel(Ni) < 2, continue; end
            for a = 1:numel(Ni)
                v = Ni(a);
                Nv = N{v};
                Ni_gtv = Ni(Ni>v);
                if isempty(Ni_gtv) || isempty(Nv), continue; end
                w = intersect(Ni_gtv, Nv); 
                k = numel(w);
                if k==0, continue; end
                if ptr + k > size(buf,1)
                    buf = [buf; zeros(chunk,3,'double')];
                end
                buf(ptr+1:ptr+k, :) = [double(i)*ones(k,1,'double'), double(v)*ones(k,1,'double'), w(:)];
                ptr = ptr + k;
            end
        end
        TriList = buf(1:ptr,:);
    end
    Skel = struct('NodeList', double(1:P)', ...
                  'EdgeList', EdgeList, ...
                  'TriList',  TriList);
end

