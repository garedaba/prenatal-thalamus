function G = GfromSurf(verts,faces,parc)

edges = [
    faces(:, [1 2]);
    faces(:, [2 3]);
    faces(:, [3 1])
];

% Sort vertex indices within each edge (ensures unique ordering)
edges = sort(edges, 2);

% Remove duplicate edges
edges_unique = unique(edges, 'rows');

% I think it still keeps edges between a vertex in a region of a
% parcellation and an unassigned region because in same cases it resulted
% in disconnected vertices which messed up things later on. This was the
% easiest fix
if nargin > 2
    edge2keep = sum(parc(edges_unique),2)>0;
    edges_unique(~edge2keep,:) = [];
end

% --- Optional: compute edge lengths ---
edge_lengths = sqrt(sum((verts(edges_unique(:,1),:) - verts(edges_unique(:,2),:)).^2, 2));

% --- Combine into weighted edge list ---
edge_list = [edges_unique, edge_lengths];

% Example: create MATLAB graph object
G = graph(edge_list(:,1), edge_list(:,2), edge_list(:,3));