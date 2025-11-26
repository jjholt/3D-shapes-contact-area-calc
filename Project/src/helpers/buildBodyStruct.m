function body = buildBodyStruct(F, V)
% buildBodyStruct
% ---------------
% Build a struct for a rigid body mesh with:
%   F          : faces (Nx3)
%   V          : vertices (Mx3)
%   triCentroid: per-face centroid
%   triNormal  : per-face unit normal
%   triArea    : per-face area
%   bbox       : 3x2 [xmin xmax; ymin ymax; zmin zmax]
%
% Also creates alias fields:
%   centroids, normals, areas
% so other functions can use more generic names.

body.F = F;
body.V = V;

% Triangle vertices
v1 = V(F(:,1),:);
v2 = V(F(:,2),:);
v3 = V(F(:,3),:);

% Centroids
triCentroid = (v1 + v2 + v3) / 3;
body.triCentroid = triCentroid;
body.centroids   = triCentroid;   % alias

% Normals and areas
e1 = v2 - v1;
e2 = v3 - v1;
n  = cross(e1, e2, 2);       % non-normalised normals
area = 0.5 * sqrt(sum(n.^2,2));

% Avoid division by zero
nz = vecnorm(n,2,2);
nz(nz==0) = 1;
triNormal = n ./ nz;

body.triNormal = triNormal;
body.normals   = triNormal;       % alias

body.triArea = area;
body.areas   = area;              % alias

% Bounding box
xmin = min(V(:,1)); xmax = max(V(:,1));
ymin = min(V(:,2)); ymax = max(V(:,2));
zmin = min(V(:,3)); zmax = max(V(:,3));

body.bbox = [xmin xmax; ymin ymax; zmin zmax];
end
