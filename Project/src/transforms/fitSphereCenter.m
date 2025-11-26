function center = fitSphereCenter(V)
% fitSphereCenter
% ----------------
% Computes the least-squares best-fit sphere centre for a point cloud V (Nx3).
%
% Solves for c = [cx, cy, cz] that minimises the error in:
%       |x - c|^2 = R^2
% in a linearised least-squares sense.
%
% INPUT:
%   V : Nx3 matrix of vertices
%
% OUTPUT:
%   center : 1x3 vector [cx cy cz]

x = V(:,1);
y = V(:,2);
z = V(:,3);

% Build linear LS system
A = [2*x, 2*y, 2*z, ones(size(x))];
b = x.^2 + y.^2 + z.^2;

params = A \ b;
center = params(1:3)';   % only return [cx cy cz]

end
