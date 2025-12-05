function [F, V] = loadStlMesh(filename, tag)
% loadStlMesh  Wrapper around MATLAB's stlread.
%
%   [F, V] = loadStlMesh(filename, tag)
%
% Returns:
%   F : (#faces x 3) face connectivity
%   V : (#verts x 3) vertices

    if nargin < 2
        tag = 'mesh';
    end

    fprintf('Loading "%s" as %s...\n', filename, tag);

    if exist(filename, 'file') ~= 2
        error('Cannot find STL file: %s', filename);
    end

    % Modern stlread returns a triangulation
    TR = stlread(filename);

    if isa(TR, 'triangulation') || isa(TR, 'delaunayTriangulation')
        F = double(TR.ConnectivityList);
        V = double(TR.Points);
    else
        % Fallback for older signature [F,V] = stlread(file)
        [F, V] = stlread(filename);
        F = double(F);
        V = double(V);
    end
end
