function [F, V] = loadStlMesh(filename, label)
% loadStlMesh
% -----------
% Loads an STL file (ASCII or binary) and returns:
%   F : faces (Nx3)
%   V : vertices (Mx3)
%
% If the STL has duplicated vertices (very common), they are merged.
%
% INPUT:
%   filename : string, path to the STL
%   label    : optional, printed label for info
%
% OUTPUT:
%   F, V

if nargin < 2
    label = 'mesh';
end

fprintf('Loading "%s" as %s...\n', filename, label);

% Use built-in stlread if available (R2022b+)
useBuiltin = false;
if exist('stlread','file') == 2
    try
        TR = stlread(filename);
        V = TR.Points;
        F = TR.ConnectivityList;
        useBuiltin = true;
    catch
        useBuiltin = false;
    end
end

if ~useBuiltin
    % ---- Manual ASCII/binary STL reader ----
    [F_raw, V_raw] = readSTLmanual(filename);
    
    % Deduplicate vertices
    [V, ~, idx] = uniquetol(V_raw, 1e-12, 'ByRows', true);
    F = idx(F_raw);
end

fprintf('Loaded: %d vertices, %d faces\n', size(V,1), size(F,1));

end


%% === Manual STL reader (ASCII or binary) ============================
function [F, V] = readSTLmanual(filename)
fid = fopen(filename, 'rb');
if fid < 0
    error('Cannot open STL file: %s', filename);
end

% Detect ASCII vs binary
firstLine = fgetl(fid);
frewind(fid);

if contains(lower(firstLine), 'solid')
    % Could still be binary, but assume ASCII and try
    try
        [F, V] = readSTL_ascii(fid);
        fclose(fid);
        return;
    catch
        % fallback to binary
        frewind(fid);
    end
end

% Binary read fallback
[F, V] = readSTL_binary(fid);
fclose(fid);
end


function [F,V] = readSTL_ascii(fid)
% Preallocate with a reasonable starting capacity and grow in chunks
initFaceCap = 10000;          % initial guess for number of faces
initVertCap = initFaceCap*3;  % 3 vertices per face

V = zeros(initVertCap,3);     % preallocate vertices
F = zeros(initFaceCap,3);     % preallocate faces

nV = 0;   % number of vertices actually stored
nF = 0;   % number of faces actually stored

while ~feof(fid)
    line = strtrim(fgetl(fid));
    if ~ischar(line) && ~isstring(line)
        break;
    end

    % Look for 'facet normal ...'
    if startsWith(line, 'facet', 'IgnoreCase', true)
        % We know a facet has exactly 3 vertices, then 'endfacet'
        % Read until we've seen 3 'vertex' lines and then 'endfacet'
        localVertIdx = zeros(3,1);   % indices of this facet's vertices
        vCount = 0;

        while ~feof(fid)
            innerLine = strtrim(fgetl(fid));
            if ~ischar(innerLine) && ~isstring(innerLine)
                break;
            end

            if startsWith(innerLine, 'vertex', 'IgnoreCase', true)
                % Parse vertex
                nums = sscanf(innerLine, 'vertex %f %f %f');
                if numel(nums) ~= 3
                    error('Malformed vertex line in ASCII STL.');
                end

                % Ensure vertex capacity
                if nV + 1 > size(V,1)
                    % grow V by a chunk
                    V = [V; zeros(size(V,1),3)]; %#ok<AGROW>
                end

                nV = nV + 1;
                V(nV,:) = nums.';
                localVertIdx(vCount+1) = nV;
                vCount = vCount + 1;

            elseif startsWith(innerLine, 'endfacet', 'IgnoreCase', true)
                % End of this facet: record the face if we got 3 vertices
                if vCount ~= 3
                    error('Facet without 3 vertices in ASCII STL.');
                end

                if nF + 1 > size(F,1)
                    % grow F by a chunk
                    F = [F; zeros(size(F,1),3)]; %#ok<AGROW>
                end

                nF = nF + 1;
                F(nF,:) = localVertIdx.';
                break;  % go back to outer while to look for next facet
            else
                % other lines (outer loop, endloop, etc.) – ignore
            end
        end
    end
end

% Trim unused preallocated rows
V = V(1:nV,:);
F = F(1:nF,:);
end
