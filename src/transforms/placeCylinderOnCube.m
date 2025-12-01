function [V_cyl_aligned, shift] = placeCylinderOnCube(V_cyl, V_cube)
% placeCylinderOnCube
% -------------------
% Translate a cylinder mesh so it sits on top of a cube mesh.
%
% ASSUMPTIONS:
%   - Cylinder axis is approximately the global Z-axis.
%   - Cylinder is "vertical" in the STL (no tilting).
%   - Cube top face is the face with maximum Z coordinate.
%
% INPUT:
%   V_cyl  : Mx3 vertices of cylinder
%   V_cube : Nx3 vertices of cube
%
% OUTPUT:
%   V_cyl_aligned : translated cylinder vertices
%   shift         : 1x3 translation vector applied to V_cyl

    % --- 1) Find cube centre (overall centroid) ---
    centerCube = mean(V_cube, 1);  % [cx, cy, cz]

    % --- 2) Find cylinder centroid (overall) ---
    centerCyl  = mean(V_cyl, 1);   % [cx, cy, cz]

    % We will align XY using centroids (good for symmetric shapes).
    shiftXY = centerCube(1:2) - centerCyl(1:2);

    % --- 3) Find cube top (max Z) and cylinder bottom (min Z) ---
    zTopCube    = max(V_cube(:,3));   % top surface of cube
    zBottomCyl  = min(V_cyl(:,3));    % bottom of cylinder

    % We want: (zBottomCyl + shiftZ) = zTopCube  ->  shiftZ = zTopCube - zBottomCyl
    shiftZ = zTopCube - zBottomCyl;

    % --- 4) Combine full translation ---
    shift = [shiftXY, shiftZ];

    % --- 5) Apply translation to cylinder vertices ---
    V_cyl_aligned = V_cyl + shift;
end
