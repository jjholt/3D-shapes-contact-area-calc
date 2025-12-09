function [femur, tibia, patella, R_global, t_global] = setupKneeFromStl(modelDir, thetaFlexDeg)
% setupKneeFromStl
% Load Femur.STL, Tibia.STL, Patella.STL (from separate files),
% put them into a common anatomical-ish frame, and pose the knee
% at a chosen flexion angle.
%
% INPUTS
%   modelDir      : folder containing Femur.STL, Tibia.STL, Patella.STL
%   thetaFlexDeg  : knee flexion angle in degrees (femur relative to tibia).
%                   Positive value -> femur flexes posteriorly.
%                   Default = 30.
%
% OUTPUTS
%   femur, tibia, patella: structs with fields .V (Nx3) and .F (Mx3)
%                          in the final posed global frame.
%   R_global, t_global   : rotation and translation used to align with tibia

    if nargin < 1 || isempty(modelDir)
        scriptDir   = fileparts(mfilename('fullpath'));
        projectRoot = fileparts(scriptDir); % assume this file is in src/
        modelDir    = fullfile(projectRoot, 'model');
    end
    if nargin < 2 || isempty(thetaFlexDeg)
        thetaFlexDeg = 30;  % nice mid-flexion to show contact
    end

    %% 1. Load STL meshes in their original frames
    fprintf('Loading STL meshes...\n');
    [femur.F,   femur.V]   = loadAnyStl(fullfile(modelDir, 'Femur.STL'));
    [tibia.F,   tibia.V]   = loadAnyStl(fullfile(modelDir, 'Tibia.STL'));
    [patella.F, patella.V] = loadAnyStl(fullfile(modelDir, 'Patella.STL'));

    %% 2. Global centering (just to avoid huge coordinates)
    allVerts = [femur.V; tibia.V; patella.V];
    t_global = mean(allVerts, 1);   % centroid of whole set

    femur.V   = femur.V   - t_global;
    tibia.V   = tibia.V   - t_global;
    patella.V = patella.V - t_global;

    %% 3. Use tibia PCA to define global axes (like before)
    Vt = tibia.V;
    Ct = cov(Vt);
    [E, D] = eig(Ct);
    [~, idx] = sort(diag(D), 'descend');
    E = E(:, idx);              % principal directions

    shaftAxis = E(:,1);         % tibial shaft
    other1    = E(:,2);
    other2    = E(:,3);

    proj1 = Vt * other1;
    proj2 = Vt * other2;
    range1 = max(proj1) - min(proj1);
    range2 = max(proj2) - min(proj2);

    if range1 >= range2
        yAxis = other1;
        xAxis = cross(yAxis, shaftAxis);
    else
        yAxis = other2;
        xAxis = cross(yAxis, shaftAxis);
    end

    xAxis = xAxis / norm(xAxis);
    yAxis = yAxis / norm(yAxis);
    zAxis = shaftAxis / norm(shaftAxis);

    R_global = [xAxis, yAxis, zAxis];
    if det(R_global) < 0
        xAxis   = -xAxis;
        R_global = [xAxis, yAxis, zAxis];
    end

    % Rotate ALL bones into this common frame
    femur.V   = femur.V   * R_global;
    tibia.V   = tibia.V   * R_global;
    patella.V = patella.V * R_global;

    %% 4. Place tibia: plateau at z = 0
    tibZ      = tibia.V(:,3);
    zPlateau  = max(tibZ);                 % "top" of tibia
    tibia.V(:,3) = tibia.V(:,3) - zPlateau;

    %% 5. Place femur just above tibial plateau
    % Small cartilage gap (same units as STL: guess ~2)
    gap = 2.0;

    [fMin, fMax] = bbox3D(femur.V);
    % Decide which end is distal (closer to tibial plateau z=0)
    distToMin = abs(fMin(3) - 0);
    distToMax = abs(fMax(3) - 0);
    if distToMin <= distToMax
        distalZ   = fMin(3);
        shiftZ    = (0 + gap) - distalZ;   % move distal surface to gap above plateau
    else
        distalZ   = fMax(3);
        shiftZ    = (0 + gap) - distalZ;
    end
    femur.V(:,3) = femur.V(:,3) + shiftZ;

    % Align femur centroid in X and Y with tibia centroid
    tibCent = mean(tibia.V,  1);
    femCent = mean(femur.V,  1);
    dXY     = tibCent(1:2) - femCent(1:2);
    femur.V(:,1:2) = femur.V(:,1:2) + dXY;

    %% 6. Place patella anterior to femur, roughly in trochlear region
    [tMin, tMax] = bbox3D(tibia.V);
    [fMin, fMax] = bbox3D(femur.V);
    [pMin, pMax] = bbox3D(patella.V);

    % X: centre with femur
    patCent = mean(patella.V, 1);
    femCent = mean(femur.V,   1);
    dx = femCent(1) - patCent(1);

    % Z: somewhere between tibial plateau and mid-femur (trochlear region)
    femDistalZ = min(fMin(3), fMax(3));   % distal side after placing above tibia
    femProxZ   = max(fMin(3), fMax(3));
    zTrochlea  = femDistalZ + 0.3*(femProxZ - femDistalZ);  % ~30% up from condyles
    dz = zTrochlea - patCent(3);

    % Y: a bit in front of femur
    femFrontY = max(fMax(2));             % anterior-most femur
    femBackY  = min(fMin(2));
    femAPth   = femFrontY - femBackY;
    anteriorOffset = 0.2 * femAPth + gap; % some clearance in front of trochlea
    dy = (femFrontY + anteriorOffset) - patCent(2);

    patella.V = patella.V + [dx dy dz];

    %% 7. Apply knee flexion: rotate femur + patella about ML axis
    % Axis: global X through a pivot near the centre of the plateau
    tibCent = mean(tibia.V, 1);
    pivot   = [tibCent(1), tibCent(2), 0];   % point on the knee joint line

    theta = deg2rad(thetaFlexDeg);
    % Rotation matrix for ROW VECTORS about X axis
    Rx = [1 0          0;
          0 cos(theta)  sin(theta);
          0 -sin(theta) cos(theta)];

    % Rotate femur
    Vrel = femur.V - pivot;
    femur.V = Vrel * Rx + pivot;

    % Rotate patella (it tracks with the femur in this simple model)
    VrelP = patella.V - pivot;
    patella.V = VrelP * Rx + pivot;

    %% 8. Visualise the posed knee
    figure; hold on; axis equal; grid on;
    patch('Faces', tibia.F,   'Vertices', tibia.V,   ...
          'FaceColor',[0.8 0.8 0.9],'EdgeColor','none','FaceAlpha',0.7);
    patch('Faces', femur.F,   'Vertices', femur.V,   ...
          'FaceColor',[0.9 0.7 0.7],'EdgeColor','none','FaceAlpha',0.7);
    patch('Faces', patella.F, 'Vertices', patella.V, ...
          'FaceColor',[0.7 0.9 0.7],'EdgeColor','none','FaceAlpha',0.7);

    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('Knee posed at %.1f° flexion', thetaFlexDeg));
    view(3); camlight; lighting gouraud; rotate3d on;

    fprintf('Done: tibia fixed, femur + patella posed at %.1f° flexion.\n', thetaFlexDeg);
end

%% Helper: axis-aligned bounding box
function [minV, maxV] = bbox3D(V)
    minV = min(V,[],1);
    maxV = max(V,[],1);
end

%% Helper: use your robust STL reader
function [F, V] = loadAnyStl(filename)
    if exist('loadStlMesh','file') == 2
        [F, V] = loadStlMesh(filename);
    else
        data = stlread(filename);
        if isa(data,'triangulation')
            F = data.ConnectivityList;
            V = data.Points;
        else
            F = data.faces;
            V = data.vertices;
        end
    end
end