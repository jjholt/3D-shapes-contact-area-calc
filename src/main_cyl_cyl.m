% main_cyl_on_cyl_static.m
% ------------------------------------------------------------
% Static cylinder-on-cylinder contact using computeContactArea_STS.
% - Loads the same cylinder STL twice (bottom + top).
% - Uses placeCylinderOnCylinder to stack the top cylinder on the bottom.
% - Computes static contact area (no time loop).
% ------------------------------------------------------------

clear; clc; close all;

%% ====== FIGURE OUT PROJECT ROOT & PATHS ======
scriptDir = fileparts(mfilename('fullpath'));

% If this file is in project root, model folder is here; if in src/, go up one level
if exist(fullfile(scriptDir, 'model'), 'dir')
    projectRoot = scriptDir;
else
    projectRoot = fileparts(scriptDir);
end

% Add src (helpers, transforms, etc.) to path
addpath(genpath(fullfile(projectRoot, 'src')));

% Model directory at project root
modelDir = fullfile(projectRoot, 'model');

%% ================= USER SETTINGS =======================
cylFile = fullfile(modelDir, 'cylinder 30mm diameter 4 it.stl');

tol = 0.25;  % contact tolerance (same units as STL)

opts.sampleThreshold   = 0.2;
opts.neighRadiusFactor = 5.0;
opts.maxNeighbours     = 30;
opts.roiExpandFactor   = 1.5;

%% ====================== LOAD STL MESHES ==========================
fprintf('Loading cylinder STL twice (bottom + top)...\n');
[Fbot, Vbot] = loadStlMesh(cylFile, 'bottom');
[Ftop, Vtop] = loadStlMesh(cylFile, 'top');

fprintf('Mesh loaded: %d vertices, %d faces\n', size(Vbot,1), size(Fbot,1));

%% ========== ALIGN TOP CYLINDER ON BOTTOM CYLINDER ================
% Use helper to place top cylinder along axis of bottom cylinder.
% This gives end-to-end "stacked" cylinders, no overlap, no rotation issues.

[Vtop_aligned, shift] = placeCylinderOnCylinder(Vtop, Vbot);
fprintf('Top cylinder translated by [%.4f  %.4f  %.4f]\n', shift);

%% =================== BUILD BODY STRUCTS ==========================
bottom = buildBodyStruct(Fbot, Vbot);
top    = buildBodyStruct(Ftop, Vtop_aligned);

fprintf('Bottom total area : %.3f\n', sum(bottom.triArea));
fprintf('Top total area    : %.3f\n', sum(top.triArea));

%% ============= CHOOSE MASTER / SLAVE SURFACES ====================
% For stacked cylinder–cylinder: bottom = master, top = slave
master = bottom;
slave  = top;

fprintf('Master surface: bottom cylinder\n');
fprintf('Slave surface : top cylinder\n');

%% ============= COMPUTE CONTACT AREA (SURFACE-TO-SURFACE) =========
fprintf('Building KD-tree on master centroids...\n');
master.kdtree = KDTreeSearcher(master.triCentroid);

fprintf('Computing contact area...\n');
[contactArea, contactMask] = computeContactArea_STS(slave, master, tol, opts);

fprintf('\nEstimated contact area: %.4f\n', contactArea);
fprintf('Slave triangles in contact: %d\n', nnz(contactMask));
fprintf('Fraction contact: %.4f\n', nnz(contactMask)/numel(contactMask));

%% ======================= VISUALISATION ===========================
figure('Color','w'); hold on; axis equal;
title(sprintf('Static stacked cylinder-on-cylinder (tol = %.3f)', tol), 'Interpreter','none');
xlabel X; ylabel Y; zlabel Z;

% Plot master mesh in light grey
patch('Faces',master.F,'Vertices',master.V,...
      'FaceColor',[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',0.3);

% Colour slave triangles by contact status
contactColors = repmat([0.2 0.2 1.0], size(slave.F,1), 1); % blue
contactColors(contactMask,:) = repmat([1.0 0.2 0.2], nnz(contactMask), 1); % red

patch('Faces',slave.F,'Vertices',slave.V,...
      'FaceVertexCData',contactColors,'FaceColor','flat',...
      'EdgeColor','k','FaceAlpha',0.9);

legend({'Bottom (master)','Top (slave: red = contact)'});
view(3); camlight; lighting gouraud; grid on;

if exist('showContactOnly','file') == 2
    showContactOnly(slave, contactMask, tol);
end