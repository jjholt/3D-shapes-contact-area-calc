% main.m
% Entry point for knee contact simulation using optical transform data.
%
% Expected folder structure:
%   3D-shapes-contact-area-calc/
%       main.m
%       model/
%           femur.stl
%           tibia.stl
%           patella.stl
%       data/
%           transforms_optical-data.mat
%       src/
%           helpers, contact code, etc.
%
% Requires on path:
%   - kneeContactFromTransforms.m
%   - loadStlMesh.m
%   - buildBodyStruct.m
%   - computeContactArea_STS.m

clear; clc; close all;

%% ====== DETECT PROJECT ROOT & SET PATHS ======
scriptDir = fileparts(mfilename('fullpath'));

if exist(fullfile(scriptDir, 'model'), 'dir')
    % main.m is in project root
    projectRoot = scriptDir;
else
    % main.m is in src/; project root is its parent
    projectRoot = fileparts(scriptDir);
end

modelDir = fullfile(projectRoot, 'model');
dataDir  = fullfile(projectRoot, 'data');

% Add src (helpers, transforms, contact, etc.)
addpath(genpath(fullfile(projectRoot, 'src')));

fprintf('Project root detected: %s\n', projectRoot);

%% ====== FILE LOCATIONS ======
% Edit these STL names if files are named differently
femurStl   = fullfile(modelDir, 'femur.stl');
tibiaStl   = fullfile(modelDir, 'tibia.stl');
patellaStl = fullfile(modelDir, 'patella.stl');

matFile    = fullfile(dataDir, 'transforms_optical-data.mat');

% Basic existence checks
assert(exist(femurStl,   'file') == 2, 'Cannot find femur STL at %s',   femurStl);
assert(exist(tibiaStl,   'file') == 2, 'Cannot find tibia STL at %s',   tibiaStl);
assert(exist(patellaStl, 'file') == 2, 'Cannot find patella STL at %s', patellaStl);
assert(exist(matFile,    'file') == 2, 'Cannot find MAT file at %s',    matFile);

%% ====== CONTACT PARAMETERS ======
% tol ~= cartilage-gap-ish in whatever units STL is (likely mm)
tol = 1.0;

opts.sampleThreshold   = 0.2;
opts.neighRadiusFactor = 5.0;
opts.maxNeighbours     = 50;
opts.debugPlot         = false;   % set true if contact code supports debug plots

%% ====== RUN CONTACT ANALYSIS PIPELINE ======
fprintf('Running kneeContactFromTransforms...\n');

results = kneeContactFromTransforms( ...
    femurStl, ...
    tibiaStl, ...
    patellaStl, ...
    matFile, ...
    tol, ...
    opts );

A_pf    = results.A_pf;     % patella–femur
A_tf    = results.A_tf;     % tibia–femur
nFrames = results.nFrames;
time    = results.time;     % [] if not present in MAT

fprintf('Completed contact computation over %d frames.\n', nFrames);

%% ====== PLOT CONTACT AREA OVER FRAMES ======
figure('Name','Knee contact area','Color','w'); hold on;

if ~isempty(time) && numel(time) == nFrames
    x = time(:);
    plot(x, A_pf, 'LineWidth', 1.5);
    plot(x, A_tf, 'LineWidth', 1.5);
    xlabel('Time');
else
    x = (1:nFrames).';
    plot(x, A_pf, 'LineWidth', 1.5);
    plot(x, A_tf, 'LineWidth', 1.5);
    xlabel('Frame');
end

ylabel('Contact area');
legend({'Patella–Femur','Tibia–Femur'}, 'Location','best');
grid on;
title('Knee contact area over motion');

%% ====== OPTIONAL: VISUALISE FIRST FRAME GEOMETRY ======
doVizFirstFrame = true;

if doVizFirstFrame
    fprintf('Visualising first frame geometry...\n');

    % Reload meshes
    [F_femur,   V_femur]   = loadStlMesh(femurStl,   'femur');
    [F_tibia,   V_tibia]   = loadStlMesh(tibiaStl,   'tibia');
    [F_patella, V_patella] = loadStlMesh(patellaStl, 'patella');

    % Load transforms and use Intact.in_global struct
    S     = load(matFile);
    Tcond = S.transforms.Intact;
    G     = Tcond.in_global;   % struct with fields femur, tibia, patella

    % 4x4xN series for each bone
    T_femur   = G.femur;
    T_tibia   = G.tibia;
    T_patella = G.patella;

    % Apply transforms to get vertex positions for all frames
    V_femur_frames   = applyTransformSeries(V_femur,   T_femur);
    V_tibia_frames   = applyTransformSeries(V_tibia,   T_tibia);
    V_patella_frames = applyTransformSeries(V_patella, T_patella);

    % Visualise first frame (k = 1)
    k = 1;
    Vf1 = V_femur_frames(:, :, k);
    Vt1 = V_tibia_frames(:, :, k);
    Vp1 = V_patella_frames(:, :, k);

    figure('Name','Knee geometry – first frame','Color','w'); hold on;

    patch('Faces', F_femur,   'Vertices', Vf1, ...
          'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    patch('Faces', F_tibia,   'Vertices', Vt1, ...
          'FaceColor', [0.5 0.8 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.8);
    patch('Faces', F_patella, 'Vertices', Vp1, ...
          'FaceColor', [0.9 0.5 0.5], 'EdgeColor', 'none', 'FaceAlpha', 0.9);

    axis equal; grid on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(sprintf('Femur, Tibia, Patella – Frame %d', k));
    view(3); camlight headlight; lighting gouraud;
end

%% ====== LOCAL HELPER: applyTransformSeries ============================
function V_frames = applyTransformSeries(V, T)
%APPLYTRANSFORMSERIES Apply a 4x4xN transform series to vertices.
%   V: nVerts x 3
%   T: 4 x 4 x N
%   V_frames: nVerts x 3 x N

    nV    = size(V, 1);
    Vh    = [V, ones(nV, 1)]';      % 4 x nVerts
    Vh_fr = pagemtimes(T, Vh);      % 4 x nVerts x N
    V_frames = permute(Vh_fr(1:3, :, :), [2 1 3]);  % nVerts x 3 x N
end
