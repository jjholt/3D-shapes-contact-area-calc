clear; clc;

scriptDir = fileparts(mfilename('fullpath'));

if exist(fullfile(scriptDir, 'model'), 'dir')
    projectRoot = scriptDir;
else
    projectRoot = fileparts(scriptDir);
end

modelDir = fullfile(projectRoot, 'model');
dataDir  = fullfile(projectRoot, 'data');

addpath(genpath(fullfile(projectRoot, 'src')));

%Load transforms from .mat
tfFile = fullfile(dataDir, 'transforms_optical-data.mat');
S      = load(tfFile);

T_fem_all = S.transforms.Intact.in_global.femur;   % 4x4xN
T_tib_all = S.transforms.Intact.in_global.tibia;   % 4x4xN
T_pat_all = S.transforms.Intact.in_global.patella; % 4x4xN

nFrames = min([size(T_fem_all,3), size(T_tib_all,3), size(T_pat_all,3)]);
fprintf('Loaded transforms, using %d frames.\n', nFrames);

%Plot frames only (no meshes)
figure; hold on; axis equal; grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Optical global frames (femur / tibia / patella)');
view(3);

scale = 10;    % axis length for drawn frames (tweak)
step  = max(1, floor(nFrames / 50));  % plot ~50 timesteps max

for k = 1:step:nFrames
    plotFrame(T_fem_all(:,:,k), scale, [1 0 0]); % femur = red
    plotFrame(T_tib_all(:,:,k), scale, [0 1 0]); % tibia = green
    plotFrame(T_pat_all(:,:,k), scale, [0 0 1]); % patella = blue
end

camlight; lighting gouraud;
title(sprintf('Optical frames only (every %dth frame)', step));

%Load meshes and overlay one frame
[F_fem, V_fem] = loadStlMesh(fullfile(modelDir, 'Femur.STL'));
[F_tib, V_tib] = loadStlMesh(fullfile(modelDir, 'Tibia.STL'));
[F_pat, V_pat] = loadStlMesh(fullfile(modelDir, 'Patella.STL'));

frameIdx = 1;          % choose a frame to inspect
useInverse = false;    % set true to test inv(T) instead

T_fem = T_fem_all(:,:,frameIdx);
T_tib = T_tib_all(:,:,frameIdx);
T_pat = T_pat_all(:,:,frameIdx);

if useInverse
    T_fem = inv(T_fem);
    T_tib = inv(T_tib);
    T_pat = inv(T_pat);
end

V_fem_world = applySingleTransform(V_fem, T_fem);
V_tib_world = applySingleTransform(V_tib, T_tib);
V_pat_world = applySingleTransform(V_pat, T_pat);

figure; hold on; axis equal; grid on;
xlabel('X'); ylabel('Y'); zlabel('Z');
if useInverse
    title(sprintf('Frame %d (using inv(T_*))', frameIdx));
else
    title(sprintf('Frame %d (using T_*)', frameIdx));
end
view(3);

% Plot meshes
patch('Faces', F_tib, 'Vertices', V_tib_world, ...
      'FaceColor',[0.8 0.8 0.9],'EdgeColor','none','FaceAlpha',0.4);
patch('Faces', F_fem, 'Vertices', V_fem_world, ...
      'FaceColor',[0.9 0.7 0.7],'EdgeColor','none','FaceAlpha',0.4);
patch('Faces', F_pat, 'Vertices', V_pat_world, ...
      'FaceColor',[0.7 0.9 0.7],'EdgeColor','none','FaceAlpha',0.4);

% Plot frames at the same time
plotFrame(T_fem, 10, [1 0 0]);
plotFrame(T_tib, 10, [0 1 0]);
plotFrame(T_pat, 10, [0 0 1]);

legend({'Tibia','Femur','Patella'}, 'Location','bestoutside');
camlight; lighting gouraud; rotate3d on;

fprintf('Showing frame %d with useInverse = %d\n', frameIdx, useInverse);


function V_out = applySingleTransform(V_in, T)
    % Apply 4x4 homogeneous transform T to N x 3 vertices
    N  = size(V_in, 1);
    Vh = [V_in, ones(N,1)];   % N x 4
    Vt = (T * Vh.').';        % N x 4
    V_out = Vt(:,1:3);
end

function plotFrame(T, s, col)
    % Plot a coordinate frame given 4x4 T
    % T: 4x4 homogeneous transform
    % s: axis length
    % col: colour for origin point
    R = T(1:3,1:3);
    p = T(1:3,4).';

    % Origin
    plot3(p(1), p(2), p(3), '.', 'Color', col, 'MarkerSize', 15);

    % Axes
    xAxis = p + s * R(:,1).';
    yAxis = p + s * R(:,2).';
    zAxis = p + s * R(:,3).';

    line([p(1) xAxis(1)], [p(2) xAxis(2)], [p(3) xAxis(3)], ...
         'Color', col, 'LineStyle','-');   % local x
    line([p(1) yAxis(1)], [p(2) yAxis(2)], [p(3) yAxis(3)], ...
         'Color', col, 'LineStyle','--');  % local y
    line([p(1) zAxis(1)], [p(2) zAxis(2)], [p(3) zAxis(3)], ...
         'Color', col, 'LineStyle',':');   % local z
end