% main.m
% Time-dependent cube–cylinder contact demo using computeContactArea_STS

clear;
clc;

%% ====== FIGURE OUT PROJECT ROOT & PATHS ======
scriptDir = fileparts(mfilename('fullpath'));

% If main.m is in project root, model folder is here; if in src/, go up one level
if exist(fullfile(scriptDir, 'model'), 'dir')
    projectRoot = scriptDir;
else
    projectRoot = fileparts(scriptDir);  % go up one level
end

% Add src (helpers, transforms, etc.) to path
addpath(genpath(fullfile(projectRoot, 'src')));

% Model directory at project root
modelDir = fullfile(projectRoot, 'model');

%% ================= USER SETTINGS =======================
cubeFile     = fullfile(modelDir, '50mm cube 4 iterations.stl');
cylinderFile = fullfile(modelDir, 'cylinder 30mm diameter 4 it.stl');

tol = 0.25;  % contact tolerance (same units as STL, e.g. mm)

opts.sampleThreshold   = 0.2;
opts.neighRadiusFactor = 5.0;
opts.maxNeighbours     = 30;
opts.roiExpandFactor   = 1.5;

N         = 50;   % number of time steps
maxShift  = 5.0;  % mm, small slide left–right relative to base pose

%% =============== LOAD & ALIGN REFERENCE MESHES =================
[Fc, Vc] = loadStlMesh(cubeFile, 'cube');
[Fs, Vs] = loadStlMesh(cylinderFile, 'cylinder');

% Align cylinder ONCE in the base (reference) pose, same way as static script
[Vs_aligned, cylShift] = placeCylinderOnCube(Vs, Vc);

cubeRef_V     = Vc;           % reference cube vertices
cylinderRef_V = Vs_aligned;   % reference cylinder vertices (already on cube)

cubeRef     = buildBodyStruct(Fc, cubeRef_V);
cylinderRef = buildBodyStruct(Fs, cylinderRef_V);

fprintf('Cube:    %d verts, %d faces\n', size(Fc,1), size(cubeRef_V,1));
fprintf('Cylinder (aligned): %d verts, %d faces\n', size(Fs,1), size(cylinderRef_V,1));

%% =============== DEFINE TIME-DEPENDENT TRANSFORMS ==============
% T_cube(:,:,k) and T_cyl(:,:,k) are world_from_body transforms
T_cube = repmat(eye(4), 1, 1, N);   % cube fixed

T_cyl  = repmat(eye(4), 1, 1, N);   % start as identity for all frames

for k = 1:N
    alpha = (k-1)/(N-1);                 % 0 -> 1
    dx    = (alpha - 0.5)*2*maxShift;    % -maxShift .. +maxShift

    T = eye(4);
    T(1,4) = dx;                         % RELATIVE slide along X only

    T_cyl(:,:,k) = T;
end

%% =============== VISUALISATION SETUP ============================
contactAreas = zeros(N,1);

% Figure 1: full geometry (cube + cylinder, contact highlighted)
figFull = figure('Color','w');
axFull  = gca; hold(axFull,'on'); axis(axFull,'equal');
xlabel(axFull,'X'); ylabel(axFull,'Y'); zlabel(axFull,'Z');
title(axFull,'Cube–cylinder contact over time');
view(axFull,3); camlight(axFull); lighting(axFull,'gouraud');
grid(axFull,'on');

hCube = patch('Faces', Fc, 'Vertices', cubeRef_V, ...
              'FaceColor', [0.8 0.8 0.8], ...
              'EdgeColor', 'none', ...
              'FaceAlpha', 0.3, ...
              'Parent', axFull);

hCyl  = patch('Faces', Fs, 'Vertices', cylinderRef_V, ...
              'FaceColor', [0.2 0.2 1.0], ...
              'EdgeColor', 'k', ...
              'FaceAlpha', 0.9, ...
              'Parent', axFull);

legend(axFull, {'Cube (master)','Cylinder (slave)'});

% Figure 2: contact patch only (cylinder)
figPatch = figure('Color','w');
axPatch  = gca; hold(axPatch,'on'); axis(axPatch,'equal');
xlabel(axPatch,'X'); ylabel(axPatch,'Y'); zlabel(axPatch,'Z');
title(axPatch, 'Contact patch on cylinder');
view(axPatch,3); camlight(axPatch); lighting(axPatch,'gouraud');
grid(axPatch,'on');

hContact = patch('Faces', [], 'Vertices', [], ...
                 'FaceColor', [1.0 0.2 0.2], ...
                 'EdgeColor', 'k', ...
                 'FaceAlpha', 1.0, ...
                 'Parent', axPatch);

%% =============== TIME LOOP: CONTACT PER FRAME ===================
for k = 1:N
    % --- 1) Transform vertices for this frame (relative to base pose) ---
    Vc_k = transformVertices(cubeRef_V,     T_cube(:,:,k));     % here: identity
    Vs_k = transformVertices(cylinderRef_V, T_cyl(:,:,k));      % sliding cyl

    % --- 2) Rebuild body structs from transformed vertices ---
    cube_k     = buildBodyStruct(Fc, Vc_k);
    cylinder_k = buildBodyStruct(Fs, Vs_k);

    % Choose master/slave: cube = master, cylinder = slave
    master = cube_k;
    slave  = cylinder_k;

    % Build KD-tree on master centroids (same as static STS code)
    master.kdtree = KDTreeSearcher(master.triCentroid);

    % --- 3) Compute contact using STS solver ---
    [A_contact, contactMask] = computeContactArea_STS(slave, master, tol, opts);
    contactAreas(k) = A_contact;

    fprintf('Frame %d/%d: contact area = %.6f, contact tris (slave) = %d\n', ...
        k, N, A_contact, nnz(contactMask));

    % --- 4) UPDATE FIGURE 1: cube + cylinder with contact highlighted ---
    if ~isvalid(hCube) || ~isvalid(hCyl)
        warning('Figure closed by user, stopping animation.');
        break;
    end

    set(hCube, 'Vertices', Vc_k, 'Faces', Fc);

    nFacesSlave   = size(slave.F,1);
    contactColors = repmat([0.2 0.2 1.0], nFacesSlave, 1);                 % blue
    contactColors(contactMask,:) = repmat([1.0 0.2 0.2], nnz(contactMask), 1); % red

    set(hCyl, 'Vertices', Vs_k, ...
              'Faces',    Fs, ...
              'FaceVertexCData', contactColors, ...
              'FaceColor', 'flat');

    title(axFull, sprintf('Frame %d/%d, contact area = %.6f', ...
                          k, N, A_contact));

    % --- 5) UPDATE FIGURE 2: contact patch only on cylinder ---
    if any(contactMask)
        F_contact = slave.F(contactMask,:);
        set(hContact, 'Faces', F_contact, 'Vertices', Vs_k);
        title(axPatch, sprintf('Contact patch (frame %d/%d, A = %.6f)', ...
                               k, N, A_contact));
    else
        set(hContact, 'Faces', [], 'Vertices', []);
        title(axPatch, sprintf('Frame %d/%d, NO contact (tol = %.3g)', ...
                               k, N, tol));
    end

    drawnow;
    % pause(0.05); % uncomment to slow down animation
end

%% =============== PLOT CONTACT AREA VS FRAME ======================
figure;
plot(1:numel(contactAreas), contactAreas, '-o');
xlabel('Frame'); ylabel('Contact area');
title('Contact area over time');
grid on;


%% ====== HELPER: APPLY 4x4 TRANSFORM TO VERTEX ARRAY =============
function V_out = transformVertices(V_in, T)
    % V_in: Nx3, vertices in local/base frame
    % T   : 4x4, homogeneous transform (world_from_local)
    Vh = [V_in, ones(size(V_in,1),1)];  % Nx4
    Vw = (T * Vh.').';                  % Nx4
    V_out = Vw(:,1:3);
end


