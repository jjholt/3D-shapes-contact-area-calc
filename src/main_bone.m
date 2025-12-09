function main_bone()
% Use optical transform data as ground-truth poses for femur, tibia,
% and patella, optionally re-orient the whole knee, then run contact.

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

    fprintf('Loading STL meshes...\n');
    [F_fem, V_fem] = loadAnyStl(fullfile(modelDir, 'Femur.STL'));
    [F_tib, V_tib] = loadAnyStl(fullfile(modelDir, 'Tibia.STL'));
    [F_pat, V_pat] = loadAnyStl(fullfile(modelDir, 'Patella.STL'));

    % Load optical transforms
    tfFile = fullfile(dataDir, 'transforms_optical-data.mat');
    fprintf('Loading transform data from %s\n', tfFile);
    S = load(tfFile);

    % Assuming structure: transforms.Intact.in_global.(femur/tibia/patella)
    Tf = S.transforms.Intact.in_global.femur;
    Tt = S.transforms.Intact.in_global.tibia;
    Tp = S.transforms.Intact.in_global.patella;

    nFrames = min([size(Tf,3), size(Tt,3), size(Tp,3)]);
    fprintf('Found %d frames of transform data.\n', nFrames);

    % Choose frame index
    frameIdx = round(nFrames/2);
    fprintf('Using frame %d as ground truth.\n', frameIdx);

    T_fem_GT = Tf(:,:,frameIdx);
    T_tib_GT = Tt(:,:,frameIdx);
    T_pat_GT = Tp(:,:,frameIdx);

    % Apply transforms to meshes
    V_fem_world = applySingleTransform(V_fem, T_fem_GT);
    V_tib_world = applySingleTransform(V_tib, T_tib_GT);
    V_pat_world = applySingleTransform(V_pat, T_pat_GT);

    % Re-orientation using tibia
    useReorientation = true;   % set false to stay in optical global

    if useReorientation
        [V_fem_plot, V_tib_plot, V_pat_plot] = ...
            reorientKneeWithTibia(V_fem_world, V_tib_world, V_pat_world);
        coordLabel = sprintf('optical frame %d (re-oriented)', frameIdx);
    else
        V_fem_plot = V_fem_world;
        V_tib_plot = V_tib_world;
        V_pat_plot = V_pat_world;
        coordLabel = sprintf('optical frame %d (raw global)', frameIdx);
    end

    % Visualise posed knee
    figure; hold on; axis equal; grid on;
    patch('Faces', F_tib, 'Vertices', V_tib_plot, ...
          'FaceColor',[0.8 0.8 0.9],'EdgeColor','none','FaceAlpha',0.7);
    patch('Faces', F_fem, 'Vertices', V_fem_plot, ...
          'FaceColor',[0.9 0.7 0.7],'EdgeColor','none','FaceAlpha',0.7);
    patch('Faces', F_pat, 'Vertices', V_pat_plot, ...
          'FaceColor',[0.7 0.9 0.7],'EdgeColor','none','FaceAlpha',0.7);

    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(['Knee at ' coordLabel]);
    view(3); camlight; lighting gouraud; rotate3d on;

    % Build body structs and compute contact
    femurBody   = buildBodyStruct(F_fem, V_fem_plot,   'femur');
    tibiaBody   = buildBodyStruct(F_tib, V_tib_plot,   'tibia');
    patellaBody = buildBodyStruct(F_pat, V_pat_plot,   'patella');

    % Make sure masters have KD-trees
    tibiaBody   = addKDTreeToBody(tibiaBody);
    patellaBody = addKDTreeToBody(patellaBody);

    tol  = 0.25;   % in STL units
    opts.sampleThreshold   = 0.2;
    opts.neighRadiusFactor = 5.0;
    opts.maxNeighbours     = 50;

    fprintf('Computing contact areas (femur–tibia and femur–patella)...\n');
    [A_FT, mask_FT] = computeContactArea_STS(femurBody, tibiaBody,   tol, opts); %#ok<NASGU>
    [A_FP, mask_FP] = computeContactArea_STS(femurBody, patellaBody, tol, opts); %#ok<NASGU>

    fprintf('Femur–Tibia  contact area: %.3f\n', A_FT);
    fprintf('Femur–Patella contact area: %.3f\n', A_FP);
end


function V_out = applySingleTransform(V_in, T)
    N  = size(V_in, 1);
    Vh = [V_in, ones(N,1)];      % N x 4
    Vt = (T * Vh.').';           % N x 4
    V_out = Vt(:,1:3);
end

function [V_fem_out, V_tib_out, V_pat_out] = ...
         reorientKneeWithTibia(V_fem, V_tib, V_pat)

    allVerts = [V_fem; V_tib; V_pat];
    centAll  = mean(allVerts, 1);

    V_fem = V_fem - centAll;
    V_tib = V_tib - centAll;
    V_pat = V_pat - centAll;

    C        = cov(V_tib);
    [E, D]   = eig(C);
    [~, idx] = sort(diag(D), 'descend');
    E        = E(:, idx);

    shaftAxis = E(:,1);
    other1    = E(:,2);
    other2    = E(:,3);

    proj1 = V_tib * other1;
    proj2 = V_tib * other2;
    if (max(proj1) - min(proj1)) >= (max(proj2) - min(proj2))
        yAxis = other1;
        xAxis = cross(yAxis, shaftAxis);
    else
        yAxis = other2;
        xAxis = cross(yAxis, shaftAxis);
    end

    xAxis = xAxis / norm(xAxis);
    yAxis = yAxis / norm(yAxis);
    zAxis = shaftAxis / norm(shaftAxis);

    R = [xAxis, yAxis, zAxis];
    if det(R) < 0
        xAxis = -xAxis;
        R     = [xAxis, yAxis, zAxis];
    end

    V_fem = V_fem * R;
    V_tib = V_tib * R;
    V_pat = V_pat * R;

    tibZ     = V_tib(:,3);
    zPlateau = max(tibZ);

    V_fem(:,3) = V_fem(:,3) - zPlateau;
    V_tib(:,3) = V_tib(:,3) - zPlateau;
    V_pat(:,3) = V_pat(:,3) - zPlateau;

    V_fem_out = V_fem;
    V_tib_out = V_tib;
    V_pat_out = V_pat;
end

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

function body = addKDTreeToBody(body)
    if isfield(body, 'kdtree')
        return;
    end

    F = body.F;
    V = body.V;
    C = ( V(F(:,1),:) + V(F(:,2),:) + V(F(:,3),:) ) / 3;
    body.kdtree = KDTreeSearcher(C);
end