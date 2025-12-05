function runDynamicContact_CylinderOnCylinder(modelDir)
% runDynamicContact_CylinderOnCylinder
% ------------------------------------
% Time-dependent demo: top cylinder sliding over a bottom cylinder,
% computing contact area at each frame.

clc; close all;

%% ================= USER SETTINGS =======================
bottomCylFile = fullfile(modelDir, 'cylinder 30mm diameter 4 it.stl');
topCylFile    = fullfile(modelDir, 'cylinder 30mm diameter 4 it.stl');  % can be same

tol = 0.25;  % contact tolerance (same units as STL, e.g. mm)

opts.sampleThreshold   = 0.2;
opts.neighRadiusFactor = 5.0;
opts.maxNeighbours     = 30;
opts.roiExpandFactor   = 1.5;

N        = 50;   % number of time steps
maxShift = 5.0;  % mm, slide left–right relative to base pose

%% =============== LOAD & ALIGN REFERENCE MESHES =================
[Fbot, Vbot] = loadStlMesh(bottomCylFile, 'cyl_bottom');
[Ftop, Vtop] = loadStlMesh(topCylFile,    'cyl_top');

fprintf('Bottom cylinder: %d faces, %d verts\n', size(Fbot,1), size(Vbot,1));
fprintf('Top cylinder:    %d faces, %d verts\n', size(Ftop,1), size(Vtop,1));

% Align top cylinder ONCE in the base (reference) pose:
%   - centres aligned in X,Y
%   - top of bottom cylinder touches bottom of top cylinder in Z
[Vtop_aligned, topShift] = placeCylinderOnCylinder(Vtop, Vbot);
fprintf('Top cylinder translated by [%.4f  %.4f  %.4f]\n', topShift);

bottomRef_V = Vbot;          % reference bottom cylinder vertices
topRef_V    = Vtop_aligned;  % reference top cylinder vertices

%% =============== DEFINE TIME-DEPENDENT TRANSFORMS ==============
% T_bot(:,:,k) and T_top(:,:,k) are world_from_body transforms
T_bot = repmat(eye(4), 1, 1, N);   % bottom cylinder fixed

T_top = repmat(eye(4), 1, 1, N);   % start as identity for all frames

for k = 1:N
    alpha = (k-1)/(N-1);                 % 0 -> 1
    dx    = (alpha - 0.5)*2*maxShift;    % -maxShift .. +maxShift

    T = eye(4);
    T(1,4) = dx;                         % slide along X

    T_top(:,:,k) = T;
end

%% =============== VISUALISATION SETUP ============================
contactAreas = zeros(N,1);

% Figure 1: full geometry (bottom + top cylinder, contact highlighted)
figFull = figure('Color','w');
axFull  = gca; hold(axFull,'on'); axis(axFull,'equal');
xlabel(axFull,'X'); ylabel(axFull,'Y'); zlabel(axFull,'Z');
title(axFull,'Cylinder–cylinder contact over time');
view(axFull,3); camlight(axFull); lighting(axFull,'gouraud');
grid(axFull,'on');

hBot = patch('Faces', Fbot, 'Vertices', bottomRef_V, ...
             'FaceColor', [0.8 0.8 0.8], ...
             'EdgeColor', 'none', ...
             'FaceAlpha', 0.3, ...
             'Parent', axFull);

hTop = patch('Faces', Ftop, 'Vertices', topRef_V, ...
             'FaceColor', [0.2 0.2 1.0], ...
             'EdgeColor', 'k', ...
             'FaceAlpha', 0.9, ...
             'Parent', axFull);

legend(axFull, {'Bottom cylinder (master)','Top cylinder (slave)'});

% Figure 2: contact patch only (top cylinder)
figPatch = figure('Color','w');
axPatch  = gca; hold(axPatch,'on'); axis(axPatch,'equal');
xlabel(axPatch,'X'); ylabel(axPatch,'Y'); zlabel(axPatch,'Z');
title(axPatch, 'Contact patch on top cylinder');
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
    Vbot_k = transformVertices(bottomRef_V, T_bot(:,:,k));  % bottom fixed
    Vtop_k = transformVertices(topRef_V,    T_top(:,:,k));  % top sliding

    % --- 2) Rebuild body structs from transformed vertices ---
    bot_k = buildBodyStruct(Fbot, Vbot_k);
    top_k = buildBodyStruct(Ftop, Vtop_k);

    % Master/slave choice: bottom = master, top = slave
    master = bot_k;
    slave  = top_k;

    % Build KD-tree on master centroids
    master.kdtree = KDTreeSearcher(master.triCentroid);

    % --- 3) Compute contact using STS solver ---
    [A_contact, contactMask] = computeContactArea_STS(slave, master, tol, opts);
    contactAreas(k) = A_contact;

    fprintf('Frame %d/%d: contact area = %.6f, contact tris (slave) = %d\n', ...
        k, N, A_contact, nnz(contactMask));

    % --- 4) UPDATE FIGURE 1: geometry with contact highlighted ---
    if ~isvalid(hBot) || ~isvalid(hTop)
        warning('Figure closed by user, stopping animation.');
        break;
    end

    set(hBot, 'Vertices', Vbot_k, 'Faces', Fbot);

    nFacesSlave   = size(slave.F,1);
    contactColors = repmat([0.2 0.2 1.0], nFacesSlave, 1);                     % blue
    contactColors(contactMask,:) = repmat([1.0 0.2 0.2], nnz(contactMask), 1); % red

    set(hTop, 'Vertices', Vtop_k, ...
              'Faces',    Ftop, ...
              'FaceVertexCData', contactColors, ...
              'FaceColor', 'flat');

    title(axFull, sprintf('Frame %d/%d, contact area = %.6f', ...
                          k, N, A_contact));

    % --- 5) UPDATE FIGURE 2: contact patch only on top cylinder ---
    if any(contactMask)
        F_contact = slave.F(contactMask,:);
        set(hContact, 'Faces', F_contact, 'Vertices', Vtop_k);
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
figure('Color','w');
plot(1:numel(contactAreas), contactAreas, '-o');
xlabel('Frame'); ylabel('Contact area');
title('Cylinder–cylinder contact area over time');
grid on;

end


%% ====== HELPER: APPLY 4x4 TRANSFORM TO VERTEX ARRAY =============
function V_out = transformVertices(V_in, T)
    % V_in: Nx3, vertices in local/base frame
    % T   : 4x4, homogeneous transform (world_from_local)
    Vh = [V_in, ones(size(V_in,1),1)];  % Nx4
    Vw = (T * Vh.').';                  % Nx4
    V_out = Vw(:,1:3);
end
