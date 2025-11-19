function ContactArea_STL_STS()
% ContactArea_STL_STS
% 
% General surface-to-surface contact area estimator between two STL meshes.
%
% - Uses structs to organise mesh data (faces, vertices, centroids, normals, areas).
% - Uses a KD-tree on master triangle centroids for fast neighbour search.
% - Uses a bounding-box Region Of Interest (ROI) to cull triangles far from contact.
% - Uses a surface-to-surface style rule:
%       For each slave triangle, sample 3 vertices + centroid,
%       find nearby master triangles,
%       compute exact point-to-triangle distance,
%       mark triangle as "in contact" if enough samples are within tolerance.
%
% USER SETTINGS:
%   - Adjust STL filenames, tolerance, and options in the block below.

clc; close all;

%% ====================== USER SETTINGS ============================
% STL files
bodyAFile = 'sphere 50mm diameter 3 it.stl';    
bodyBFile = 'sphere cap opening angle 120 degrees 3 it.stl';  

% Contact tolerance (same units as STL, e.g. mm)
tol = 0.5;

% Options for the algorithm, variables
opts.sampleThreshold   = 0.5;    % fraction of samples within tol to call triangle "in contact"
opts.neighRadiusFactor = 5.0;    % neighbour search radius = factor * tol
opts.maxNeighbours     = 30;     % max number of candidate master tris per sample point
opts.roiExpandFactor   = 1.5;    % expand ROI box by this * tol

% Visualisation toggle
doVisualise = true;

% Debug / tester toggles
doDebugPrints = true;   % turn on summary prints

%% ====================== LOAD STL MESHES ==========================
fprintf('Loading STL meshes...\n');
[FA, VA] = loadStlMesh(bodyAFile, 'bodyA');  % assume this is the full sphere
[FB, VB] = loadStlMesh(bodyBFile, 'bodyB');  % assume this is the spherical cap

fprintf('Body A: %d vertices, %d faces\n', size(VA,1), size(FA,1));
fprintf('Body B: %d vertices, %d faces\n', size(VB,1), size(FB,1));

%% ========== ALIGN SPHERICAL CAP CENTRE TO SPHERE CENTRE ==========
% Fit best-fit sphere centres to each mesh
centerA = fitSphereCenter(VA);   % sphere centre (body A)
centerB = fitSphereCenter(VB);   % cap's spherical centre (body B)

% Translate cap so its centre matches the sphere's centre
shift = centerA - centerB;
VB = VB + shift;

if doDebugPrints
    % Re-check cap centre after alignment
    centerB_after = fitSphereCenter(VB);
    fprintf('Centre alignment check: |centerA - centerB| = %.6f\n', ...
            norm(centerA - centerB_after));
end

%% =================== BUILD BODY STRUCTS ==========================
bodyA = buildBodyStruct(FA, VA);
bodyB = buildBodyStruct(FB, VB);

fprintf('Body A total surface area: %.3f\n', sum(bodyA.triArea));
fprintf('Body B total surface area: %.3f\n', sum(bodyB.triArea));

if doDebugPrints
    fprintf('Body A: mean tri area = %.6f, min = %.6f, max = %.6f\n', ...
        mean(bodyA.triArea), min(bodyA.triArea), max(bodyA.triArea));
    fprintf('Body B: mean tri area = %.6f, min = %.6f, max = %.6f\n', ...
        mean(bodyB.triArea), min(bodyB.triArea), max(bodyB.triArea));

    % Check normals are unit length on average
    lenA = sqrt(sum(bodyA.triNormal.^2,2));
    lenB = sqrt(sum(bodyB.triNormal.^2,2));
    fprintf('Body A normals: mean |n| = %.6f (should be ~1)\n', mean(lenA));
    fprintf('Body B normals: mean |n| = %.6f (should be ~1)\n', mean(lenB));
end

%% ============= CHOOSE MASTER / SLAVE SURFACES ====================
% Heuristic from Abaqus guidelines:
%   - If one surface is smaller -> slave.
%   - Otherwise, slave = finer mesh (more faces).
%
% Here we base it on total area first, then face count.

areaA = sum(bodyA.triArea);
areaB = sum(bodyB.triArea);

if areaA <= areaB
    slave = bodyA;
    master = bodyB;
    slaveName  = 'Body A (slave)';
    masterName = 'Body B (master)';
else
    slave = bodyB;
    master = bodyA;
    slaveName  = 'Body B (slave)';
    masterName = 'Body A (master)';
end

fprintf('Slave surface: %s\n', slaveName);
fprintf('Master surface: %s\n', masterName);

%% ============= COMPUTE CONTACT AREA (SURFACE-TO-SURFACE) =========
fprintf('Building KD-tree on master centroids...\n');
master.kdtree = KDTreeSearcher(master.triCentroid);

fprintf('Computing contact area...\n');

[contactArea, contactMask] = computeContactArea_STS(slave, master, tol, opts);

fprintf('\nEstimated contact area: %.4f (same units^2 as STL)\n', contactArea);
fprintf('Number of slave triangles in contact: %d\n', nnz(contactMask));

if doDebugPrints
    fracContactTris = nnz(contactMask) / numel(contactMask);
    fprintf('Fraction of slave triangles in contact: %.4f\n', fracContactTris);

    % Rough scale check: contact area should be <= total slave area
    totalSlaveArea = sum(slave.triArea);
    fprintf('Contact area / total slave area = %.6f\n', ...
        contactArea / max(totalSlaveArea, eps));
end

%% ======================= VISUALISATION ===========================
if doVisualise
    figure('Color','w'); hold on; axis equal;
    title(sprintf('Contact patch (tol = %.3f)', tol), 'Interpreter','none');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    
    % Plot master mesh in light grey
    patch('Faces',master.F,'Vertices',master.V,...
          'FaceColor',[0.8 0.8 0.8],'EdgeColor','none','FaceAlpha',0.3);
    
    % Colour slave triangles by contact status
    contactColors = repmat([0.2 0.2 1.0], size(slave.F,1), 1);  % default blue-ish
    contactColors(contactMask,:) = repmat([1.0 0.2 0.2], nnz(contactMask), 1); % red for contact
    
    patch('Faces',slave.F,'Vertices',slave.V,...
          'FaceVertexCData',contactColors,...
          'FaceColor','flat','EdgeColor','k','FaceAlpha',0.9);
    
    legend({'Master','Slave (blue = no contact, red = contact)'});
    view(3); camlight; lighting gouraud;

    % Also show a separate figure with only the contact patch on the slave
    showContactOnly(slave, contactMask, tol);
end

end % main function


%% =================================================================
%                       HELPER FUNCTIONS
%% =================================================================

function [A_contact, contactMask] = computeContactArea_STS(slave, master, tol, opts)
% computeContactArea_STS
% ----------------------
% Surface-to-surface style contact area estimate.
%
% For each slave triangle in a Region Of Interest (bounding box intersection),
% we:
%   - sample its 3 vertices + centroid,
%   - find nearby master triangles via KD-tree,
%   - compute exact point–triangle distances,
%   - classify triangle as "in contact" if fraction of samples within tol
%     exceeds sampleThreshold.
%
% INPUTS:
%   slave, master : structs from buildBodyStruct()
%   tol           : contact tolerance
%   opts          : options struct
%
% OUTPUTS:
%   A_contact     : total contact area on slave surface
%   contactMask   : logical array over slave.F (true = in contact)

% Unpack options with defaults
if ~isfield(opts,'sampleThreshold'),   opts.sampleThreshold   = 0.5; end
if ~isfield(opts,'neighRadiusFactor'), opts.neighRadiusFactor = 5.0; end
if ~isfield(opts,'maxNeighbours'),     opts.maxNeighbours     = 30; end
if ~isfield(opts,'roiExpandFactor'),   opts.roiExpandFactor   = 1.5; end

sampleThreshold   = opts.sampleThreshold;
neighRadius       = opts.neighRadiusFactor * tol;
maxNeighbours     = opts.maxNeighbours;
roiExpandFactor   = opts.roiExpandFactor;

numSlaveTris = size(slave.F,1);
contactMask  = false(numSlaveTris,1);

% ----------------- STEP 1: ROI Bounding Box -----------------------
roiMin = max(slave.bbox(:,1), master.bbox(:,1));
roiMax = min(slave.bbox(:,2), master.bbox(:,2));

% Expand ROI slightly by tol
roiMin = roiMin - roiExpandFactor * tol;
roiMax = roiMax + roiExpandFactor * tol;

% If boxes do not overlap at all, no contact
if any(roiMin > roiMax)
    A_contact = 0.0;
    return;
end

% Select slave triangles whose centroids lie in ROI
C = slave.triCentroid;
inROI = C(:,1) >= roiMin(1) & C(:,1) <= roiMax(1) & ...
        C(:,2) >= roiMin(2) & C(:,2) <= roiMax(2) & ...
        C(:,3) >= roiMin(3) & C(:,3) <= roiMax(3);

candidateIdx = find(inROI);
fprintf('Slave triangles in ROI: %d out of %d\n', numel(candidateIdx), numSlaveTris);

if isempty(candidateIdx)
    A_contact = 0.0;
    return;
end

% ---------------- STEP 2: Loop over candidate slave tris ----------
V_slave = slave.V;
F_slave = slave.F;
C_slave = slave.triCentroid;

for kk = 1:numel(candidateIdx)
    iTri = candidateIdx(kk);

    % Sample points: 3 vertices + centroid
    vertsIdx = F_slave(iTri,:);
    p1 = V_slave(vertsIdx(1),:);
    p2 = V_slave(vertsIdx(2),:);
    p3 = V_slave(vertsIdx(3),:);
    pc = C_slave(iTri,:);

    samplePts = [p1; p2; p3; pc];
    nSamples  = size(samplePts,1);
    inContactSamples = false(nSamples,1);

    % Approximate slave triangle normal (surface-to-surface style)
    % (already stored, but we don't actually need it for distance-only check
    %  unless you want projection along normal; here we use pure distance.)
    % normal_s = slave.triNormal(iTri,:);

    for s = 1:nSamples
        ps = samplePts(s,:);   % 1x3

        % KD-tree neighbour search in centroid space
        [idxCell, distCell] = rangesearch(master.kdtree, ps, neighRadius);

        idxList  = idxCell{1};
        distList = distCell{1};

        if isempty(idxList)
            % No nearby master triangles -> definitely no contact at this point
            continue;
        end

        % Limit number of neighbours to maxNeighbours (closest ones)
        if numel(idxList) > maxNeighbours
            [~, sortIdx] = sort(distList, 'ascend');
            idxList = idxList(sortIdx(1:maxNeighbours));
        end

        % Compute exact point–triangle distances to candidate master tris
        minDist = inf;
        for mIdx = idxList(:)'
            triIdx = master.F(mIdx,:);
            a = master.V(triIdx(1),:);
            b = master.V(triIdx(2),:);
            c = master.V(triIdx(3),:);

            d = pointTriangleDistance(ps, a, b, c);
            if d < minDist
                minDist = d;
            end
            if minDist <= tol
                % Early exit for this sample point
                break;
            end
        end

        if minDist <= tol
            inContactSamples(s) = true;
        end
    end

    % Surface-to-surface style: average over samples
    fracContact = sum(inContactSamples) / nSamples;

    if fracContact >= sampleThreshold
        contactMask(iTri) = true;
    end
end

% ---------------- STEP 3: Sum contact area ------------------------
A_contact = sum(slave.triArea(contactMask));


end 