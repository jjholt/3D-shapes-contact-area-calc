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
    %  unless projection wanted along normal; here we use pure distance.)
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