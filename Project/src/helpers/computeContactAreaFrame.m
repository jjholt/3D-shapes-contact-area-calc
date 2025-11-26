function [contactArea, contactMaskSlave, contactMaskMaster] = ...
         computeContactAreaFrame(masterBody, slaveBody, tol, opts)
% computeContactAreaFrame
% -----------------------
% Compute contact area between two already-transformed meshes
% (masterBody and slaveBody) for a single time frame.
%
% INPUTS:
%   masterBody : struct with fields at least:
%                - V (Nm x 3) vertices
%                - F (Mm x 3) faces (1-based indices)
%
%   slaveBody  : struct with fields at least:
%                - V (Ns x 3) vertices
%                - F (Ks x 3) faces (1-based indices)
%
%   tol        : contact tolerance (same units as V, e.g. mm)
%   opts       : options struct with fields:
%                - sampleThreshold   (fraction of sample points
%                                     within tol to call a slave
%                                     triangle "in contact")
%                - neighRadiusFactor (multiplier for tol when
%                                     defining neighbor search radius)
%
% OUTPUTS:
%   contactArea       : scalar, sum of areas of slave faces in contact
%   contactMaskSlave  : (Ns x 1) logical, which slave faces are in contact
%   contactMaskMaster : (Nm x 1) logical, which master faces are in contact

    % ---------- Options with defaults ----------
    if ~isfield(opts, 'sampleThreshold')
        opts.sampleThreshold = 0.5;
    end
    if ~isfield(opts, 'neighRadiusFactor')
        opts.neighRadiusFactor = 5.0;
    end

    sampleThreshold = opts.sampleThreshold;
    neighRadius     = opts.neighRadiusFactor * tol;

    % ---------- Triangle geometry (master) ----------
    Fm = masterBody.F;
    Vm = masterBody.V;

    v1m = Vm(Fm(:,1),:);
    v2m = Vm(Fm(:,2),:);
    v3m = Vm(Fm(:,3),:);

    centroidsM = (v1m + v2m + v3m) / 3;
    % areasM not strictly needed unless you want to use it later
    % e1m = v2m - v1m;
    % e2m = v3m - v1m;
    % nM  = cross(e1m, e2m, 2);
    % areasM = 0.5 * vecnorm(nM, 2, 2);

    % ---------- Triangle geometry (slave) ----------
    Fs = slaveBody.F;
    Vs = slaveBody.V;

    v1s = Vs(Fs(:,1),:);
    v2s = Vs(Fs(:,2),:);
    v3s = Vs(Fs(:,3),:);

    centroidsS = (v1s + v2s + v3s) / 3;

    e1s = v2s - v1s;
    e2s = v3s - v1s;
    nS  = cross(e1s, e2s, 2);
    areasS = 0.5 * vecnorm(nS, 2, 2);  % slave triangle areas

    nSlave  = size(Fs, 1);
    nMaster = size(Fm, 1);

    contactMaskSlave  = false(nSlave, 1);
    contactMaskMaster = false(nMaster, 1);

    % ---------- KD-tree on master centroids ----------
    kdtree = KDTreeSearcher(centroidsM);

    % neighbour indices for each slave centroid
    idxNeighbors = rangesearch(kdtree, centroidsS, neighRadius);

    % ---------- Loop over slave triangles ----------
    for i = 1:nSlave
        masterIdx = idxNeighbors{i};

        if isempty(masterIdx)
            continue;  % no nearby master triangles -> no contact
        end

        % Sample points on slave triangle: 3 vertices + centroid
        triVerts  = Vs(Fs(i, :), :);     % 3 x 3
        centroid  = centroidsS(i, :);    % 1 x 3
        samplePts = [triVerts; centroid]; % 4 x 3
        nSamples  = 4;

        nWithinTol = 0;

        % Loop over candidate master triangles
        for j = 1:numel(masterIdx)
            mID = masterIdx(j);

            % Master triangle vertices
            mVerts = Vm(Fm(mID, :), :);   % 3 x 3
            a = mVerts(1,:);
            b = mVerts(2,:);
            c = mVerts(3,:);

            % Check distances from each sample point to this triangle
            for s = 1:nSamples
                d = pointTriangleDistance(samplePts(s, :), a, b, c);  % <-- 4 args

                if d <= tol
                    nWithinTol = nWithinTol + 1;
                    contactMaskMaster(mID) = true;
                end
            end
        end

        fracInTol = nWithinTol / nSamples;

        if fracInTol >= sampleThreshold
            contactMaskSlave(i) = true;
        end
    end

    % ---------- Contact area (slave side) ----------
    contactArea = sum(areasS(contactMaskSlave));

end
