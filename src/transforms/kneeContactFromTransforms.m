function results = kneeContactFromTransforms(femurStl, tibiaStl, patellaStl, matFile, tol, opts)
%KNEECONTACTFROMTRANSFORMS
%   Use optical transform data (transforms.Intact.in_global) to move bone
%   STL models (femur, tibia, patella) and compute contact area over time.

    %% --------- Defaults ----------
    if nargin < 5 || isempty(tol)
        tol = 1.0;
    end
    if nargin < 6 || isempty(opts)
        opts.sampleThreshold   = 0.2;
        opts.neighRadiusFactor = 5.0;
        opts.maxNeighbours     = 50;
        opts.debugPlot         = false;
    end

    %% --------- Load STL models FIRST ----------
    fprintf('Loading STL meshes...\n');
    [F_femur,   V_femur]   = loadStlMesh(femurStl,   'femur');
    [F_tibia,   V_tibia]   = loadStlMesh(tibiaStl,   'tibia');
    [F_patella, V_patella] = loadStlMesh(patellaStl, 'patella');

    %% --------- Load transform data ----------
    fprintf('Loading transform data from %s\n', matFile);
    S = load(matFile);

    if ~isfield(S, 'transforms')
        error('MAT file does not contain a ''transforms'' struct at the top level.');
    end
    if ~isfield(S.transforms, 'Intact')
        error('Expected S.transforms.Intact to exist. Check MAT file structure.');
    end

    Tcond = S.transforms.Intact;

    if ~isfield(Tcond, 'in_global')
        error('Expected transforms.Intact.in_global. Check MAT file structure.');
    end

    G = Tcond.in_global;   % struct with fields: femur, tibia, patella

    if ~all(isfield(G, {'femur','tibia','patella'}))
        error('in_global must have fields femur, tibia, patella. Found: %s', ...
              strjoin(fieldnames(G), ', '));
    end

    T_femur   = G.femur;    % 4 x 4 x N
    T_tibia   = G.tibia;    % 4 x 4 x N
    T_patella = G.patella;  % 4 x 4 x N

    % Sanity checks
    assert(ndims(T_femur)   == 3 && size(T_femur,1)   == 4 && size(T_femur,2)   == 4, 'T_femur must be 4x4xN');
    assert(ndims(T_tibia)   == 3 && size(T_tibia,1)   == 4 && size(T_tibia,2)   == 4, 'T_tibia must be 4x4xN');
    assert(ndims(T_patella) == 3 && size(T_patella,1) == 4 && size(T_patella,2) == 4, 'T_patella must be 4x4xN');

    nFrames = size(T_femur, 3);
    if size(T_tibia,3) ~= nFrames || size(T_patella,3) ~= nFrames
        error('Femur, tibia and patella transforms must have same number of frames.');
    end

    fprintf('Found %d frames of transform data.\n', nFrames);

    % Optional time vector if present
    time = [];
    if isfield(Tcond, 'time')
        time = Tcond.time;
    elseif isfield(Tcond, 't')
        time = Tcond.t;
    end

    %% --------- Apply transforms to all frames ----------
    fprintf('Applying transform series to bone meshes...\n');
    V_femur_frames   = applyTransformSeries(V_femur,   T_femur);
    V_tibia_frames   = applyTransformSeries(V_tibia,   T_tibia);
    V_patella_frames = applyTransformSeries(V_patella, T_patella);

    %% --------- Allocate outputs ----------
    A_pf = zeros(nFrames, 1);  % patella–femur
    A_tf = zeros(nFrames, 1);  % tibia–femur

    %% --------- Main frame loop ----------
    fprintf('Computing contact area over %d frames...\n', nFrames);
    for k = 1:nFrames
    V_femur_k   = V_femur_frames(:, :, k);
    V_patella_k = V_patella_frames(:, :, k);
    V_tibia_k   = V_tibia_frames(:, :, k);

    % Build body structs from faces + transformed vertices
    femur_k   = buildBodyStruct(F_femur,   V_femur_k);
    patella_k = buildBodyStruct(F_patella, V_patella_k);
    tibia_k   = buildBodyStruct(F_tibia,   V_tibia_k);

    % ----- Ensure femur_k has a KD-tree for neighbour search -----
    if ~isfield(femur_k, 'kdtree')
      % Try to guess the point set field used in code
    if isfield(femur_k, 'verts')
        masterPts = femur_k.verts;
      elseif isfield(femur_k, 'V')
         masterPts = femur_k.V;
      else
            % Fall back to raw vertices V_femur_k
         masterPts = V_femur_k;
     end

     % Build KD-tree on master points (Statistics & ML toolbox)
      femur_k.kdtree = KDTreeSearcher(masterPts);
    end

% Patellofemoral contact: patella (slave) vs femur (master)
[A_pf(k), ~] = computeContactArea_STS(patella_k, femur_k, tol, opts);

% Tibiofemoral contact: tibia (slave) vs femur (master)
[A_tf(k), ~] = computeContactArea_STS(tibia_k, femur_k, tol, opts);

        if mod(k, max(1, round(nFrames/10))) == 0 || k == 1 || k == nFrames
            fprintf('  Frame %4d/%4d: A_pf = %.2f, A_tf = %.2f\n', ...
                    k, nFrames, A_pf(k), A_tf(k));
        end
    end

    %% --------- Pack results ----------
    results = struct();
    results.A_pf       = A_pf;
    results.A_tf       = A_tf;
    results.nFrames    = nFrames;
    results.time       = time;
    results.tol        = tol;
    results.opts       = opts;
    results.transforms = S.transforms;
end

% ---------- local helper ----------
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
