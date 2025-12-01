function bodyOut = applyTransformToBody(bodyIn, T)
    % bodyIn: struct with at least .V and .F (and possibly .centroids, .normals, etc.)
    % T: 4x4 homogeneous transform (world_from_body)

    V = bodyIn.V;                     % Nx3
    Vh = [V, ones(size(V,1),1)];      % Nx4
    Vw = (T * Vh.').';                % Nx4 -> world
    Vw = Vw(:,1:3);

    bodyOut = bodyIn;
    bodyOut.V = Vw;

    % If your body struct stores normals/centroids, update them too:
    if isfield(bodyIn, 'centroids')
        C = bodyIn.centroids;
        Ch = [C, ones(size(C,1),1)];
        Cw = (T * Ch.').';
        bodyOut.centroids = Cw(:,1:3);
    end

    if isfield(bodyIn, 'normals')
        R = T(1:3,1:3);
        bodyOut.normals = (R * bodyIn.normals.').';
    end
end
