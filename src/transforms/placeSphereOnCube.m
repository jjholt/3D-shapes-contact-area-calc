function [V_sphere_out, shift] = placeSphereOnCube(V_sphere_in, V_cube)
% placeSphereOnCube
% -----------------
% Translate a sphere mesh so that:
%   1) Its projected centre in XY is aligned with the cube centre in XY.
%   2) Its lowest point in Z just touches the cube's top face.
%
% INPUTS:
%   V_sphere_in : Ns x 3 vertices of the sphere (arbitrary pose)
%   V_cube      : Nc x 3 vertices of the cube
%
% OUTPUTS:
%   V_sphere_out : Ns x 3 transformed sphere vertices
%   shift        : 1 x 3 translation applied

    % Cube centre and top Z
    cubeCenter = mean(V_cube, 1);        % [xc, yc, zc]
    cubeTopZ   = max(V_cube(:,3));       % highest z of cube

    % Sphere centre and bottom Z
    sphCenter = mean(V_sphere_in, 1);
    sphBottomZ = min(V_sphere_in(:,3));

    % We want:
    %   - sphere centre XY -> cube centre XY
    %   - sphere bottom Z  -> cube top Z
    dx = cubeCenter(1) - sphCenter(1);
    dy = cubeCenter(2) - sphCenter(2);
    dz = cubeTopZ      - sphBottomZ;

    shift = [dx, dy, dz];

    V_sphere_out = V_sphere_in + shift;
end
