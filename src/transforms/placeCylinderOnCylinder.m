function [V_top_out, shift] = placeCylinderOnCylinder(V_top_in, V_bottom)
% placeCylinderOnCylinder
% -----------------------
% Translate a "top" cylinder mesh so that:
%   1) Its projected centre in XY is aligned with the bottom cylinder centre in XY.
%   2) Its lowest point in Z just touches the highest point of the bottom cylinder.
%
% Assumes both cylinders are upright with their axes along Z.

    % Bottom cylinder centre and top Z
    bottomCenter = mean(V_bottom, 1);        % [xb, yb, zb]
    bottomTopZ   = max(V_bottom(:,3));       % highest z of bottom cylinder

    % Top cylinder centre and bottom Z
    topCenter  = mean(V_top_in, 1);
    topBottomZ = min(V_top_in(:,3));

    % We want:
    %   - top centre XY -> bottom centre XY
    %   - top bottom Z  -> bottom top Z
    dx = bottomCenter(1) - topCenter(1);
    dy = bottomCenter(2) - topCenter(2);
    dz = bottomTopZ     - topBottomZ;

    shift = [dx, dy, dz];

    V_top_out = V_top_in + shift;
end
