function showContactOnly(slave, contactMask, tol)
% showContactOnly
% ----------------
% Plots only the contact triangles on the slave surface in a new figure.

if nargin < 3
    tol = [];
end

if ~any(contactMask)
    warning('No contact triangles detected (contactMask is all false).');
end

figure('Color','w');
hold on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');

titleStr = 'Contact patch (slave surface only)';
if ~isempty(tol)
    titleStr = sprintf('%s, tol = %.3g', titleStr, tol);
end
title(titleStr, 'Interpreter','none');

% Extract only the faces in contact
F_contact = slave.F(contactMask,:);
V_slave   = slave.V;

if ~isempty(F_contact)
    patch('Faces',F_contact,'Vertices',V_slave, ...
          'FaceColor',[1.0 0.2 0.2], ...  % red
          'EdgeColor','k', ...
          'FaceAlpha',1.0);
end

view(3); camlight; lighting gouraud;
grid on;

end
