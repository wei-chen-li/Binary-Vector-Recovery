clear
clc

figure
hold on
axis equal off

Model = myModel;
phi_range = [pi/2 2*pi];

drawPipe(Model.pipe, Model.comsol.truncate_length, phi_range)
drawCoil(Model.coil)
drawSensors(Model.sensors)
drawMesh(Model.mesh, phi_range)
% drawSensitivity(Model)

view(-30,25)
camtarget([0 0 0])
camlight
camva(1.9)



function drawPipe(pipe, pipe_length, phi_range, color, alpha)

if nargin < 3, phi_range = [pi/2 2*pi]; end
if nargin < 4, color = [0.7 0.7 0.7]; end
if nargin < 5, alpha = 1; end

drawCylinder([pipe.inner_diameter/2 pipe.outer_diameter/2], phi_range, [-pipe_length/2 pipe_length/2], color, alpha)

end


function drawCoil(coil, color, alpha)

if nargin < 3, color = [1 0 0]; end
if nargin < 4, alpha = 0.5; end

drawCylinder([coil.inner_radius coil.outer_radius], [0 2*pi], [ coil.spacing/2  coil.spacing/2+coil.height], color, alpha)
drawCylinder([coil.inner_radius coil.outer_radius], [0 2*pi], [-coil.spacing/2 -coil.spacing/2-coil.height], color, alpha)

end


function drawSensors(sensors, color, alpha)

if nargin < 2, color = [0 0 0]; end
if nargin < 3, alpha = 1.0; end

rho = sensors.placement_rho;
phi = sensors.placement_phi;

nodes = table2array(combinations([-1 1] * 2e-4, [-1 1] * 6e-4, [-1 1] * 3e-4))';
faces = [[1 3 7 5]; [2 6 8 4]; [1 5 6 2]; [3 4 8 7]; [1 2 4 3]; [5 7 8 6]];

for i = 1:numel(rho)
    pos = [0 -rho(i)*cos(phi(i)) rho(i)*sin(phi(i))];
    axis = [0 sin(phi(i)) cos(phi(i))];

    p2 = reshape(axis,3,1) / norm(axis);
    p1 = [0; -p2(3); p2(2)]; p1 = p1 / norm(p1);
    p3 = cross(p1, p2);
    T = [p1 p2 p3];
    nodes_ = T * nodes + reshape(pos,3,1);
    
    patch('Faces',faces, 'Vertices',nodes_', ...
      'FaceColor', color, ...
      'FaceAlpha', alpha, ...
      'LineStyle', 'none')
end

end


function drawCylinder(rho_range, phi_range, z_range, color, alpha)

rho1 = rho_range(1); rho2 = rho_range(2);
phi1 = phi_range(1); phi2 = phi_range(2);
z1   = z_range(1);   z2   = z_range(2);

phi = linspace(phi1, phi2);

y_in  = -rho1 * cos(phi);
y_out = -rho2 * cos(phi);
z_in  = rho1 * sin(phi);
z_out = rho2 * sin(phi);
x_bot = -z2 * ones(size(phi));
x_top = -z1 * ones(size(phi));

bot = surf([x_bot;x_bot], [y_out;y_in ], [z_out;z_in ]);
top = surf([x_top;x_top], [y_out;y_in ], [z_out;z_in ]);
out = surf([x_bot;x_top], [y_out;y_out], [z_out;z_out]);
in  = surf([x_bot;x_top], [y_in; y_in ], [z_in; z_in ]);
surfaces = [bot, top, out, in];

if abs(wrapToPi(phi1-phi2)) > 1e-6
    I = ones(1,2);
    n = length(phi);
    surfaces(end+1) = surf([x_bot(1)*I;x_top(1)*I]', [y_out(1)*I;y_in(1)*I], [z_out(1)*I;z_in(1)*I]);
    surfaces(end+1) = surf([x_bot(n)*I;x_top(n)*I]', [y_out(n)*I;y_in(n)*I], [z_out(n)*I;z_in(n)*I]);
end

set(surfaces, 'FaceColor',color, 'FaceAlpha',alpha, 'LineStyle','none');

end


function drawMesh(mesh, phi_range)

if nargin < 2

voxel_corners = mesh.voxel_corners;
[x,y,z] = evalCyl2simCart(voxel_corners(1,:,:), voxel_corners(2,:,:), voxel_corners(3,:,:));
voxel_corners = cat(1, x,y,z);

gauss_nodes = mesh.gauss_nodes;
[x,y,z] = evalCyl2simCart(gauss_nodes(1,:,:), gauss_nodes(2,:,:), gauss_nodes(3,:,:));
gauss_nodes = cat(1, x,y,z);

for n = 1:size(voxel_corners,2)
    x = voxel_corners(1,n,:);
    y = voxel_corners(2,n,:);
    z = voxel_corners(3,n,:);
    plot3([x(1) x(2)], [y(1) y(2)], [z(1) z(2)], 'r')
    plot3([x(1) x(4)], [y(1) y(4)], [z(1) z(4)], 'g')
    plot3([x(1) x(5)], [y(1) y(5)], [z(1) z(5)], 'b')
end
for n = 5:100:size(voxel_corners,2)
for q = 1:size(gauss_nodes,3)
    plot3(gauss_nodes(1,n,q), gauss_nodes(2,n,q), gauss_nodes(3,n,q), 'k.', 'MarkerSize',3)
end
end

else % nargin >= 2

phi_range(1) = phi_range(1) - pi * 1e-3;
phi_range(2) = phi_range(2) + pi * 1e-3;

voxel_corners = mesh.voxel_corners;
rho = unique(voxel_corners(1,:,:));
phi = unique(voxel_corners(2,:,:));
z   = unique(voxel_corners(3,:,:));

phi = phi(phi_range(1) <= phi & phi <= phi_range(2));
z1 = z(1:length(z)/2);
z2 = z(length(z)/2+1:end);

for k = 1:length(z)
    [x_,y_,z_] = evalCyl2simCart([rho(1) rho(end)], [phi(1) phi(1)], [z(k) z(k)]);
    plot3(x_,y_,z_, 'k')
    [x_,y_,z_] = evalCyl2simCart([rho(1) rho(end)], [phi(end) phi(end)], [z(k) z(k)]);
    plot3(x_,y_,z_, 'k')

    phi_fine = linspace(phi(1), phi(end));
    [x_,y_,z_] = evalCyl2simCart(rho(1)*ones(size(phi_fine)), phi_fine, z(k)*ones(size(phi_fine)));
    plot3(x_,y_,z_, 'k')
    [x_,y_,z_] = evalCyl2simCart(rho(end)*ones(size(phi_fine)), phi_fine, z(k)*ones(size(phi_fine)));
    plot3(x_,y_,z_, 'k')
end

for k = 1:length(rho)
    [x_,y_,z_] = evalCyl2simCart([rho(k) rho(k)], [phi(1) phi(1)], [z1(1) z1(end)]);
    plot3(x_,y_,z_, 'k')
    [x_,y_,z_] = evalCyl2simCart([rho(k) rho(k)], [phi(1) phi(1)], [z2(1) z2(end)]);
    plot3(x_,y_,z_, 'k')
    [x_,y_,z_] = evalCyl2simCart([rho(k) rho(k)], [phi(end) phi(end)], [z1(1) z1(end)]);
    plot3(x_,y_,z_, 'k')
    [x_,y_,z_] = evalCyl2simCart([rho(k) rho(k)], [phi(end) phi(end)], [z2(1) z2(end)]);
    plot3(x_,y_,z_, 'k')
end

for k = 1:length(phi)
    [x_,y_,z_] = evalCyl2simCart([rho(1) rho(1)], [phi(k) phi(k)], [z1(1) z1(end)]);
    plot3(x_,y_,z_, 'k')
    [x_,y_,z_] = evalCyl2simCart([rho(1) rho(1)], [phi(k) phi(k)], [z2(1) z2(end)]);
    plot3(x_,y_,z_, 'k')
    [x_,y_,z_] = evalCyl2simCart([rho(end) rho(end)], [phi(k) phi(k)], [z1(1) z1(end)]);
    plot3(x_,y_,z_, 'k')
    [x_,y_,z_] = evalCyl2simCart([rho(end) rho(end)], [phi(k) phi(k)], [z2(1) z2(end)]);
    plot3(x_,y_,z_, 'k')
end

end % nargin

end


function drawSensitivity(Model)

addpath("functions\")

freq = 16000;
mesh_size = 5;

num_sensors = length(Model.sensors.placement_phi);
shift = 15;

rho = (0.95 * Model.pipe.inner_diameter + 0.05 * Model.pipe.outer_diameter) / 2;
phi = linspace(0, 2*pi, num_sensors*shift+1); phi = phi(1:end-1);
z = linspace(-Model.comsol.truncate_length/2, Model.comsol.truncate_length/2, 401);

[phi,z,rho] = ndgrid(phi,z,rho);

Model = ComsolInterface('Build', Model);
ComsolInterface('Solve', Model, freq, mesh_size)
[Ephi1, Ez1] = ComsolInterface('EvaluateField', Model, {'Ephi' 'Ez'}, rho,phi,z);
E1 = cat(3, Ephi1, Ez1);

Model = ComsolInterface('Build', Model, 'sensor1');
ComsolInterface('Solve', Model, freq, mesh_size)
[Ephi2, Ez2] = ComsolInterface('EvaluateField', Model, {'Ephi' 'Ez'}, rho,phi,z);
E2 = cat(3, Ephi2, Ez2);

S_sigma = dot(E1, E2, 3);

C = 0;
for i = 1:num_sensors
    C = C + abs(S_sigma).^2;
    S_sigma = circshift(S_sigma, shift, 1);
end
C = sqrt(C);

[x_,y_,z_] = evalCyl2simCart(rho,phi,z);
x_ = [x_; x_(1,:)];
y_ = [y_; y_(1,:)];
z_ = [z_; z_(1,:)];
C  = [C ; C(1,:) ];

surf(x_,y_,z_,C, 'LineStyle','none', 'FaceColor','flat')

end


function [x_,y_,z_] = evalCyl2simCart(rho,phi,z)
x_ = -z;
y_ = -rho .* cos(phi);
z_ =  rho .* sin(phi);
end