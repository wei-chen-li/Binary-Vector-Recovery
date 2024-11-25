clear
clc

figure
hold on
axis equal off

Model = myModel;
drawSensitivity(Model)
% drawPipe(Model.pipe, 100e-3)
% drawCoil(Model.coil)
% drawSensors(Model.sensors)
% drawMesh(Model.mesh)

view(-30,30)
camtarget([0 0 0])
% camlight
camva(6)


function drawCoil(coil, color, alpha)

if nargin < 3, color = [1 0 0]; end
if nargin < 4, alpha = 0.5; end

drawCylinder(coil.inner_radius, coil.outer_radius, 0, 2*pi,  coil.spacing/2,  coil.spacing/2+coil.height, color, alpha)
drawCylinder(coil.inner_radius, coil.outer_radius, 0, 2*pi, -coil.spacing/2, -coil.spacing/2-coil.height, color, alpha)

end


function drawPipe(pipe, pipe_length, color, alpha)

if nargin < 3, color = [0.7 0.7 0.7]; end
if nargin < 4, alpha = 1; end

drawCylinder(pipe.inner_diameter/2, pipe.outer_diameter/2, pi/2, 2*pi, -pipe_length/2, pipe_length/2, color, alpha)

end


function drawSensors(sensors, color, alpha)

if nargin < 2, color = [0 0 0]; end
if nargin < 3, alpha = 0.7; end

rho = sensors.placement_rho;
phi = sensors.placement_phi;

nodes = table2array(combinations([-1 1] * 3e-4, [-1 1] * 6e-4, [-1 1] * 2e-4))';
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


function drawCylinder(rho1, rho2, phi1, phi2, z1, z2, color, alpha)

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


function drawMesh(mesh)

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

end


function drawSensitivity(Model)

addpath("functions\")

freq = 1000;
mesh_size = 5;

num_sensors = length(Model.sensors.placement_phi);
shift = 15;

rho = (Model.pipe.inner_diameter + Model.pipe.outer_diameter) / 4;
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