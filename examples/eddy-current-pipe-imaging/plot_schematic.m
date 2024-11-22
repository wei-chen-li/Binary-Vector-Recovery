clear
clc

figure
hold on
axis equal off

model = myModel;
drawPipe(model.pipe, 100e-3)
drawCoil(model.coil)
drawSensors(model.sensors)

view(-30,30)
camtarget([0 0 0])
camlight
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
    disp(pos)
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