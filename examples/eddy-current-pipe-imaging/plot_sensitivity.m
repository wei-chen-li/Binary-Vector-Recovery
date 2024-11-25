function plot_sensitivity(filename, sensor_num)

if nargin < 2
    sensor_num = [];
end

if nargin < 1
    filename = 'data\Phi.mat';
end
load(filename, 'Phi','freqs','voxels_corner')

rho = unique(voxels_corner(1,:));

num_cols = length(freqs);
num_rows = length(rho) - 1;

figure
t = tiledlayout(num_rows,num_cols, 'TileSpacing','tight', 'Padding','none');

for k = 1:num_cols
    VisualizeSinglePhi(Phi, k, sensor_num, voxels_corner, t)
end

end


function VisualizeSinglePhi(Phi, k, sensor_num, voxels_corner, tiled_chart_layout)

if isempty(sensor_num)
    C = vecnorm(Phi(:,:,k), 2,1);
else
    C = vecnorm(Phi(sensor_num,:,k), 2,1);
end

rho = unique(voxels_corner(1,:));
phi = unique(voxels_corner(2,:));
z   = unique(voxels_corner(3,:));

z1 = z(1:length(z)/2);
z2 = z(length(z)/2+1:end);

[phi1,z1] = ndgrid(phi,z1);
[phi2,z2] = ndgrid(phi,z2);

for i = 1:length(rho)-1
    mask = voxels_corner(1,:,1) == rho(i);
    C_layer = reshape(C(mask), length(phi)-1, []);
    C_layer_1 = C_layer(:,1:length(z)/2-1);
    C_layer_2 = C_layer(:,length(z)/2:end);

    [x1_,y1_,z1_] = evalCyl2simCart(rho(i), phi1, z1);
    [x2_,y2_,z2_] = evalCyl2simCart(rho(i), phi2, z2);

    nexttile(tilenum(tiled_chart_layout, i,k))
    hold on

    surf(x1_, y1_, z1_, C_layer_1, 'LineStyle','none')
    surf(x2_, y2_, z2_, C_layer_2, 'LineStyle','none')
    axis equal off
    view(-30,30)
    camtarget([0 0 0])
end

end


function [x_,y_,z_] = evalCyl2simCart(rho,phi,z)
x_ = -z;
y_ = -rho .* cos(phi);
z_ =  rho .* sin(phi);
end
