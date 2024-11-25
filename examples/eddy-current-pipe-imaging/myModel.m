function Model = myModel

coil = {};
coil.inner_radius = 18e-3;
coil.outer_radius = 20e-3;
coil.height = 4e-3;
coil.spacing = 16e-3;
coil.turns = 20;

pipe = struct('sigma',3.774e7, 'mur',1);
pipe.inner_diameter = 45e-3;
pipe.outer_diameter = 50e-3;

num_sensors = 20;
phi = linspace(0, 2*pi, num_sensors+1); phi = phi(1:end-1);
sensors.placement_phi = phi;
sensors.placement_rho = 19e-3 * ones(size(phi));

comsol.truncate_length = 200e-3;
comsol.sensor_size = 1e-3;

rho = [22.5 25] * 1e-3;
phi = linspace(0, 2*pi, 101);
z = (-9:1.5:-3) * 1e-3;
[voxel_corners_1, voxel_volumes_1, gauss_nodes_1] = buildVoxels(rho, phi, z);
[voxel_corners_2, voxel_volumes_2, gauss_nodes_2] = buildVoxels(rho, phi, flip(-z));
voxel_corners = [voxel_corners_1 voxel_corners_2];
voxel_volumes = [voxel_volumes_1 voxel_volumes_2];
gauss_nodes   = [gauss_nodes_1   gauss_nodes_2  ];
mesh = struct('voxel_corners',voxel_corners, 'voxel_volumes',voxel_volumes, 'gauss_nodes',gauss_nodes);

Model.coil    = coil;
Model.pipe    = pipe;
Model.sensors = sensors;
Model.comsol  = comsol;
Model.mesh    = mesh;

end


function [voxel_corners, voxel_volumes, gauss_nodes] = buildVoxels(rho_, phi_, z_)

[rho, phi, z] = ndgrid(rho_, phi_, z_);

num_voxels = (length(rho_)-1) * (length(phi_)-1) * (length(z_)-1);

voxel_corners = zeros(3, num_voxels, 8);
voxel_volumes = zeros(1, num_voxels);
gauss_nodes   = zeros(3, num_voxels, 8);

a = 1 / sqrt(3);
xis = [-a  a -a  a -a  a -a  a
       -a -a  a  a -a -a  a  a
       -a -a -a -a  a  a  a  a];

n = 0;
for k = 1:length(z_)-1
for j = 1:length(phi_)-1
for i = 1:length(rho_)-1
    n = n + 1;

    voxel_corners(:,n,1) = [rho(i,  j  ,k  )  phi(i,  j,  k  )  z(i,  j,  k  )];
    voxel_corners(:,n,2) = [rho(i+1,j,  k  )  phi(i+1,j,  k  )  z(i+1,j,  k  )];
    voxel_corners(:,n,3) = [rho(i+1,j+1,k  )  phi(i+1,j+1,k  )  z(i+1,j+1,k  )];
    voxel_corners(:,n,4) = [rho(i,  j+1,k  )  phi(i,  j+1,k  )  z(i,  j+1,k  )];
    voxel_corners(:,n,5) = [rho(i,  j  ,k+1)  phi(i,  j,  k+1)  z(i,  j,  k+1)];
    voxel_corners(:,n,6) = [rho(i+1,j,  k+1)  phi(i+1,j,  k+1)  z(i+1,j,  k+1)];
    voxel_corners(:,n,7) = [rho(i+1,j+1,k+1)  phi(i+1,j+1,k+1)  z(i+1,j+1,k+1)];
    voxel_corners(:,n,8) = [rho(i,  j+1,k+1)  phi(i,  j+1,k+1)  z(i,  j+1,k+1)];

    voxel_volumes(n) = pi * (rho_(i+1)^2 - rho_(i)^2) * (phi_(j+1) - phi_(j)) / (2*pi) * (z_(k+1) - z_(k));
    voxel_volumes(n) = abs(voxel_volumes(n));

    for q = 1:size(xis,2)
        N = ShapeFunctionHex(xis(:,q));
        N = reshape(N, 1,1,8);
        gauss_nodes(:,n,q) = sum(voxel_corners(:,n,:) .* N, 3);
    end
end
end
end

end


function N = ShapeFunctionHex(xi)

xi1 = xi(1); xi2 = xi(2); xi3 = xi(3);

N = zeros(8,1);
N(1) = (1 - xi1) * (1 - xi2) * (1 - xi3) / 8;
N(2) = (1 + xi1) * (1 - xi2) * (1 - xi3) / 8;
N(3) = (1 + xi1) * (1 + xi2) * (1 - xi3) / 8;
N(4) = (1 - xi1) * (1 + xi2) * (1 - xi3) / 8;
N(5) = (1 - xi1) * (1 - xi2) * (1 + xi3) / 8;
N(6) = (1 + xi1) * (1 - xi2) * (1 + xi3) / 8;
N(7) = (1 + xi1) * (1 + xi2) * (1 + xi3) / 8;
N(8) = (1 - xi1) * (1 + xi2) * (1 + xi3) / 8;

end