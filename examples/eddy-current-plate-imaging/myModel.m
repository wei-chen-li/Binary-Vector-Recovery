function Model = myModel

comsol_coil_fname = 'myCoil.mph';
cache_coil_fname = 'myCoil_cache.mat';

if isfile(cache_coil_fname) && dir(cache_coil_fname).datenum > dir(comsol_coil_fname).datenum
    load(cache_coil_fname, 'r', 'J', 'dv')
else
    [r, J, dv] = MeshCoil(comsol_coil_fname);
    save(cache_coil_fname, 'r', 'J', 'dv')
end
coil.type = 'arbitrary';
coil.B  = @(x,y,z) B_func (r,J,dv, x,y,z);
coil.Bz = @(x,y,z) Bz_func(r,J,dv, x,y,z);

layers = {};
layers{end+1} = struct('sigma',3.774e7, 'mur',1, 'thickness',3e-3); % Aluminum
layers{end+1} = struct('sigma',0      , 'mur',1, 'thickness',Inf ); % air

x = linspace(-14e-3, 14e-3, 8);
y = linspace(-14e-3, 14e-3, 8);
z = 2e-3;
[x,y,z] = ndgrid(x, y, z);
sensors.positions = [reshape(x,1,[]); reshape(y,1,[]); reshape(z,1,[])];
sensors.axes = zeros(size(sensors.positions));
sensors.axes(2,:) = 1;

x = linspace(-15e-3, 15e-3, 61);
y = linspace(-15e-3, 15e-3, 61);
z = linspace(-3e-3, 0, 7);
mesh = BuildHexMesh(x, y, z, 1);

d = 0;
sigma = zeros(1, size(mesh.elems_center,2));
mur   =  ones(1, size(mesh.elems_center,2));
for k = 1:length(layers)
    d2 = d - layers{k}.thickness;
    z = mesh.elems_center(3,:);
    sigma(d >= z & z > d2) = layers{k}.sigma;
    mur  (d >= z & z > d2) = layers{k}.mur;
    d = d2;
end

comsol.truncate_width = 150e-3; % (m)
comsol.truncate_maxz  =  60e-3; % (m)
comsol.coil_filename = comsol_coil_fname;

Model.layers  = layers;
Model.coil    = coil;
Model.sensors = sensors;
Model.comsol  = comsol;
Model.voxels  = struct('center',mesh.elems_center, 'sigma',sigma, 'mur',mur, 'nodes',mesh.nodes, 'elems',mesh.elems);
Model.Sensitivity = @(f)ComputeSensitivity(layers, mesh.gauss_points, mesh.gauss_weights, coil, sensors.positions, sensors.axes, f);
clear global coilAbovePlate

end


function [r, J, dv] = MeshCoil(comsol_coil_filename)

model = mphopen(comsol_coil_filename);
comp = model.component('comp1');

try    mesh = comp.mesh.create('mesh1');
catch; mesh = comp.mesh('mesh1');
end
mesh.autoMeshSize(4);

study = model.study.create('std1');
study.create('ccc', 'CoilCurrentCalculation');
study.run;

pd = mpheval(model, {'x' 'y' 'z' 'mf.Jex' 'mf.Jey' 'mf.Jez' 'dvol' 'meshtype'}, 'pattern','gauss');
r = [pd.d1; pd.d2; pd.d3];
J = [pd.d4; pd.d5; pd.d6];
meshtype = pd.d8;
assert(all(meshtype == 6)); % all tetrahedrons
dv = pd.d7 * (1/6); % gauss weight for tetrahedron is 1/6

end

function B = B_func(r, J, dv, x,y,z)

assert(isscalar(x))
mu0_over_4pi = 1e-7;

r_eval = [x y z].';
delta_r = r_eval - r;
dB = mu0_over_4pi * cross(J, delta_r ./ vecnorm(delta_r).^3) .* dv;
B = sum(dB, 2);

end

function Bz = Bz_func(r, J, dv, x,y,z)

B = B_func(r, J, dv, x,y,z);
Bz = B(3);

end

function mesh = BuildHexMesh(x_verts, y_verts, z_verts, gauss_order)

assert(isvector(x_verts))
assert(isvector(y_verts))
assert(isvector(z_verts))

n_x = length(x_verts) - 1;
n_y = length(y_verts) - 1;
n_z = length(z_verts) - 1;

[x,y,z] = ndgrid(x_verts, y_verts, z_verts);
nodes = [reshape(x,1,[]); reshape(y,1,[]); reshape(z,1,[])];

elems = zeros(8, n_x*n_y*n_z, 'uint32');
p = 1;
for k = 1:n_z
for j = 1:n_y
for i = 1:n_x
    elems(1,p) = (k-1) * (n_x+1) * (n_y+1) + (j-1) * (n_x+1) + i;
    elems(2,p) = elems(1,p) + 1;
    elems(3,p) = elems(2,p) + (n_x+1);
    elems(4,p) = elems(3,p) - 1;
    elems(5:8,p) = elems(1:4,p) + (n_x+1) * (n_y+1);
    p = p + 1;
end
end
end

num_elems = size(elems,2);
elems_center  = zeros(3, num_elems);
gauss_points  = zeros(3, 0, num_elems);
gauss_weights = zeros(1, 0, num_elems);
for j = 1:size(elems,2)
    elem = elems(:,j);
    elem(elem==0) = [];
    elem_nodes = nodes(:,elem);
    elems_center(:,j) = mean(elem_nodes,2);

    [xis,weights] = GaussQuadratureHex(gauss_order);
    for l = 1:length(weights)
        xi = xis(:,l);
        N     = ShapeFunctionHex(xi);
        dNdxi = ShapeFunctionHex(xi, 'diff');
        x     = elem_nodes * N;
        dxdxi = elem_nodes * dNdxi;
        gauss_points(:,l,j) = x;
        gauss_weights(:,l,j) = det(dxdxi) * weights(l);
    end
end

z = unique(gauss_points(3,:,:));
for k = 2:length(z)
    if z(k) - z(k-1) < 1e-10
        mask = gauss_points(3,:,:) == z(k);
        gauss_points(3,mask) = z(k-1);
        z(k) = z(k-1);
    end
end

mesh.nodes = nodes;
mesh.elems = elems;
mesh.elems_center  = elems_center;
mesh.gauss_points  = gauss_points;
mesh.gauss_weights = gauss_weights;

end

function [xis, wts] = GaussQuadratureHex(gauss_order)

switch gauss_order
    case 0
        xis = [0
               0
               0];
        wts =  8 ;
    case 1
        a = 0.5773502692;
        xis = [-a  a -a  a -a  a -a  a
               -a -a  a  a -a -a  a  a
               -a -a -a -a  a  a  a  a];
        wts = [ 1  1  1  1  1  1  1  1];
    case 2
        a = 0.7745966692;
        w1 = 0.171467763;
        w2 = 0.274348421;
        w3 = 0.438957474;
        w4 = 0.702331959;
        xis = [-a  0  a -a  0  a -a  0  a -a  0  a -a  0  a -a  0  a -a  0  a -a  0  a -a  0  a
               -a -a -a  0  0  0  a  a  a -a -a -a  0  0  0  a  a  a -a -a -a  0  0  0  a  a  a
               -a -a -a -a -a -a -a -a -a  0  0  0  0  0  0  0  0  0  a  a  a  a  a  a  a  a  a];
        wts = [w1 w2 w1 w2 w3 w2 w1 w2 w1 w2 w3 w2 w3 w4 w3 w2 w3 w2 w1 w2 w1 w2 w3 w2 w1 w2 w1];
end

end

function val = ShapeFunctionHex(xi, varargin)

xi1 = xi(1); xi2 = xi(2); xi3 = xi(3);
is_diff = ~isempty(varargin);
if ~is_diff
    N = zeros(8, 1);
    N(1) = (1 - xi1) * (1 - xi2) * (1 - xi3) / 8;
    N(2) = (1 + xi1) * (1 - xi2) * (1 - xi3) / 8;
    N(3) = (1 + xi1) * (1 + xi2) * (1 - xi3) / 8;
    N(4) = (1 - xi1) * (1 + xi2) * (1 - xi3) / 8;
    N(5) = (1 - xi1) * (1 - xi2) * (1 + xi3) / 8;
    N(6) = (1 + xi1) * (1 - xi2) * (1 + xi3) / 8;
    N(7) = (1 + xi1) * (1 + xi2) * (1 + xi3) / 8;
    N(8) = (1 - xi1) * (1 + xi2) * (1 + xi3) / 8;
    val = N;
else
    dNdxi = zeros(8, 3);
    dNdxi(1,1) = -(1 - xi2) * (1 - xi3) / 8;
    dNdxi(1,2) = -(1 - xi1) * (1 - xi3) / 8;
    dNdxi(1,3) = -(1 - xi1) * (1 - xi2) / 8;
    dNdxi(2,1) =  (1 - xi2) * (1 - xi3) / 8;
    dNdxi(2,2) = -(1 + xi1) * (1 - xi3) / 8;
    dNdxi(2,3) = -(1 + xi1) * (1 - xi2) / 8;
    dNdxi(3,1) =  (1 + xi2) * (1 - xi3) / 8;
    dNdxi(3,2) =  (1 + xi1) * (1 - xi3) / 8;
    dNdxi(3,3) = -(1 + xi1) * (1 + xi2) / 8;
    dNdxi(4,1) = -(1 + xi2) * (1 - xi3) / 8;
    dNdxi(4,2) =  (1 - xi1) * (1 - xi3) / 8;
    dNdxi(4,3) = -(1 - xi1) * (1 + xi2) / 8;
    dNdxi(5,1) = -(1 - xi2) * (1 + xi3) / 8;
    dNdxi(5,2) = -(1 - xi1) * (1 + xi3) / 8;
    dNdxi(5,3) =  (1 - xi1) * (1 - xi2) / 8;
    dNdxi(6,1) =  (1 - xi2) * (1 + xi3) / 8;
    dNdxi(6,2) = -(1 + xi1) * (1 + xi3) / 8;
    dNdxi(6,3) =  (1 + xi1) * (1 - xi2) / 8;
    dNdxi(7,1) =  (1 + xi2) * (1 + xi3) / 8;
    dNdxi(7,2) =  (1 + xi1) * (1 + xi3) / 8;
    dNdxi(7,3) =  (1 + xi1) * (1 + xi2) / 8;
    dNdxi(8,1) = -(1 + xi2) * (1 + xi3) / 8;
    dNdxi(8,2) =  (1 - xi1) * (1 + xi3) / 8;
    dNdxi(8,3) =  (1 - xi1) * (1 + xi2) / 8;
    val = dNdxi;
end

end

function [S_sigma, S_mu] = ComputeSensitivity(layers, gauss_points, gauss_weights, coil, sensor_positions, sensor_axes, freq)

omega = 2*pi * freq;
x = gauss_points(1,:,:);
y = gauss_points(2,:,:);
z = gauss_points(3,:,:);
num_sensors = size(sensor_positions, 2);
num_elems = size(gauss_points, 3);

global coilAbovePlate
if isempty(coilAbovePlate)
    coilAbovePlate = CoilAbovePlate3D(coil, layers, 'TruncateWidth',0.15);
end

if nargout >= 1
    [Ex,Ey] = coilAbovePlate.EvaluateField({'Ex' 'Ey'}, x,y,z, omega);
    E1 = [Ex;Ey];
    S_sigma = zeros(num_sensors, num_elems);
end
if nargout >= 2
    [Hx,Hy,Hz] = coilAbovePlate.EvaluateField({'Hx' 'Hy' 'Hz'}, x,y,z, omega);
    H1 = [Hx;Hy;Hz];
    S_mu = zeros(num_sensors, num_elems);
end

for i = 1:num_sensors
    sensor.type = 'MagneticSensor';
    sensor.position = [0 0 sensor_positions(3,i)];
    sensor.axis = sensor_axes(:,i);
    sensorAbovePlate = CoilAbovePlate3D(sensor, layers);

    if nargout >= 1
        [Ex,Ey] = sensorAbovePlate.EvaluateField({'Ex' 'Ey'}, ...
            x-sensor_positions(1,i), ...
            y-sensor_positions(2,i), z, omega);
        E2 = [Ex;Ey];
        S_sigma(i,:) = sum(dot(E1,E2) .* gauss_weights, 2);
    end

    if nargout >= 2
        [Hx,Hy,Hz] = sensorAbovePlate.EvaluateField({'Hx' 'Hy' 'Hz'}, ...
            x-sensor_positions(1,i), ...
            y-sensor_positions(2,i), z, omega);
        H2 = [Hx;Hy;Hz];    
        S_mu(i,:) = sum(dot(H1,H2) .* gauss_weights, 2) * -1j * omega;
    end
end

end
