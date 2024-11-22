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

phi = linspace(0, 2*pi, 21); phi = phi(1:end-1);
sensors.placement_phi = phi;
sensors.placement_rho = 19e-3 * ones(size(phi));

comsol.truncate_length = 200e-3;
comsol.sensor_size = 1e-3;

Model.coil    = coil;
Model.pipe    = pipe;
Model.sensors = sensors;
Model.comsol  = comsol;

end
