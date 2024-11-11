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

sensors.placement_radius = 19e-3;
sensors.placement_angle = linspace(0, 2*pi, 21); sensors.placement_angle = sensors.placement_angle(1:end-1);

comsol.truncate_length = 200e-3;

Model.coil    = coil;
Model.pipe    = pipe;
Model.sensors = sensors;
Model.comsol  = comsol;

end
