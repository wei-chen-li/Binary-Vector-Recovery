clear
clc

figure

model = myModel;
drawCoil()
drawPlate(model)
drawEddyCurrent(model)
% drawMesh(model.voxels.nodes)

for k = 1:size(model.sensors.positions,2)
    drawSensor(model.sensors.positions(:,k), model.sensors.axes(:,k))
end

view(-30,20)
camtarget([0 0 0])
camlight
camva(1)


function drawCoil(alpha)

if nargin < 1, alpha = 0.3; end

coilFiles = findFiles('coil\d*\.(step|stp|stl)');
for k = 1:length(coilFiles)
    gm = importGeometry(coilFiles{k});
    h = pdegplot(gm, 'Lighting','off');
    h(1).FaceColor = [1 0 0];
    h(1).FaceAlpha = alpha;
    h(2).Color = 'none';
    hold on; axis off;
end

originQuiver = findobj(gca,'Type','Quiver');
set(originQuiver,'Visible','off');
originText = findobj(gca,'String','x','-or','String','y','-or','String','z');
set(originText,'Visible','off');

end


function drawPlate(model, color, alpha)

if nargin < 2, color = [0.6 0.6 0.6]; end
if nargin < 3, alpha = 0.5; end

w = model.comsol.truncate_width;
nodes = table2array(combinations([-w/2 w/2], [-w/2 w/2], [-model.plate.thickness 0]))';

faces = [[1 5 6 2]; [1 2 4 3]];
patch('Faces',faces, 'Vertices',nodes', 'FaceColor',color, 'LineStyle','none', 'FaceAlpha',1)

faces = [[2 6 8 4]];
patch('Faces',faces, 'Vertices',nodes', 'FaceColor',color, 'LineStyle','none', 'FaceAlpha',alpha)

end


function drawEddyCurrent(model, linewidth)

addpath("functions\")

if nargin < 2, linewidth = 1.5; end

w = model.comsol.truncate_width;

layers = {};
layers{end+1} = model.plate;
layers{end+1} = struct('sigma',0, 'mur',1, 'thickness',Inf); % air
coilAbovePlate = CoilAbovePlate(model.coil, layers, 'TruncateWidth',w);

x = linspace(-w/2, w/2, 61);
[x,y] = ndgrid(x,x);
z = -model.plate.thickness / 2 * ones(size(x));

[Jx, Jy] = coilAbovePlate.EvaluateField({'Jx' 'Jy'}, x,y,z, 2*pi*1000);

quiver(x, y, imag(Jx), imag(Jy), 'Color',[1 0 0], 'LineWidth',linewidth)

end


function drawSensor(position, axis) 

nodes = table2array(combinations([-1 1] * 3e-4, [-1 1] * 6e-4, [-1 1] * 2e-4))';
faces = [[1 3 7 5]; [2 6 8 4]; [1 5 6 2]; [3 4 8 7]; [1 2 4 3]; [5 7 8 6]];

p2 = reshape(axis,3,1) / norm(axis);
p1 = [p2(2); -p2(1); 0]; p1 = p1 / norm(p1);
p3 = cross(p1, p2);
T = [p1 p2 p3];
nodes = T * nodes + position;

patch('Faces',faces, 'Vertices',nodes', ...
      'FaceColor', [0 0 0], ...
      'FaceAlpha', 0.7, ...
      'LineStyle', 'none')

end


function drawMesh(nodes, color, linewidth)

if nargin < 2, color = [0 0 0]; end
if nargin < 3, linewidth = 0.01; end

x = unique(nodes(1,:));
y = unique(nodes(2,:));
z = unique(nodes(3,:));

I = [1 1];

[ys,zs] = ndgrid(y,z);
for k = 1:numel(ys)
    if ~(ys(k) == min(y) || zs(k) == max(z)), continue; end
    plot3([min(x) max(x)], ys(k)*I, zs(k)*I, 'Color',color, 'LineWidth',linewidth)
end
[zs,xs] = ndgrid(z,x);
for k = 1:numel(zs)
    if ~(xs(k) == min(x) || zs(k) == max(z)), continue; end
    plot3(xs(k)*I, [min(y) max(y)], zs(k)*I, 'Color',color, 'LineWidth',linewidth)
end
[xs,ys] = ndgrid(x,y);
for k = 1:numel(xs)
    if ~(xs(k) == min(x) || ys(k) == min(y)), continue; end
    plot3(xs(k)*I, ys(k)*I, [min(z) max(z)], 'Color',color, 'LineWidth',linewidth)
end

end


function matchingFiles = findFiles(pattern)

searchDir = 'assets/';
matchingFiles = {};

fileList = dir(fullfile(searchDir, '*'));

for i = 1:length(fileList)
    if fileList(i).isdir
        continue;
    end
    
    [~, name, ext] = fileparts(fileList(i).name);
    fileName = strcat(name, ext);
    
    if ~isempty(regexp(fileName, pattern, 'once'))
        matchingFiles{end+1} = fullfile(searchDir, fileList(i).name);
    end
end

end