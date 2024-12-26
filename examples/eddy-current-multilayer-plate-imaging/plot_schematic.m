clear
clc

figure

model = myModel;
drawCoil()
drawPlate(model)
drawMesh(model.voxels.nodes, [0 0 0.3])

for k = 1:size(model.sensors.positions,2)
    drawSensor(model.sensors.positions(:,k), model.sensors.axes(:,k))
end

view(-30,20)
axis equal off
camtarget([0 0 0])
camlight
camva(1)


function drawCoil(alpha)

if nargin < 1, alpha = 0.4; end

coilFiles = findFiles('coil\d*\.(step|stp|stl)');
for k = 1:length(coilFiles)
    gm = importGeometry(coilFiles{k});
    gm = scale(gm, 1e-3);
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


function drawPlate(model, alpha)

if nargin < 2, alpha = 0.6; end

w = model.comsol.truncate_width;

faces = [[1 5 6 2]; [1 2 4 3]; [2 6 8 4]];

nodes = table2array(combinations([-w/2 w/2], [-w/2 w/2], [-1e-3  0e-3]))';
patch('Faces',faces, 'Vertices',nodes', 'FaceColor',[0.62 0.62 0.62], 'LineStyle','none', 'FaceAlpha',alpha)

nodes = table2array(combinations([-w/2 w/2], [-w/2 w/2], [-2e-3 -1e-3]))';
patch('Faces',faces, 'Vertices',nodes', 'FaceColor',[0.60 0.50 0.30], 'LineStyle','none', 'FaceAlpha',alpha)

end


function drawSensor(position, axis) 

nodes = table2array(combinations([-1 1] * 0.75e-3, [-1 1] * 0.75e-3, [-1 1] * 0.225e-3))';
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