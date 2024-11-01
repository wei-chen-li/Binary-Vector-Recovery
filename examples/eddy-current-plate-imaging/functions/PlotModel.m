function PlotModel(model)

figure

coilFiles = findFiles('coil\d*\.(step|stp|stl)');
for k = 1:length(coilFiles)
    gm = importGeometry(coilFiles{k});
    h = pdegplot(gm, 'Lighting','off');
    h(1).FaceColor = [1 0 0];
    h(1).FaceAlpha = 0.5;
    h(2).Color = 'none';
    hold on; axis off;
end

strucFiles = findFiles('(plate|pipe)\d*\.(step|stp|stl)');
for k = 1:length(strucFiles)
    gm = importGeometry(strucFiles{k});
    h = pdegplot(gm, 'Lighting','off');
    h(1).FaceColor = [0.8 0.8 0.8];
    h(1).FaceAlpha = 0.5;
    h(2).Color = 'none';
    hold on; axis off;
end

originQuiver = findobj(gca,'Type','Quiver');
set(originQuiver,'Visible','off');
originText = findobj(gca,'String','x','-or','String','y','-or','String','z');
set(originText,'Visible','off');

for k = 1:size(model.sensors.positions,2)
    drawSensors(model.sensors.positions(:,k), model.sensors.axes(:,k))
end

drawMesh(model.voxels.nodes)

end


function drawSensors(position, axis)

a = 0.8e-3;
N = 100;
theta = linspace(0, 2*pi, N);
nodes = [a*cos(theta) a*cos(theta)
         a*sin(theta) a*sin(theta)
         a*ones(1,N)  -a*ones(1,N)];
p3 = reshape(axis,3,1) / norm(axis);
p2 = [-p3(2); p3(1); 0]; p2 = p2 / norm(p2);
p1 = cross(p2, p3);
T = [p1 p2 p3];
nodes = T * nodes + position;

top = patch(nodes(1,  1:  N), nodes(2,  1:  N), nodes(3,  1:  N), '');
bot = patch(nodes(1,N+1:2*N), nodes(2,N+1:2*N), nodes(3,N+1:2*N), '');
out = surf(reshape(nodes(1,:),N,2), reshape(nodes(2,:),N,2), reshape(nodes(3,:),N,2));

set([bot,top,out], 'FaceColor',[0 0 0], 'FaceAlpha',0.5, 'LineStyle','none');

end


function drawMesh(nodes)

x = unique(nodes(1,:));
y = unique(nodes(2,:));
z = unique(nodes(3,:));

I = [1 1];

[ys,zs] = ndgrid(y,z);
for k = 1:numel(ys)
    plot3([min(x) max(x)], ys(k)*I, zs(k)*I, 'k')
end
[zs,xs] = ndgrid(z,x);
for k = 1:numel(zs)
    plot3(xs(k)*I, [min(y) max(y)], zs(k)*I, 'k')
end
[xs,ys] = ndgrid(x,y);
for k = 1:numel(xs)
    plot3(xs(k)*I, ys(k)*I, [min(z) max(z)], 'k')
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