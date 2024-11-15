function plot_defects(defects, voxels_center)

addpath("functions\")

if nargin < 2
    load('data\Phi.mat', 'voxels_center')
end
if nargin < 1
    [~, ~, defects] = LoadSaveData("data\delta-B_circ2mmx4.txt");
end

x = unique(voxels_center(1,:));
y = unique(voxels_center(2,:));
z = unique(voxels_center(3,:));
x = [x(1)-(x(2)-x(1))/2 x(end)+(x(2)-x(1))/2];
y = [y(1)-(y(2)-y(1))/2 y(end)+(y(2)-y(1))/2];
z = [z(1)-(z(2)-z(1))/2 z(end)+(z(2)-z(1))/2];

fig_size = 600;

figure('Position',[100 150 fig_size fig_size])
tiledlayout(1,1, 'TileSpacing','none', 'Padding','none')
nexttile
xticks([]); yticks([]); xlim(x); ylim(y); axis square; box on; hold on;

figure('Position',[150+fig_size 150 fig_size fig_size/(y(2)-y(1))*(z(2)-z(1))])
tiledlayout(1,1, 'TileSpacing','none', 'Padding','none')
nexttile
xticks([]); yticks([]); xlim(y); ylim(z); box on; hold on;

figure('Position',[100 50 fig_size fig_size/(x(2)-x(1))*(z(2)-z(1))])
tiledlayout(1,1, 'TileSpacing','none', 'Padding','none')
nexttile
xticks([]); yticks([]); xlim(x); ylim(z); box on; hold on;

for p = 1:length(defects)
    defect = defects{p};
    switch lower(defect.type)
        case 'rectangular'
            x = [defect.c_x-defect.dx/2 defect.c_x+defect.dx/2 defect.c_x+defect.dx/2 defect.c_x-defect.dx/2];
            y = [defect.c_y-defect.dy/2 defect.c_y-defect.dy/2 defect.c_y+defect.dy/2 defect.c_y+defect.dy/2];
            x = [x x(1)]; y = [y y(1)];
        case 'circular'
            theta = linspace(0, 2*pi);
            x = defect.c_x + defect.d/2 * cos(theta);
            y = defect.c_y + defect.d/2 * sin(theta);
        case 'polygonal'
            x = defect.xs;
            y = defect.ys;
        otherwise
            error('Unknown defect shape "%s"', defect.type)
    end

    figure(1)
    if max(defect.z1, defect.z2) < -1e-5, linestyle = '--'; else, linestyle = '-'; end
    plot(x, y, 'k', 'LineStyle',linestyle)

    figure(2)
    a = [defect.z1 defect.z2 defect.z2 defect.z1];
    b = [min(y) min(y) max(y) max(y)];
    a = [a a(1)]; b = [b b(1)];
    plot(b, a, 'k', 'LineStyle','--')

    figure(3)
    a = [min(x) min(x) max(x) max(x)];
    b = [defect.z1 defect.z2 defect.z2 defect.z1];
    a = [a a(1)]; b = [b b(1)];
    plot(a, b, 'k', 'LineStyle','--')
end

end