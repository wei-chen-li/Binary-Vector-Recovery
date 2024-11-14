function plot_sensitivity(filename)

if nargin < 1
    filename = 'data\Phi.mat';
end
load(filename, 'Phi','freqs','voxels_center')

num_cols = length(freqs);
num_rows = length(unique(voxels_center(3,:)));

figure
t = tiledlayout(num_rows,num_cols, 'TileSpacing','tight', 'Padding','none');

for k = 1:num_cols
    VisualizeSinglePhi(Phi, k, voxels_center, t)
end

end


function VisualizeSinglePhi(Phi, k, voxels_center, tiled_chart_layout)

C = vecnorm(Phi(:,:,k), 2,1);

z = unique(voxels_center(3,:));
z = sort(z, 'descend');

for i = 1:length(z)
    mask = voxels_center(3,:) == z(i);
    C_layer = reshape(C(mask), sqrt(sum(mask)),[]).';

    x = [min(voxels_center(1,mask)) max(voxels_center(1,mask))] * 1e3;
    y = [min(voxels_center(2,mask)) max(voxels_center(2,mask))] * 1e3;

    nexttile(tilenum(tiled_chart_layout, i,k))

    imagesc(x, y, C_layer, [0 max(C)])
    axis square xy
end

end
