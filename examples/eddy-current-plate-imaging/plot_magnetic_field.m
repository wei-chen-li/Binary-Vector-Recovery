function plot_magnetic_field(filename)

if nargin < 1
    filename = "data\delta-B_circ1mmx4.txt";
end

delta_B = readmatrix(filename,  'CommentStyle','%');
num_freqs = size(delta_B, 1);

figure
tiledlayout(2,num_freqs, 'TileSpacing','tight', 'Padding','none')

for i = 1:2
for j = 1:num_freqs
    nexttile
    if i==1, part = @real; else, part = @imag; end
    M = size(delta_B,2);
    C = rot90(reshape(part(delta_B(j,:)), sqrt(M),sqrt(M)));
    C_max = max(abs(C), [], 'all');
    imagesc(C, [-C_max C_max])
    axis square off
end
end