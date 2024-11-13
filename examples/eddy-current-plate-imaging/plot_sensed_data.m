clear
clc

delta_B = readmatrix("data\delta-B_circ1mmx4.txt",  'CommentStyle','%');
num_freqs = size(delta_B, 1);

figure
tiledlayout(2,num_freqs, 'TileSpacing','tight', 'Padding','none')

for i = 1:2
for j = 1:num_freqs
    nexttile
    if i==1, part = @real; else, part = @imag; end
    imagesc(rot90(reshape(part(delta_B(j,:)), 8,8)))
    axis square off
end
end