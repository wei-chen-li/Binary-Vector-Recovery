function plot_sensor_data(filename)

addpath("functions\")

if nargin < 1
    filename = "data\delta-B_circ1mmx4.txt";
end

[B, freqs, ~] = LoadSaveData(filename);

figure
tiledlayout(2, length(freqs), 'TileSpacing','tight', 'Padding','none')

for i = 1:2
for k = 1:length(freqs)
    nexttile
    if i==1, part = @real; else, part = @imag; end
    M = size(B,1);
    C = rot90(reshape(part(B(:,:,k)), sqrt(M),sqrt(M)));
    C_max = max(abs(C), [], 'all');
    imagesc(C, [-C_max C_max])
    axis square off
end
end