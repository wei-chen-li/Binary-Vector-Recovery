function varargout = LoadSaveData(filename, B, freqs, defects)

if ~(nargin == 1 || nargin == 4)
    error("Incorrect number of argments")
end

if nargin == 4
    saveData(filename, B, freqs, defects);
else
    varargout = cell(1,3);
    [varargout{:}] = loadData(filename);
end

end


function saveData(filename, B, freqs, defects)

if size(B,2) ~= 1
    error("B sould be column vectors")
end

fid = fopen(filename, 'w');

for i = 1:length(defects)
    fprintf(fid, '%% defect %s\n', jsonencode(defects{i}));
end
for i = 1:length(freqs)
    fprintf(fid, '%% frequency %g Hz\n', freqs(i));
    fprintf(fid, '%s\n', num2str(B(:,:,i).'));
end

fclose(fid);

end


function [B, freqs, defects] = loadData(filename)

B = zeros(0,0,0);
freqs = [];
defects = {};

lines = readlines(filename);

for k = 1:length(lines)
    str = strtrim(lines(k));
    str = convertStringsToChars(str);
    if isempty(str)
        continue
    end

    if ~startsWith(str, '%')
        B(:,:,end+1) = str2num(str).';
    else
        str = strtrim(str(2:end));
        if startsWith(str, 'frequency')
            tokens = split(str);
            freqs(end+1) = str2num(tokens{2});
            assert(strcmp(tokens{3}, "Hz"))
        elseif startsWith(str, 'defect')
            i = strfind(str, ' ');
            defects{end+1} = jsondecode(str(i:end));
        else
            error('Unknown "%s"', str)
        end
    end
end

end