% ===============================
% User-defined parameters
% ===============================
n = 256;  % Total number of colors in colormap
threshold = 1.5;

% Binary colors at 0 +/- threshold
negative_bin_color = [0.5, 0.0, 0.5];  % Purple at -threshold  % Dark blue at -threshold
positive_bin_color = [1.0, 0.5, 0.0];  % Orange at +threshold  % Dark red at +threshold

% Gradient colors for beyond-threshold values
neg_grad_start = [0.2, 0.2, 0.6];  % Dark blue (at -threshold)
neg_grad_end   = [0.1, 0.1, 0.3];  % Very dark blue (beyond -threshold)

pos_grad_start = [0.6, 0.2, 0.2];  % Dark red (at +threshold)
pos_grad_end   = [0.9, 0.6, 0.6];  % Light red (beyond +threshold)

% ===============================
% Generate the colormap
% ===============================
abs_thresh = abs(threshold);

% Define segments: [neg_bin | neg_grad | pos_grad | pos_bin]
segment_ratio = [0.1, 0.4, 0.4, 0.1];  % Adjustable
segment_counts = round(segment_ratio * n);
segment_counts(end) = n - sum(segment_counts(1:end-1));  % Adjust to sum to n

% Segment 1: Negative bin at -threshold
c_neg_bin = repmat(negative_bin_color, segment_counts(1), 1);

% Segment 2: Negative gradient (beyond -threshold)
c_neg_grad = [linspace(neg_grad_start(1), neg_grad_end(1), segment_counts(2))', ...
              linspace(neg_grad_start(2), neg_grad_end(2), segment_counts(2))', ...
              linspace(neg_grad_start(3), neg_grad_end(3), segment_counts(2))'];

% Segment 3: Positive gradient (beyond +threshold)
c_pos_grad = [linspace(pos_grad_start(1), pos_grad_end(1), segment_counts(3))', ...
              linspace(pos_grad_start(2), pos_grad_end(2), segment_counts(3))', ...
              linspace(pos_grad_start(3), pos_grad_end(3), segment_counts(3))'];

% Segment 4: Positive bin at +threshold
c_pos_bin = repmat(positive_bin_color, segment_counts(4), 1);

% Combine all
cmap = [ c_neg_grad;c_neg_bin; c_pos_bin;c_pos_grad];

% ===============================
% Save colormap to file
% ===============================
filename = '/home/mac/mlin2/cluster_analysis/brainnet config/custom_colormap.txt';
[filepath, ~, ~] = fileparts(filename);
if ~isfolder(filepath)
    error('Directory does not exist: %s', filepath);
end

fileID = fopen(filename, 'w');
if fileID == -1
    error('Failed to open file: %s', filename);
end

fprintf(fileID, '[\n');
for i = 1:n
    fprintf(fileID, '  %.6f %.6f %.6f;\n', cmap(i, :));
end
fprintf(fileID, ']\n');
fclose(fileID);

% ===============================
% Optional: Visualize colormap
% ===============================
figure;
colormap(cmap);
colorbar;
title('Custom Colormap');

% Confirmation
fprintf('Colormap saved to %s\n', filename);

%% mri
% ===============================
% User-defined parameters
% ===============================
n = 256;  % Total number of colors in colormap

% Anchor colors (approximate from parula)
start_color = [0.5081, 0.2663, 0.6292];  % Purple-ish (start
mid_color   = [0.1, 0.4, 0.6];  % Blue-ish (middle)
end_color   = [0.2886, 0.7604, 0.4862];  % Green-ish (end
%end_color   = [0.6, 0.2, 0.8];    % Vivid purple/magenta

% ===============================
% Generate the continuous colormap
% ===============================
n1 = round(n / 2);
n2 = n - n1;

% Gradient 1: Green to Blue
c1 = [linspace(start_color(1), mid_color(1), n1)', ...
      linspace(start_color(2), mid_color(2), n1)', ...
      linspace(start_color(3), mid_color(3), n1)'];

% Gradient 2: Blue to Purple
c2 = [linspace(mid_color(1), end_color(1), n2)', ...
      linspace(mid_color(2), end_color(2), n2)', ...
      linspace(mid_color(3), end_color(3), n2)'];

% Combine into full colormap
cmap = [c1; c2];

% ===============================
% Save colormap to file
% ===============================
filename = '/home/mac/mlin2/cluster_analysis/brainnet config/mri_bgp.txt';
[filepath, ~, ~] = fileparts(filename);
if ~isfolder(filepath)
    error('Directory does not exist: %s', filepath);
end

fileID = fopen(filename, 'w');
if fileID == -1
    error('Failed to open file: %s', filename);
end

fprintf(fileID, '[\n');
for i = 1:n
    fprintf(fileID, '  %.6f %.6f %.6f;\n', cmap(i, :));
end
fprintf(fileID, ']\n');
fclose(fileID);

% ===============================
% Optional: Visualize colormap
% ===============================
figure;
colormap(cmap);
colorbar;
title('MRI B-G-P Colormap');

% Confirmation
fprintf('Colormap saved to %s\n', filename);

