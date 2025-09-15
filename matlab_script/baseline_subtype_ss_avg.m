%% basic setting
% Correct usage of fullfile to construct paths
root = '/home/mac/mlin2/cluster_analysis';
filename = '/type assignment/sense_pt8_threshold.csv';%final_inter_c2mean_2_bl_ss%inter_c2mean_ss_fname_bl
% Configuration for visualization paths
surface = 'ICBM';
surface_path = fullfile('/home/mac/mlin2/Downloads/utilities/BrainNetViewer_20191031/Data/SurfTemplate', 'BrainMesh_ICBM152_smoothed.nv');
cfig = fullfile(root,'brainnet config', 'cet_rbgyrm1_pt5_3.mat');
% Example of fullfile to create a complete file path for the CSV file
fdata_path = fullfile(root, filename);
fdata = readtable(fdata_path, 'ReadVariableNames', true);
disp('read');

%% subtype mean
% Parameters
subtype_label = '3x';  % Only used in file naming
output_dir = fullfile(root, 'roi_inter_c2means_thrept8');
fcol = fdata.fname;
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Directory created: %s\n', output_dir);
else
    fprintf('Directory already exists: %s\n', output_dir);
end

% Get unique subtype clusters (excluding 99/missing)
unique_clusters = unique(fdata.subtype);
unique_clusters = unique_clusters(unique_clusters ~= 99);

% Initialize nii_files cell
nii_files = cell(1, max(unique_clusters));

% ===== Generate .nii Files =====
for cluster_num = unique_clusters'
    % Find indices for this subtype cluster
    cluster_indices = fdata.subtype == cluster_num;

    % Get image paths directly from fcol
    suvr_image_paths = fcol(cluster_indices);

    % Convert to cellstr in case it's not
    suvr_image_paths = cellstr(suvr_image_paths);

    % Prepare output filename
    output_filename = sprintf('mean_%s_s%d.nii', subtype_label, cluster_num);

    % Build SPM imcalc batch
    matlabbatch{1}.spm.util.imcalc.input = suvr_image_paths;
    matlabbatch{1}.spm.util.imcalc.output = output_filename;
    matlabbatch{1}.spm.util.imcalc.outdir = {output_dir};
    matlabbatch{1}.spm.util.imcalc.expression = 'mean(X)';
    matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
    matlabbatch{1}.spm.util.imcalc.options.mask = 0;
    matlabbatch{1}.spm.util.imcalc.options.interp = 1;
    matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

    try
        spm('defaults', 'fmri');
        spm_jobman('run', matlabbatch);
        disp(['Mean image saved as: ', fullfile(output_dir, output_filename)]);
    catch ME
        warning('Failed to run SPM job for cluster %d: %s', cluster_num, ME.message);
    end

    % Save the .nii file path
    nii_files{1, cluster_num} = fullfile(output_dir, output_filename);

    % Clear batch
    clear matlabbatch;
end

% ===== Generate BrainNet Visualizations =====
for cluster_num = unique_clusters'
    nii_file = nii_files{1, cluster_num};
    save_path = fullfile(output_dir, sprintf('brainnet_mean_%s_s%d.jpg', subtype_label, cluster_num));

    try
        BrainNet_MapCfg(surface_path, nii_file, save_path, cfig);
        disp(['BrainNet image saved as: ', save_path]);
    catch ME
        warning('Failed to generate BrainNet image for cluster %d: %s', cluster_num, ME.message);
    end
    display(save_path);
end





%% subtype and stage mean
% Parameters
subtype_label = '3x';  % Used in output file names
output_dir = fullfile(root, 'roi_inter_c2means_rainbow');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
    fprintf('Directory created: %s\n', output_dir);
else
    fprintf('Directory already exists: %s\n', output_dir);
end

% Define stage groupings
stage_groups = {
    [0:8], 
    [9:11], 
    [12:14], 
    [15:20]
};

% Get all valid subtypes (excluding 99/missing)
unique_subtypes = unique(fdata.subtype);
unique_subtypes = unique_subtypes(unique_subtypes ~= 99 & ~isnan(unique_subtypes));

% Store output paths
nii_files = containers.Map;

% ===== Generate Subtype + Stage Group Mean .nii Files =====
for subtype = unique_subtypes'
    for sg_idx = 1:length(stage_groups)
        stage_range = stage_groups{sg_idx};

        % Find subjects matching subtype and stage group
        idx = fdata.subtype == subtype & ismember(fdata.stage, stage_range);
        if ~any(idx)
            continue;
        end

        % Get image paths directly from fdata.fname
        suvr_image_paths = fdata.fname(idx);
        suvr_image_paths = cellstr(suvr_image_paths);  % enforce cell array format

        % Output filename and path
        output_filename = sprintf('mean_%s_s%d_stage%d-%d.nii', subtype_label, subtype, min(stage_range), max(stage_range));
        output_filepath = fullfile(output_dir, output_filename);

        % Build SPM ImCalc batch
        matlabbatch{1}.spm.util.imcalc.input = suvr_image_paths;
        matlabbatch{1}.spm.util.imcalc.output = output_filename;
        matlabbatch{1}.spm.util.imcalc.outdir = {output_dir};
        matlabbatch{1}.spm.util.imcalc.expression = 'mean(X)';
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 1;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

        try
            spm('defaults', 'fmri');
            spm_jobman('run', matlabbatch);
            disp(['Saved mean image: ', output_filepath]);
            nii_files(sprintf('s%d_stage%d-%d', subtype, min(stage_range), max(stage_range))) = output_filepath;
        catch ME
            warning('SPM failed for subtype %d, stage %d–%d: %s', subtype, min(stage_range), max(stage_range), ME.message);
        end

        clear matlabbatch;
    end
end

% ===== Generate BrainNet Visualizations =====
for key = keys(nii_files)
    key_str = key{1};
    nii_path = nii_files(key_str);
    save_path = fullfile(output_dir, sprintf('brainnet_mean_%s_%s.jpg', subtype_label, key_str));

    try
        BrainNet_MapCfg(surface_path, nii_path, save_path, cfig);
        disp(['BrainNet image saved: ', save_path]);
    catch ME
        warning('BrainNet failed for %s: %s', key_str, ME.message);
    end
end
