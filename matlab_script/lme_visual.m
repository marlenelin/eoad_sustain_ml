%%
subtype = {1, 2, 3};
mask_name = 'rs16Sub_harvard_bin01.nii';
root_dir = '/home/mac/mlin2/cluster_analysis/';
og_mask = '/home/mac/mlin2/cluster_analysis/TPM_GM_Mask_10pc.nii';
suffix = 'ref_amy';
%%
% Loop through the subtypes
for subtype_num = 1:length(subtype)
    % Define paths
    t_path = fullfile(root_dir, 'voxelstat analysis/',...
        ['full_s',num2str(subtype_num),suffix], '/day_to_baseline.nii');
    est_path = fullfile(root_dir, 'voxelstat analysis/',...
        ['full_s',num2str(subtype_num),suffix], '/day_to_baseline_est.nii');
    mask_path = fullfile(root_dir, mask_name);
    c_struct_path = fullfile(root_dir, 'voxelstat analysis/',...
        ['full_s',num2str(subtype_num),suffix],  ...
        ['/uncorrected_full_s',num2str(subtype_num),suffix,'_long.mat']);
    % Load c_struct and extract day_to_baseline estimate
    if exist(c_struct_path, 'file')
        data = load(c_struct_path, 'c_struct'); % Assumes the MAT file contains 'c_struct'
        day2bl_est = data.c_struct.eValues.day_to_baseline;
        day2bl_est = day2bl_est * 365.25; % Cast day2bl_est to match mask data type
        % Save day_to_baseline estimate as a NIfTI file
        VoxelStatsWriteNifti(day2bl_est,est_path,og_mask);
        %niftiwrite(day2bl_est, day2bl_path, mask_info); % Use mask_info for consistent metadata
        disp(['NIfTI file saved from c_struct to: ', est_path]);
    else
        error(['c_struct file not found: ', c_struct_path]);
    end

    % Read the NIfTI files
    t_info = niftiinfo(t_path); % Metadata for the day_to_baseline file
    t_data = niftiread(t_path); % Read the day_to_baseline data
    mask_data = niftiread(mask_path); % Read the mask data
    est_info = niftiinfo(est_path); % Metadata for the est image
    est_data = niftiread(est_path); % Read the est image data

    % Ensure dimensions match between mask and data
    if ~isequal(size(t_data), size(mask_data)) || ~isequal(size(est_data), size(mask_data))
        error('Dimensions do not match between day2bl, mask, and est files.');
    end

    % Convert data types to double for compatibility
    t_data = double(t_data);
    mask_data = double(~mask_data); 
    % Create a combined mask by multiplying day_to_baseline and harvard mask
    combined_mask = t_data & mask_data;
    % Apply the combined mask to the est image
    est_data = double(est_data); % Ensure est_data is also double
    est_data(combined_mask == 0) = 0;

    % Convert est_data back to the original data type before saving
    est_data = cast(est_data, est_info.Datatype);

    % Save the modified est image
    niftiwrite(est_data, est_path, est_info);
    disp(['Modified est image saved to: ', est_path]);
end
%%

warning('on');
subtype = {1, 2, 3}; % List of subtypes 
surface_path = fullfile('/home/mac/mlin2/Downloads/BrainNetViewer_20191031/Data/SurfTemplate', 'BrainMesh_ICBM152_smoothed.nv');
cfig = fullfile('/home/mac/mlin2/cluster_analysis/brainnet config', 'mri_d2blest.mat');


% Loop through the subtypes
for subtype_num = 1:length(subtype)
    % Define paths
    est_path = fullfile(root_dir, 'voxelstat analysis/', ...
        ['full_s',num2str(subtype_num),suffix], '/day_to_baseline_est.nii');
    estimg_path = fullfile(root_dir, 'voxelstat analysis/',...
        ['full_s',num2str(subtype_num),suffix], '/day_to_baseline_est.jpg');
   
    % Check if the required NIfTI file exists
    if ~isfile(est_path)
        warning(['NIfTI file does not exist for subtype ', num2str(subtype{subtype_num}), ': ', est_path]);
        continue;
    end

    % Generate BrainNet visualization
    try
        BrainNet_MapCfg(surface_path, est_path, estimg_path, cfig);
        disp(['BrainNet image saved as: ', estimg_path]);
    catch ME
        warning(['Error generating BrainNet visualization for subtype ', num2str(subtype{subtype_num}), ': ', ME.message]);
    end
end

%% combined interaction map
% Define subtype names 
subtypen = {'one', 'two', 'three'};
cpfig = fullfile('/home/mac/mlin2/cluster_analysis/brainnet config', 'tmap_afni.mat');
% Root directory where files are located
root_dir = '/home/mac/mlin2/cluster_analysis/';

% Loop through subtypes (outer loop defines the current subtype of interest)
for subtype_num = 1:length(subtypen)
        nii_paths = {};
        for subtype_comp = 1:length(subtypen)
            if subtype_comp ~= subtype_num
                ref_nii = fullfile(root_dir, ...
                    ['voxelstat analysis/full_s', num2str(subtype_comp), suffix], ...
                    ['/subtype', subtypen{subtype_num}, 'day_to_baseline.nii']);
                display(ref_nii);
                nii_paths{end+1} = ref_nii;
            end
        end
        % Load the two NIfTI files
        nii1 = load_nii(nii_paths{1});
        nii2 = load_nii(nii_paths{2});
        nii1.img(nii1.img > 0) = 3;    % Set positive values to 1
        nii1.img(nii1.img < 0) = -3;   % Set negative values to -1
        nii1.img(nii1.img == 0) = 0;   % Ensure zero values remain 0
        
        nii2.img(nii2.img > 0) = 6;    % Set positive values to 1
        nii2.img(nii2.img < 0) = -6;   % Set negative values to -1
        nii2.img(nii2.img == 0) = 0;   % Ensure zero values remain 0

        combined_nii = nii1;
        combined_nii.img = nii1.img + nii2.img;
        
        % Save the combined NIfTI file
        combined_nii_path = fullfile(root_dir, ...
            ['voxelstat analysis/output/'], ...
            ['combined_',suffix,  subtypen{subtype_num}, '.nii']);
        save_nii(combined_nii, combined_nii_path);
        
        % Define the output path for the visualization
        savep = fullfile(root_dir, ...
            ['voxelstat analysis/output/'], ...
            ['combined_', suffix,subtypen{subtype_num}, '.jpg']);
        
        % Visualize the combined NIfTI file
        try
            BrainNet_MapCfg(surface_path, combined_nii_path, savep, cpfig);
            disp(['Visualization saved for ', subtypen{subtype_comp}, ' vs ', subtypen{subtype_num}, ' at ', savep]);
        catch ME
            warning('Failed to visualize NIfTI file: %s' );
        end
end
 

%% joint binarize
% Define folder path containing the NIfTI files
nii_folder = '/home/mac/mlin2/cluster_analysis/voxelstat analysis/output/';

% Define NIfTI file names
nii_files = {
    fullfile(nii_folder, ['combined_',suffix,'one.nii']), ...
    fullfile(nii_folder, ['combined_',suffix,'two.nii']), ...
    fullfile(nii_folder, ['combined_',suffix,'three.nii'])
};

% Define output files
output_bin_nii = fullfile(nii_folder, ['combined_',suffix,'_binary.nii']); % Output binary NIfTI
output_brainnet_image = fullfile(nii_folder, ['combined_',suffix,'_brainnet.jpg']); % BrainNet image

% Path to BrainNet Viewer surface and config files (adjust if needed)
surface_path = fullfile('/home/mac/mlin2/Downloads/BrainNetViewer_20191031/Data/SurfTemplate', 'BrainMesh_ICBM152_smoothed.nv');
cfig = fullfile('/home/mac/mlin2/cluster_analysis/brainnet config', 'combined_darkred.mat');


% Initialize empty variable for sum
nii_data_sum = [];

% Process each NIfTI file
for i = 1:length(nii_files)
    % Load the NIfTI file
    nii_ori = spm_vol(nii_files{i});
    nii_data = spm_read_vols(nii_ori);

    % Take absolute value
    nii_data = abs(nii_data);

    % Sum voxel values across all images
    if isempty(nii_data_sum)
        nii_data_sum = nii_data;
    else
        nii_data_sum = nii_data_sum + nii_data;
    end
end

% Binarize: Set all nonzero values to 1
nii_data_bin = nii_data_sum > 0;
% Apply filter 
mask_data = niftiread(fullfile(root_dir, mask_name)); % Read the mask data
mask_data = double(~mask_data); 
% Create a combined mask by multiplying day_to_baseline and harvard mask
nii_data_bin = nii_data_bin & mask_data;
% Save the binarized NIfTI file
nii_bin = nii_ori; % Copy header from original
nii_bin.fname = output_bin_nii;
spm_write_vol(nii_bin, nii_data_bin);

% Visualize with BrainNet Viewer
BrainNet_MapCfg(surface_path, output_bin_nii, output_brainnet_image, cfig);

% Display confirmation message
disp(['Binarized NIfTI saved as: ', output_bin_nii]);
disp(['BrainNet visualization saved as: ', output_brainnet_image]);
