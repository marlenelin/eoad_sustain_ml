%% Define input parameters
modality = 'amy'; % Set this to your desired modality
adjust_name = '3s'; % Set this, or change as needed
subtypes = [1, 2, 3]; % The three subtype numbers

% Define base folder and output folder
base_folder = '/home/mac/mlin2/cluster_analysis/comparison 1vr/'; % Change this if needed
output_folder = fullfile(base_folder, sprintf('%s_fullnocdr_%s', modality, adjust_name));

% Create output folder if it doesn't exist
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

% Loop through each subtype
for subtype = subtypes
    % Define input file names
    fwep_file = fullfile(base_folder, sprintf('%s_fullnocdr_%s_fwep_ttest_x3_s%d/combined_contrast.nii', modality, adjust_name, subtype));
    nofwep_file = fullfile(base_folder, sprintf('%s_fullnocdr_%s_nofwep_ttest_x3_s%d/combined_contrast.nii', modality, adjust_name, subtype));
    
    % Check if both files exist
    if ~exist(fwep_file, 'file') || ~exist(nofwep_file, 'file')
        fprintf('Skipping subtype %d: Missing files.\n', subtype);
        continue;
    end
    
    % Load NIfTI images
    fwep_nii = niftiread(fwep_file);
    nofwep_nii = niftiread(nofwep_file);
    nii_info = niftiinfo(fwep_file); % Use fwep as reference for output
    
    % Step 1: Binarize noFWEp file
    binarized_nofwep = zeros(size(nofwep_nii));
    binarized_nofwep(nofwep_nii > 0) = 1;
    binarized_nofwep(nofwep_nii < 0) = -1;
    
    % Step 2: Filter out values in binarized noFWEp where FWEp has nonzero values
    binarized_nofwep(fwep_nii ~= 0) = 0;
    
    % Step 3: Add filtered noFWEp to FWEp
    combined_p = fwep_nii + binarized_nofwep;
    
    % Define output file name
    output_file = fullfile(output_folder, sprintf('%s_%s_doublep_s%d.nii', modality, adjust_name, subtype));
    
    % Check if the file already exists and delete it
    if exist(output_file, 'file')
        delete(output_file);
    end
    
    % Save combined NIfTI file
    niftiwrite(combined_p, output_file, nii_info);
    
    fprintf('Processed subtype %d and saved to %s\n', subtype, output_file);
end

disp('Processing complete.');

%% visualize
surface = 'ICBM'; 
view = 'lm';
surface_path = '/home/mac/mlin2/Downloads/utilities/BrainNetViewer_20191031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv'; % BrainMesh_Ch2_smoothed.nv';
cfig = 'brainnet config/doublep_thre3.mat'; % Configuration file for BrainNet viewer
folder_path = '/home/mac/mlin2/cluster_analysis/comparison 1vr/'; % Folder to search within
disp('generating images');
% Get list of all nii files in the folder and subfolders that contain 'combined_contrast' in their names
nii_files = dir(fullfile(output_folder, "*.nii"));

% Iterate over each nii fil]\
for i = 1:length(nii_files)
    % Get the target file name without the extension
    [~, target, ~] = fileparts(nii_files(i).name);
    
    % Construct the paths
    fpath = fullfile(nii_files(i).folder, [target, '.nii']);  % Path to the NIfTI file
    save_path = fullfile(nii_files(i).folder, sprintf('%s_%s_%s.jpg', surface, view, target));  % Path for the output image
    
    % Call BrainNet_MapCfg with the constructed paths
    BrainNet_MapCfg(surface_path, fpath, save_path, cfig);
end




