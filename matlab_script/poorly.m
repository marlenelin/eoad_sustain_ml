%% Configuration
% Load data from the CSV file
filename = 'type assignment/roi_inter_c2mean_mtl_bl_ss.csv';%'intersection_c2mean_max_e6_typedf.csv';  % Change to your file if needed
fdata = readtable(filename, 'ReadVariableNames', true);

% Path to the directory where the SUVR images are stored
suvr_image_folder = 'scans/';  % Change this to the path where your SUVR images are located

% Subtype columns to iterate through
subtype_columns = {'x2','x3','x4','x5','x6'};  % Subtype columns in fdata with x99
type = 'roi'; 
  
% Paths for BrainNet Viewer visualization
surface_path = '/home/mac/mlin2/Downloads/BrainNetViewer_20191031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';  % Surface template path
cfig = fullfile('brainnet config','lm_0_1pt75_jet.mat');  % Change this to your specific configuration file

%% Iterate through each subtype column
for col_idx = 1:length(subtype_columns)
    column_name = subtype_columns{col_idx};
    
    % Directory for the new folder to copy files into
    new_folder = fullfile(root_comparison_dir, [column_name, '_value99_files']);  % New folder name

    % Create the new folder if it doesn't exist
    if ~exist(new_folder, 'dir')
        mkdir(new_folder);
        disp(['Created folder: ', new_folder]);
    end

    % Find all rows where the value in the current column is 99
    value_99_indices = fdata.(column_name) == 99;

    % Extract the corresponding file names (assuming fdata.fname holds file names)
    file_names = fdata.fname(value_99_indices);

    % Loop through the files associated with value 99 and copy them
    for i = 1:length(file_names)
        fname = file_names{i};
        
        % Subset the file name to exclude the last four characters
        fname_clean = fname(1:end-4);  % Removes the last 4 characters (e.g., '.nii')
        
        % Construct the NIfTI file path
        nii_file = fullfile(suvr_image_folder, ['s', fname_clean, '.nii']);
        
        % Construct the JPG file path using the first 13 characters of the NIfTI file name
        jpg_pattern = fullfile(suvr_image_folder, ['*', fname_clean(3:12), '*.jpg']);  % Use wildcard to match any JPG file that starts with the first 13 characters
        
        % Search for matching JPG files
        jpg_files = dir(jpg_pattern);
        
        % Copy the NIfTI file if it exists
        if exist(nii_file, 'file')
            copyfile(nii_file, fullfile(new_folder, ['s', fname_clean, '.nii']));
            disp(['Copied NIfTI file: ', nii_file]);
        else
            disp(['NIfTI file not found: ', nii_file]);
        end
        
        % Copy any JPG file that matches the pattern
        if ~isempty(jpg_files)
            for j = 1:length(jpg_files)
                jpg_file = fullfile(jpg_files(j).folder, jpg_files(j).name);
                copyfile(jpg_file, fullfile(new_folder, jpg_files(j).name));
                disp(['Copied JPG file: ', jpg_file]);
            end
        else
            disp(['No matching JPG files found for pattern: ', jpg_pattern]);
        end
    end

    % Visualization using BrainNet Viewer
    % Path to the customizable config file
    
    % Get all NIfTI files in the new folder for visualization
    nii_files = dir(fullfile(new_folder, '*.nii'));

    % Iterate over each NIfTI file for visualization
    for k = 1:length(nii_files)
        % Construct the paths
        fpath = fullfile(nii_files(k).folder, nii_files(k).name);  % Path to the NIfTI file
        save_path = fullfile(new_folder, sprintf('%s.jpg', nii_files(k).name));  % Path for the output image
        
        % Call BrainNet_MapCfg with the constructed paths
        BrainNet_MapCfg(surface_path, fpath, save_path, cfig);
        disp(['Generated BrainNet image for: ', fpath]);
    end
end
