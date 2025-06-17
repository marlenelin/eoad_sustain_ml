%% Define input parameters
modality = 'tau'; % Set this to your desired modality
adjust_name = 'no'; % this does not include stage
subtypes = [1, 2, 3]; % The three subtype numbers
add_mask = double(~niftiread('/home/mac/mlin2/cluster_analysis/rs16Sub_harvard_bin01.nii'));
og_mask = '/home/mac/mlin2/cluster_analysis/TPM_GM_Mask_10pc.nii';
combined_mask = niftiread(og_mask) & add_mask;
reshape_combined_mask = reshape(combined_mask,[17545,121]);
% Define base folder and output folder
base_folder = '/home/mac/mlin2/cluster_analysis/voxelstat analysis/'; % Change this if needed
output_folder = fullfile(base_folder, sprintf('%s_%s_doublep', modality, adjust_name));

% Create output folder if it doesn't exist
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
%% Loop through each subtype
for subtype = subtypes
    % Define input file names
    if strcmp(modality, "tau")
        mod_name = "";
    else
        mod_name = "_" + modality;
    end
    
    fwep_file = fullfile(base_folder, sprintf('%s_s%dref%s/corrected_%s_s%dref%s_long.mat', adjust_name, subtype, mod_name, adjust_name, subtype, mod_name));
    nofwep_file = fullfile(base_folder, sprintf('%s_s%dref%s/uncorrected_%s_s%dref%s_long.mat', adjust_name, subtype, mod_name, adjust_name, subtype, mod_name));
    
    % Check if both files exist
    if ~exist(fwep_file, 'file') || ~exist(nofwep_file, 'file')
        fprintf('Skipping subtype %d: Missing files.\n', subtype);
        continue;
    end
    
    % Load data
    fwept = load(fwep_file).corrected_struct.day_to_baseline;
    c_struct = load(nofwep_file).c_struct;
    nofwept = c_struct.tValues.day_to_baseline;
    est = c_struct.eValues.day_to_baseline * 365.25;
    
    % Apply threshold to nofwept
    nofwept(abs(nofwept) < 3.3065) = 0;
    
    % Filter nofwept where fwept has values
    filtered_nofwept = nofwept;
    filtered_nofwept(fwept ~= 0) = 0;
    filtered_nofwept(filtered_nofwept ~= 0) = 1;
    
    % Save fwept as NIfTI
    fwep_filename = fullfile(base_folder, sprintf('%s_%s_s%d_fwept.nii', modality, adjust_name, subtype));
    % use write nifti 
    
    % Save filtered_fwep 
    filtered_nofwep_filename = fullfile(base_folder, sprintf('%s_%s_s%d_filtered_nofwep.nii', modality, adjust_name, subtype));
    % use wrift nifti
    
    % Combine fwept and filtered nofwept
    combined_map = fwept + filtered_nofwept;
    % Apply og mask
    combined_map(reshape_combined_mask == 0) = 0;
    
    % Define output filename
    output_filename = fullfile(base_folder, sprintf('%s_%s_s%d_d2blest_double.nii', modality, adjust_name, subtype));
    
    % Apply to estimate image 
    est(combined_map == 1) = -est(combined_map == 1);
    est(combined_map == 0) = 0;

    % Save as NIfTI file
    %final_nii = make_nii(combined_map);
    VoxelStatsWriteNifti(est,output_filename,og_mask);
    fprintf('Saved: %s\n', output_filename);
end

disp('Processing complete.'); 
%% combine contrast 
% Define a map for numeric subtypes to their spelled-out versions
subtype_map = {'one', 'two', 'three'};

% Initialize combined maps
fwep_combined = zeros(121, 145, 121);
nofwep_combined = zeros(121, 145, 121);

for subtype = subtypes
    % Get the spelled-out subtype name (e.g., "one" for subtype 1)
    spelled_subtype = subtype_map{subtype}; 
    subtype_folder = fullfile(base_folder, sprintf('%s_s%dref%s', adjust_name, subtype, mod_name));
    
    % List all .nii files in the folder containing 'day_to_baseline' and 'subtype' in their name
    nii_files = dir(fullfile(subtype_folder, '*.nii'));
    
    % Filter files that contain both "day_to_baseline" and the spelled subtype name in any order
    fwep_files = nii_files(contains({nii_files.name}, 'day_to_baseline') & ...
                           contains({nii_files.name}, 'subtype'));
    % If no relevant files are found, skip this subtype
    if isempty(fwep_files)
        fprintf('Skipping subtype %d: No relevant fwep contrast files found.\n', subtype);
        continue;
    end
    
    % Loop over all found fwep contrast files and load them
    for fwep_idx = 1:length(fwep_files)
        % Construct full path to the fwep file
        fwep_file = fullfile(fwep_files(fwep_idx).folder, fwep_files(fwep_idx).name);
        % Load the fwep contrast NIfTI file
        contrast_fwep = load_nii(fwep_file).img;
        % Add the processed fwep contrast to the combined map
        fwep_combined = fwep_combined + abs(contrast_fwep);
    end 
    % Locate the uncorrected contrast t-map from the struct for the respective subtype
    uncorrected_file = fullfile(subtype_folder, sprintf('uncorrected_%s_s%dref%s_long.mat', adjust_name,subtype,mod_name));
    
    if ~exist(uncorrected_file, 'file')
        fprintf('Skipping subtype %d: Missing uncorrected file.\n', subtype);
        continue;
    end
     
   % Load the uncorrected contrast struct
    c_struct = load(uncorrected_file).c_struct;
    
    % Get all field names from tValues
    field_names = fieldnames(c_struct.tValues);
    
    % Define regex pattern to match any order of "day_to_baseline" and the spelled subtype
    pattern = sprintf('^(day_to_baseline.*subtype|subtype.*day_to_baseline)$');
    
    % Find all matching fields
    matching_field_idxs = find(~cellfun('isempty', regexp(field_names, pattern, 'once')));
    
    if isempty(matching_field_idxs)
        warning('No matching field found for subtype: %s', spelled_subtype);
        continue; % Skip if no match is found
    end
    
    % Loop through all matching fields and add them to nofwep_combined
    for i = 1:length(matching_field_idxs)
        matching_field = field_names{matching_field_idxs(i)};
        contrast_nofwep = c_struct.tValues.(matching_field);
        contrast_nofwep(abs(contrast_nofwep) < 3.3065) = 0;
        % Accumulate values into nofwep_combined
        nofwep_combined = nofwep_combined + reshape(abs(contrast_nofwep), [121, 145, 121]);
    end
  
end

% Process the fwep contrast (set non-zero values to 3)
fwep_combined(fwep_combined ~= 0) = 3;

% Process the nofwep contrast (set non-zero values to 1)
nofwep_combined(nofwep_combined ~= 0) = 1;

% Set nofwep values where fwep is not zero to 0 (i.e., filter out where fwep is non-zero)
nofwep_combined(fwep_combined ~= 0) = 0;
   
% Save the final combined contrast maps
final_combined = fwep_combined + nofwep_combined;
final_combined(combined_mask == 0) = 0;
final_output_filename = fullfile(base_folder, sprintf('%s_%s_doublep_combined.nii', modality, adjust_name));
VoxelStatsWriteNifti(final_combined, final_output_filename,og_mask);
fprintf('Saved: %s\n', final_output_filename);

disp('Processing complete.');


%% visualize
% combine contrast
surface = 'ICBM'; 
view = 'lm';
surface_path = '/home/mac/mlin2/Downloads/BrainNetViewer_20191031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv'; % BrainMesh_Ch2_smoothed.nv';
cfig = 'brainnet config/doublep_binarize.mat'; % Configuration file for BrainNet viewer
folder_path = '/home/mac/mlin2/cluster_analysis/voxelstat analysis/'; % Folder to search within
disp('generating images');
save_path = fullfile(base_folder, sprintf('%s_%s_doublep_combined.jpg', modality, adjust_name));
% Get list of all nii files in the folder and subfolders that contain 'combined_contrast' in their names
BrainNet_MapCfg(surface_path, final_output_filename, save_path, cfig);

% baseline estimate
bl_fig = "brainnet config/pt12_vi.mat";
%"brainnet config/mri_d2blest_pt25.mat";
%""brainnet config/amy_d2blest_pos.mat";
for subtype = subtypes
    blest_path = fullfile(base_folder, sprintf('%s_%s_s%d_d2blest_double.nii', modality, adjust_name, subtype));
    blest_spath = fullfile(base_folder, sprintf('%s_%s_s%d_d2blest_double.jpg', modality, adjust_name, subtype));
    BrainNet_MapCfg(surface_path, blest_path, blest_spath, bl_fig);
end



