%% LME
% Input
imageType = 'nifti';
stringModel = 'fname ~ subtypetwo*day_to_baseline + subtypeone*day_to_baseline  + (1 + day_to_baseline| subj)';
data_file = '/home/mac/mlin2/cluster_analysis/type assignment/lme_long.csv';
mask_file = '/home/mac/mlin2/cluster_analysis/TPM_GM_Mask_10pc.nii';
multivalueVariables = {'fname'};
categoricalVars = {''}; %none
save_fname = 'voxelstat analysis/uncorrected_subtype_3ref_long.mat';
includeString = '';
%optional: multiVarOperationMap  
% Run LME
[c_struct, slices_p, image_height_p, image_width_p, coeff_vars, voxel_num, df, voxel_dims] = ...
    VoxelStatsLME(imageType, stringModel, data_file, mask_file, multivalueVariables,categoricalVars,includeString);
% Save uncorrected tvalueMap (and eValues, seValues)
if 1
    save(save_fname, 'c_struct');
end
%% correct t-value map
% Input
fieldn = fieldnames(c_struct.tValues);
image_dims(1) = slices_p;
image_dims(2) = image_height_p;
image_dims(3) = image_width_p;
search_vol = voxel_num*prod(voxel_dims);
num_voxels = voxel_num;
peak_pval = 0.05; 
clus_th = 0.001; %default specified in stat_threshold
% To specify
fwhm = 20;
corrected_struct = struct();
output_path = '/home/mac/mlin2/cluster_analysis/voxelstat analysis/subtype*day long (ref3)/';
if ~exist(output_path, 'dir')
    mkdir(output_path);
end
% Run MCC
for f =1:length(fieldn) % iterate through variable
    field = fieldn{f};
    stat_mat = eval(['c_struct.tValues.' field]);
    [corrected_stats_mat] = VoxelStatsDoRFT(stat_mat, image_dims, ...
    search_vol, num_voxels, fwhm, df, peak_pval,clus_th);
    corrected_struct.(field) = corrected_stats_mat;
    nii_file_path = fullfile(output_path, [field, '.nii']);
    VoxelStatsWriteNifti(corrected_stats_mat,nii_file_path,mask_file);
    fprintf('Saved %s to %s\n', field, nii_file_path);
end
% Save correcte tValues struct as mat
if 1
    save_corrected_fname = fullfile(output_path,'corrected_subtype_long_ref3.mat');
    save(save_corrected_fname, 'corrected_struct');
end 
%% Visualize with BrainNet
surface_path = fullfile('/home/mac/mlin2/Downloads/BrainNetViewer_20191031/Data/SurfTemplate', 'BrainMesh_ICBM152_smoothed.nv');
cfig = fullfile('/home/mac/mlin2/cluster_analysis/brainnet config', 'tmap_afni.mat');
nii_files = dir(fullfile(output_path, '*.nii'));
for i = 1:length(nii_files) 
    nii_file = fullfile(output_path, nii_files(i).name);
    [~, file_name, ~] = fileparts(nii_files(i).name);  
    save_path = fullfile(output_path, [file_name, '.jpg']);
    try
        BrainNet_MapCfg(surface_path, nii_file, save_path, cfig);
        disp(['BrainNet image saved as: ', save_path]);
    catch ME
        warning('Failed to generate BrainNet image for %s: %s', nii_file, ME.message);
    end
end
%%
if 0
    % Step 1: Read the CSV file into a table
    input_file = '/home/mac/mlin2/cluster_analysis/type assignment/trial_lm_s1.csv';
    data = readtable(input_file);
    
    % Step 2: Convert the 'subtype' column to string
    if ismember('subtype', data.Properties.VariableNames)
        data.subtype = string(data.subtype);  % Convert 'subtype' column to string
    else
        error('The column "subtype" does not exist in the CSV file.');
    end
end
 