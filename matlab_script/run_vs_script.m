function run_vs_script(stringModel, modelName)
    %% LME Analysis Script

    % Input paths
    imageType = 'nifti';
    data_file = '/home/mac/mlin2/cluster_analysis/type assignment/full_combined_v3_amy.csv';
    mask_file = '/home/mac/mlin2/cluster_analysis/TPM_GM_Mask_10pc.nii';
    multivalueVariables = {'mname'};%mname %fname %'fbb_file
    categoricalVars = {''}; % none
    includeString = '';

    % Output directory
    output_path = fullfile('/home/mac/mlin2/cluster_analysis/voxelstat analysis', modelName);
    if ~exist(output_path, 'dir')
        mkdir(output_path);
    end
    save_fname = fullfile(output_path, ['uncorrected_' modelName '_long_mri.mat']);

    % Run LME
    [c_struct, slices_p, image_height_p, image_width_p, coeff_vars, voxel_num, df, voxel_dims] = ...
        VoxelStatsLME(imageType, stringModel, data_file, mask_file, multivalueVariables, categoricalVars, includeString);

    % Save uncorrected results
    save(save_fname, 'c_struct');

    %% Multiple Comparison Correction
    fieldn = fieldnames(c_struct.tValues);
    image_dims = [slices_p, image_height_p, image_width_p];
    search_vol = voxel_num * prod(voxel_dims);
    num_voxels = voxel_num;
    peak_pval = 0.05;
    clus_th = 0.001;
    fwhm = 20;
    corrected_struct = struct();

    for f = 1:length(fieldn)
        field = fieldn{f};
        stat_mat = eval(['c_struct.tValues.' field]);
        corrected_stats_mat = VoxelStatsDoRFT(stat_mat, image_dims, search_vol, num_voxels, fwhm, df, peak_pval, clus_th);
        corrected_struct.(field) = corrected_stats_mat;

        % Save as NIfTI
        nii_file_path = fullfile(output_path, [field, '.nii']);
        VoxelStatsWriteNifti(corrected_stats_mat, nii_file_path, mask_file);
        fprintf('Saved %s to %s\n', field, nii_file_path);
    end

    % Save corrected results
    save_corrected_fname = fullfile(output_path, ['corrected_' modelName '_long_mri.mat']);
    save(save_corrected_fname, 'corrected_struct');

    %% Visualization with BrainNet
    surface_path = fullfile('/home/mac/mlin2/Downloads/BrainNetViewer_20191031/Data/SurfTemplate', 'BrainMesh_ICBM152_smoothed.nv');
    cfig = fullfile('/home/mac/mlin2/cluster_analysis/brainnet_config', 'tmap_afni.mat');
    nii_files = dir(fullfile(output_path, '*.nii'));

    for j = 1:length(nii_files)
        nii_file = fullfile(output_path, nii_files(j).name);
        [~, file_name, ~] = fileparts(nii_files(j).name);
        save_path = fullfile(output_path, [file_name, '.jpg']);

        try
            BrainNet_MapCfg(surface_path, nii_file, save_path, cfig);
            disp(['BrainNet image saved as: ', save_path]);
        catch ME
            warning('Failed to generate BrainNet image for %s: %s', nii_file, ME.message);
        end
    end
end
