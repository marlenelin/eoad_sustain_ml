%% data and config 
filename = '/home/mac/mlin2/cluster_analysis/fixed_ss_exclude_bl.csv';%fixed_ % type and stage csv file
fdata = readtable(filename, 'ReadVariableNames', true); 
suvr_image_folder = 'scans/baseline';   
subtype_columns = {'x3'};  %, 'x5', 'x6'
stage_columns = {'x3s'}; 
type = 'fixed'; % threshold type (used in naming)
root_comparison_dir = '/home/mac/mlin2/cluster_analysis/comparison pw';  % Root folder for comparisons
%% run comparison (only  nifti)
% number of split
for col_idx = 1:length(subtype_columns)
    column_name = subtype_columns{col_idx};
    % Get the unique cluster numbers (subtypes) in the current column
    unique_clusters = unique(fdata.(column_name));
    % Exclude subtype 99 - poorly fit with subtype prob < 0.5
    unique_clusters = unique_clusters(unique_clusters ~= 99);
    % Perform pairwise comparisons between each cluster
    for i = 1:length(unique_clusters)
        for j = i+1:length(unique_clusters)
            cluster_num1 = unique_clusters(i);
            cluster_num2 = unique_clusters(j);
            % Get the list of subjects in cluster 1
            cluster_indices1 = fdata.(column_name) == cluster_num1;
            file_names_in_cluster1 = fdata.fname(cluster_indices1);
            % Get the list of subjects in cluster 2
            cluster_indices2 = fdata.(column_name) == cluster_num2;
            file_names_in_cluster2 = fdata.fname(cluster_indices2);
            % Collect the corresponding SUVR image file paths for cluster 1
            suvr_image_paths_cluster1 = {};
            for k = 1:length(file_names_in_cluster1)
                fname = file_names_in_cluster1{k};
                suvr_image_paths_cluster1{end+1, 1} = fullfile(suvr_image_folder, strcat('s', fname));  % Column vector
            end
            % Collect the corresponding SUVR image file paths for cluster 2
            suvr_image_paths_cluster2 = {};
            for k = 1:length(file_names_in_cluster2)
                fname = file_names_in_cluster2{k};
                suvr_image_paths_cluster2{end+1, 1} = fullfile(suvr_image_folder, strcat('s', fname));  % Column vector
            end
            % Define a specific folder for this pairwise comparison
            comparison_folder = fullfile(root_comparison_dir, sprintf('%s_ttest_%s_s%dv%d', type, column_name, cluster_num1, cluster_num2));
            % Skip if folder already exists, otherwise create it
            if ~exist(comparison_folder, 'dir')
                mkdir(comparison_folder);  % Create the folder if it doesn't exist
            else
                disp(['Folder already exists: ', comparison_folder, '. Skipping calculation.']);
                continue;
            end
            clear matlabbatch;
            % Factorial design for two-sample t-test
            matlabbatch{1}.spm.stats.factorial_design.dir = {comparison_folder};  % Output directory
            matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = suvr_image_paths_cluster1;  % Subtype 1
            matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = suvr_image_paths_cluster2;  % Subtype 2
            matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
            matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
            matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
            matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
            % No masking
            matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
            matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
            % Global calculation settings
            matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
            matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

            % Model estimation
            matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(comparison_folder, 'SPM.mat')}; %what's the spm for estimation?
            matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
            matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

            % Contrast specification
            matlabbatch{3}.spm.stats.con.spmmat = {fullfile(comparison_folder, 'SPM.mat')};
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'higher';
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
            matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
            matlabbatch{3}.spm.stats.con.delete = 0;

            % Run the SPM job
            spm_jobman('run', matlabbatch);
            disp(['Two-sample t-test for ', sprintf('%s_ttest_%s_s%dv%d', type, column_name, cluster_num1, cluster_num2), ' completed.']);
         
            % T-threshold calculation and TabDat extraction
            % Clear existing TabDat and spm batch
            if exist('TabDat', 'var')
                clear TabDat;
            end
            clear matlabbatch;
            cd('/home/mac/mlin2/cluster_analysis/');
            % Set up and run the results batch to extract T values
            matlabbatch{1}.spm.stats.results.spmmat = {fullfile(comparison_folder, 'SPM.mat')};
            matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
            matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
            matlabbatch{1}.spm.stats.results.conspec.threshdesc = 'FWE';
            matlabbatch{1}.spm.stats.results.conspec.thresh = 0.05;
            matlabbatch{1}.spm.stats.results.conspec.extent = 20;
            matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
            matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
            matlabbatch{1}.spm.stats.results.units = 1;
            
            % Run the results job to generate TabDat
            spm_jobman('run', matlabbatch);
            
            % Access and save the FWE-corrected T-value threshold
            fwep_value = TabDat.ftr{5, 2}(1);
            disp(['FWEp value for ', sprintf('%s_ttest_%s_s%dv%d', type, column_name, cluster_num1, cluster_num2), ': ', num2str(fwep_value)]);    
            % Save TabDat
            save(fullfile(comparison_folder, 'TabDat_data.mat'), 'TabDat');
            % Thresholding the NIfTI file based on the T critical value
            nii = spm_vol(fullfile(comparison_folder, 'spmT_0001.nii'));
            [data, ~] = spm_read_vols(nii);
            % Apply the FWE-corrected T-value threshold
            threshold = fwep_value;
            % Retain values above/below pos/neg threshold, rest = 0
            positive_data = data;
            positive_data(positive_data < threshold) = 0;
            negative_data = data;
            negative_data(negative_data > -threshold) = 0;
            negative_data(negative_data >= 0) = 0;  
            % Combine the positive and negative contrasts
            combined_data = positive_data + negative_data;
            % Save the combined data as a new NIfTI file
            nii_combined = nii;
            nii_combined.fname = fullfile(comparison_folder, 'combined_contrast.nii');
            spm_write_vol(nii_combined, combined_data);  
            cd('/home/mac/mlin2/cluster_analysis/');
        end
        cd('/home/mac/mlin2/cluster_analysis/');
    end
end

%% visualize h/l combined (brainnet) separately
% Specify constants
surface = 'ICBM';
view = 'full';
surface_path = '/home/mac/mlin2/Downloads/BrainNetViewer_20191031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv'; % BrainMesh_Ch2_smoothed.nv';
cfig = 'tmap_jet.mat'; % Configuration file for BrainNet viewer
folder_path = 'comparison pw'; % Folder to search within

% Get list of all nii files in the folder and subfolders that contain 'combined_contrast' in their names
nii_files = dir(fullfile(folder_path, '**', '*combined_contrast.nii'));

% Iterate over each nii file
for i = 1:length(nii_files)
    % Get the target file name without the extension
    [~, target, ~] = fileparts(nii_files(i).name);
    
    % Construct the paths
    fpath = fullfile(nii_files(i).folder, [target, '.nii']);  % Path to the NIfTI file
    save_path = fullfile(nii_files(i).folder, sprintf('%s_%s_%s.jpg', surface, view, target));  % Path for the output image
    
    % Call BrainNet_MapCfg with the constructed paths
    BrainNet_MapCfg(surface_path, fpath, save_path, cfig);
end




%% combined pw (for 3x only)
% Define subtype names and paths
subtypen = {'one', 'two', 'three'};
cpfig = fullfile('/home/mac/mlin2/cluster_analysis/brainnet config', 'tmap_afni.mat');
root_dir = '/home/mac/mlin2/cluster_analysis';  % Root directory
surface_path = '/home/mac/mlin2/Downloads/BrainNetViewer_20191031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';  % BrainNet surface path

% Iterate over each subtype for combined contrast
for subtype_num = 1:length(subtypen)
    % Define the paths for the two NIfTI files to combine
    ref1_nii = fullfile(root_dir, ...
        'voxelstat analysis/', ['subtype*day (ref1)/'], ...
        ['subtype', subtypen{subtype_num}, 'day_to_baseline.nii']);
    
    ref2_nii = fullfile(root_dir, ...
        'voxelstat analysis/', ['subtype*day (ref2)/'], ...
        ['subtype', subtypen{subtype_num}, 'day_to_baseline.nii']);
    
    % Check if both files exist
    if ~isfile(ref1_nii) || ~isfile(ref2_nii)
        warning('Missing NIfTI file(s) for subtype: %s', subtypen{subtype_num});
        continue;
    end
    
    % Load the two NIfTI files
    nii1 = spm_vol(ref1_nii);
    nii2 = spm_vol(ref2_nii);
    [data1, ~] = spm_read_vols(nii1);
    [data2, ~] = spm_read_vols(nii2);
    
    % Combine the data using your specified logic
    combined_data = data1 + data2 * 3;  % Adjust scaling for visibility
    
    % Save the combined NIfTI file
    combined_nii = nii1;  % Use the header from the first file
    combined_nii_fname = fullfile(root_dir, ...
        ['voxelstat analysis/output/'], ...
        ['combined_subtype_', subtypen{subtype_num}, '.nii']);
    spm_write_vol(combined_nii, combined_data);
    
    % Visualize the combined NIfTI file
    savep = fullfile(root_dir, ...
        ['voxelstat analysis/output/'], ...
        ['combined_subtype_', subtypen{subtype_num}, '_view.jpg']);
    
    try
        BrainNet_MapCfg(surface_path, combined_nii_fname, savep, cpfig);
        disp(['Visualization saved for subtype: ', subtypen{subtype_num}, ' at ', savep]);
    catch ME
        warning('Failed to visualize NIfTI file for subtype %s: %s', subtypen{subtype_num}, ME.message);
    end
end
