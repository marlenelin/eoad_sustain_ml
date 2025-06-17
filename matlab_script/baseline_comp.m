%% basic setting
% Correct usage of fullfile to construct paths
root = '/home/mac/mlin2/cluster_analysis';
filename = '/type assignment/sense_pt8_threshold.csv';%inter_c2mean_ss_fname_bl
% Configuration for visualization paths
surface = 'ICBM';
surface_path = fullfile('/home/mac/mlin2/Downloads/BrainNetViewer_20191031/Data/SurfTemplate', 'BrainMesh_ICBM152_smoothed.nv');
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


%% 1vr 
warning('on');
subtype_column = 'subtype';
type = 'tau_inter_c2mean_4s_nofwep';
root_comparison_dir = '/home/mac/mlin2/cluster_analysis/comparison 1vr/';
column_name = '4x';
disp('Reading data...');
fbb_file = fdata.fname;
% Get unique subtypes excluding missing values (99)
unique_clusters = unique(fdata.(subtype_column)(fdata.(subtype_column) ~= 99));

for i = 1:length(unique_clusters)
        cluster_num1 = unique_clusters(i);
        cluster_indices1 = fdata.(subtype_column) == cluster_num1;
        cluster_indices_rest = fdata.(subtype_column) ~= cluster_num1 & fdata.(subtype_column) ~= 99;
        
        % Extract subject lists
        subjects_in_cluster1 = fdata.subj(cluster_indices1);
        subjects_in_cluster_rest = fdata.subj(cluster_indices_rest);
        
        % Extract and combine covariatess
        covariate_names = {'sex', 'years_education', 'cdr_sb_baseline', 'stage_baseline', 'age_at_pet', 'TIV_in_mL','Centiloid_baseline'};
        for cov_idx = 1:length(covariate_names)
            cov_name = covariate_names{cov_idx};
            combined_covariate = [fdata.(cov_name)(cluster_indices1); fdata.(cov_name)(cluster_indices_rest)];   
            eval(['combined_', cov_name, ' = combined_covariate;']);
        end
        
        % Define a specific folder for this 1 vs. rest comparison
        comparison_folder = fullfile(root_comparison_dir, sprintf('%s_ttest_%s_s%d', type, subtype_column, cluster_num1));
        if ~exist(comparison_folder, 'dir')
            mkdir(comparison_folder);
        else
            disp(['Folder already exists: ', comparison_folder, '. Skipping calculation.']);
            continue;
        end
        
        % Prepare image paths
        suvr_image_paths_cluster1 = fbb_file(cluster_indices1);
        suvr_image_paths_cluster_rest = fbb_file(cluster_indices_rest);
        
        cd('/home/mac/mlin2/cluster_analysis/');
        % SPM two-sample t-test factorial design
        clear matlabbatch;

        % Factorial design for two-sample t-test
        matlabbatch{1}.spm.stats.factorial_design.dir = {comparison_folder};
        matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = suvr_image_paths_cluster1;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = suvr_image_paths_cluster_rest;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
 
        % Covariate
        stage = 1;
        demo = 1;
        cov_weights = [1 -1];
        if stage && demo  
            matlabbatch{1}.spm.stats.factorial_design.cov(1).c = combined_stage_baseline;
            matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'stage';
            matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
            matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 5; 
            matlabbatch{1}.spm.stats.factorial_design.cov(2).c = combined_sex;
            matlabbatch{1}.spm.stats.factorial_design.cov(2).cname = 'sex';
            matlabbatch{1}.spm.stats.factorial_design.cov(2).iCFI = 1;
            matlabbatch{1}.spm.stats.factorial_design.cov(2).iCC = 5; 
            matlabbatch{1}.spm.stats.factorial_design.cov(3).c = combined_age_at_pet;
            matlabbatch{1}.spm.stats.factorial_design.cov(3).cname = 'age';
            matlabbatch{1}.spm.stats.factorial_design.cov(3).iCFI = 1;
            matlabbatch{1}.spm.stats.factorial_design.cov(3).iCC = 5; 
            matlabbatch{1}.spm.stats.factorial_design.cov(4).c = combined_years_education;
            matlabbatch{1}.spm.stats.factorial_design.cov(4).cname = 'educ';
            matlabbatch{1}.spm.stats.factorial_design.cov(4).iCFI = 1;
            matlabbatch{1}.spm.stats.factorial_design.cov(4).iCC = 5; 
            matlabbatch{1}.spm.stats.factorial_design.cov(5).c = combined_cdr_sb_baseline;
            matlabbatch{1}.spm.stats.factorial_design.cov(5).cname = 'cdrsb';
            matlabbatch{1}.spm.stats.factorial_design.cov(5).iCFI = 1;
            matlabbatch{1}.spm.stats.factorial_design.cov(5).iCC = 5; 
            matlabbatch{1}.spm.stats.factorial_design.cov(6).c = combined_Centiloid_baseline;
            matlabbatch{1}.spm.stats.factorial_design.cov(6).cname = 'centiloid';
            matlabbatch{1}.spm.stats.factorial_design.cov(6).iCFI = 1;
            matlabbatch{1}.spm.stats.factorial_design.cov(6).iCC = 5; 
           % matlabbatch{1}.spm.stats.factorial_design.cov(7).c = combined_TIV_in_mL;
           %  matlabbatch{1}.spm.stats.factorial_design.cov(7).cname = 'TIV';
           %  matlabbatch{1}.spm.stats.factorial_design.cov(7).iCFI = 1;
           %  matlabbatch{1}.spm.stats.factorial_design.cov(7).iCC = 5; 
            
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            cov_weights = [1 -1 0 0 0 0 0 0]; %0];
        elseif stage
            matlabbatch{1}.spm.stats.factorial_design.cov(1).c = combined_stage_baseline;
            matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'stage';
            matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
            matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 5; 
            %matlabbatch{1}.spm.stats.factorial_design.cov.c = combined_TIV_in_mL;
            % matlabbatch{1}.spm.stats.factorial_design.cov.cname = 'TIV';
            % matlabbatch{1}.spm.stats.factorial_design.cov.iCFI = 1;
            % matlabbatch{1}.spm.stats.factorial_design.cov.iCC = 5; 
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            cov_weights = [1 -1 0];
        else
            matlabbatch{1}.spm.stats.factorial_design.cov = struct([]);
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        end
        % No masking
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};

        % Global calculation settings
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

        % Model estimation
        matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(comparison_folder, 'SPM.mat')};
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

        % Contrast specification (Cluster 1 > Rest)
        matlabbatch{3}.spm.stats.con.spmmat = {fullfile(comparison_folder, 'SPM.mat')};
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'higher w cov';
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = cov_weights;
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.delete = 0;

        % masking for contrast
        matlabbatch{3}.spm.stats.results.conspec.mask.image.name = {'TPM_GM_Mask_10pc.nii,1'};
        matlabbatch{3}.spm.stats.results.conspec.mask.image.mtype = 0;

        % Run the SPM job
        spm_jobman('run', matlabbatch);

        disp(['Two-sample t-test for ', sprintf('%s_ttest_%s_s%d_vs_rest', type, column_name, cluster_num1), ' completed.']);
        cd('/home/mac/mlin2/cluster_analysis/');
        % T-threshold calculation and TabDat extraction
        clear matlabbatch;
        % Step 1: Clear existing TabDat
        if exist('TabDat', 'var')
            clear TabDat;
        end
        spmmat_path = fullfile(comparison_folder, 'SPM.mat');
        if ~exist(spmmat_path, 'file')
            error('SPM.mat file not found in %s', comparison_folder);
        end
        
        fwep = 0;
        fwe_config = 'none';
        if fwep
            fwe_config = 'FWE';
            pth = 0.05;
        else
            pth = 0.001;
        end
        % Set up and run the results batch to extract T values
        matlabbatch{1}.spm.stats.results.spmmat = {spmmat_path};
        matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
        matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
        matlabbatch{1}.spm.stats.results.conspec.threshdesc = fwe_config;%'FWE';
        matlabbatch{1}.spm.stats.results.conspec.thresh = pth;
        matlabbatch{1}.spm.stats.results.conspec.extent = 20;
        matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
        matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
        matlabbatch{1}.spm.stats.results.units = 1;
        
        % Run the results job to generate TabDat
        spm_jobman('run', matlabbatch);
        
        % Access and save the FWE-corrected T-value threshold
        if fwep
            fwep_value = TabDat.ftr{5, 2}(1);
            disp(['FWE threshold for ', sprintf('%s_ttest_%s_s%d_vs_rest', type, column_name, cluster_num1), ': ', num2str(fwep_value)]);
        else
            fwep_value = TabDat.ftr{1,2}(1);
            disp(['Unadjusted threshold for ', sprintf('%s_ttest_%s_s%d_vs_rest', type, column_name, cluster_num1), ': ', num2str(fwep_value)]);

        end
        % Save TabDat
        save(fullfile(comparison_folder, 'TabDat_data.mat'), 'TabDat');
        
        % Thresholding the NIfTI file based on the T critical value
        nii = spm_vol(fullfile(comparison_folder, 'spmT_0001.nii'));
        [data, ~] = spm_read_vols(nii);

        % Apply the FWE-corrected T-value threshold
        threshold = fwep_value;

        % Retain positive and negative values above/below threshold
        positive_data = data;
        positive_data(positive_data < threshold) = 0;

        negative_data = data;
        negative_data(negative_data > -threshold) = 0;
        negative_data(negative_data >= 0) = 0; % Keep only negative values

        % Combine the positive and negative contrasts
        combined_data = positive_data + negative_data;

        % Save the combined data as a new NIfTI file
        nii_combined = nii;
        nii_combined.fname = fullfile(comparison_folder, 'combined_contrast.nii');
        spm_write_vol(nii_combined, combined_data);
        cd('/home/mac/mlin2/cluster_analysis/');
end 
%% visualization
% Specify constants
surface = 'ICBM'; 
view = 'lm';
surface_path = '/home/mac/mlin2/Downloads/BrainNetViewer_20191031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv'; % BrainMesh_Ch2_smoothed.nv';
cfig = 'brainnet config/tmap_afni.mat'; % Configuration file for BrainNet viewer
folder_path = '/home/mac/mlin2/cluster_analysis/comparison 1vr/'; % Folder to search within
disp('generating images');
% Get list of all nii files in the folder and subfolders that contain 'combined_contrast' in their names
nii_files = dir(fullfile(folder_path, '*subtype*', '*combined_contrast.nii'));

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





%% range check
% Load the NIfTI file
if 0
    nii = niftiread('/home/mac/mlin2/cluster_analysis/voxelstat analysis/full_s3ref_amy/subtypeoneday_to_baseline.nii');
    nii = nii(nii ~= 0);
    maxVal = max(nii(:));
    minVal = min(nii(:)); % Minimum intensity value
    disp(['Min: ', num2str(minVal), ', Max: ', num2str(maxVal)]);
    p5 = prctile(nii(:), 5); % 5th percentile
    p95 = prctile(nii(:), 95); % 95th percentile
    disp(['5%: ', num2str(p5), ', 95%: ', num2str(p95)]);
    % Convert to double (if needed)
    nii = double(nii(:)); 
    
    % Plot histogram
    figure;
    histogram(nii, 100); % Use 100 bins for better granularity
    xlabel('Intensity Value');
    ylabel('Frequency');
    title('Histogram of NIfTI Image Values');
    grid on;
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
