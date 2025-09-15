%% Basic settings
warning('on');
% Paths
root = '/home/mac/mlin2/cluster_analysis';
filename = '/type assignment/final_inter_c2mean_2_bl_ss_with_fpath_covariates.csv'; % inter_c2mean_ss_fname_bl
surface = 'ICBM';
surface_path = fullfile('/home/mac/mlin2/Downloads/BrainNetViewer_20191031/Data/SurfTemplate', 'BrainMesh_ICBM152_smoothed.nv');
cfig = fullfile(root, 'brainnet config', 'cet_rbgyrm1_pt5_3.mat');

% Read data
fdata_path = fullfile(root, filename);
fdata0 = readtable(fdata_path, 'ReadVariableNames', true);
disp('read');

% Normalize sex to 0/1 if needed
if iscell(fdata0.sex) || isstring(fdata0.sex) || iscategorical(fdata0.sex)
    sx = string(fdata0.sex);
    fdata0.sex = double(ismember(lower(sx), ["male", "m", "1"])); % males=1 else 0
end
% Toggles
% Adjustment 
fwep = 0; % 0: uncorrected p-threshold (0.001), 1: FWE 0.05
% Covariates
stage = 1;
demo = 1;
% Modalities
file = fdata0.fbb_file;  % paths fpath for tau

% Loop over i = 2:6
for i = 3
    fprintf('\n=== Processing i = %d ===\n', i);
    subtype_column = sprintf('x%d', i);
    type = sprintf('amy_fullnocdr_%ds_fwep', i);
    root_comparison_dir = '/home/mac/mlin2/cluster_analysis/comparison 1vr/';
    column_name = 'subtype';

    % Work on a copy per-iteration
    fdata = fdata0;
    stage_colname = sprintf('x%ds', i);
    % Derive stage column (assumes baseline stage is in x5s; keep as in your code)
    if ismember(stage_colname, fdata.Properties.VariableNames)
        fdata.stage_baseline = fdata.(stage_colname);
    else
        error('Expected column x s for stage_baseline is missing.');
    end

    % Required variables for this iteration (include current subtype col)
    req_vars = {'sex','years_education','age_at_pet','Centiloid_baseline','fpath', subtype_column};
    missing_mask = false(height(fdata),1);
    for rv = 1:numel(req_vars)
        if ~ismember(req_vars{rv}, fdata.Properties.VariableNames)
            error('Required variable "%s" not found in the table.', req_vars{rv});
        end
    end
    fdata = rmmissing(fdata, 'DataVariables', req_vars);

    % Unique subtypes, excluding 99
    if ~ismember(subtype_column, fdata.Properties.VariableNames)
        error('Subtype column %s not found in table after rmmissing.', subtype_column);
    end
    subvals = fdata.(subtype_column);
    unique_clusters = unique(subvals(subvals ~= 99));

    % Iterate through each subtype: 1 vs rest
    for k = 1:length(unique_clusters)
        cluster_num1 = unique_clusters(k);
        fprintf('Subtype %s == %d: two-sample t-test vs rest\n', subtype_column, cluster_num1);

        cluster_indices1   = fdata.(subtype_column) == cluster_num1;
        cluster_indices_rest = fdata.(subtype_column) ~= cluster_num1 & fdata.(subtype_column) ~= 99;

        % Extract and combine covariates
        covariate_names = {'sex', 'years_education', 'TIV_in_mL_baseline', 'stage_baseline', 'age_at_pet', 'Centiloid_baseline'};
        for cov_idx = 1:length(covariate_names)
            cov_name = covariate_names{cov_idx};
            combined_covariate = [fdata.(cov_name)(cluster_indices1); fdata.(cov_name)(cluster_indices_rest)];
            eval(['combined_', cov_name, ' = combined_covariate;']);
        end

        % Define per-comparison output folder
        comparison_folder = fullfile(root_comparison_dir, sprintf('%s_ttest_%s_s%d', type, subtype_column, cluster_num1));
        if ~exist(comparison_folder, 'dir')
            mkdir(comparison_folder);
        else
            disp(['Folder already exists: ', comparison_folder, '. Skipping calculation.']);
            continue;
        end

        % Prepare image paths (cellstr for SPM)
        suvr_image_paths_cluster1   = cellstr(file(cluster_indices1));
        suvr_image_paths_cluster_rest = cellstr(file(cluster_indices_rest));

        % Move to working dir
        cd('/home/mac/mlin2/cluster_analysis/');

        % --- SPM two-sample t-test design ---
        clear matlabbatch;
        matlabbatch{1}.spm.stats.factorial_design.dir = {comparison_folder};
        matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = suvr_image_paths_cluster1;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = suvr_image_paths_cluster_rest;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
        
        cov_weights = [1 -1]; % base for group difference
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

            %matlabbatch{1}.spm.stats.factorial_design.cov(5).c = combined_Centiloid_baseline;
            %matlabbatch{1}.spm.stats.factorial_design.cov(5).cname = 'centiloid';
            %matlabbatch{1}.spm.stats.factorial_design.cov(5).iCFI = 1;
            %matlabbatch{1}.spm.stats.factorial_design.cov(5).iCC = 5;

            %matlabbatch{1}.spm.stats.factorial_design.cov(6).c = combined_TIV_in_mL_baseline;
            %matlabbatch{1}.spm.stats.factorial_design.cov(6).cname = 'tiv';
            %matlabbatch{1}.spm.stats.factorial_design.cov(6).iCFI = 1;
            %matlabbatch{1}.spm.stats.factorial_design.cov(6).iCC = 5;


            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            cov_weights = [1 -1 0 0 0 0]; % pad zeros for covariates
        elseif stage
            matlabbatch{1}.spm.stats.factorial_design.cov(1).c = combined_stage_baseline;
            matlabbatch{1}.spm.stats.factorial_design.cov(1).cname = 'stage';
            matlabbatch{1}.spm.stats.factorial_design.cov(1).iCFI = 1;
            matlabbatch{1}.spm.stats.factorial_design.cov(1).iCC = 5;
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
            cov_weights = [1 -1 0];
        else
            matlabbatch{1}.spm.stats.factorial_design.cov = struct([]);
            matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        end

        % Masking and globals (unchanged)
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

        % Model estimation
        matlabbatch{2}.spm.stats.fmri_est.spmmat = {fullfile(comparison_folder, 'SPM.mat')};
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

        % Contrast: cluster > rest (with covariates)
        matlabbatch{3}.spm.stats.con.spmmat = {fullfile(comparison_folder, 'SPM.mat')};
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'higher w cov';
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = cov_weights;
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        matlabbatch{3}.spm.stats.con.delete = 0;

        % Results mask for contrast
        matlabbatch{3}.spm.stats.results.conspec.mask.image.name = {'TPM_GM_Mask_10pc.nii,1'};
        matlabbatch{3}.spm.stats.results.conspec.mask.image.mtype = 0;

        % Run the SPM job
        spm_jobman('run', matlabbatch);
        disp(['Two-sample t-test for ', sprintf('%s_ttest_%s_s%d_vs_rest', type, column_name, cluster_num1), ' completed.']);

        % --- Results / TabDat / thresholding ---
        cd('/home/mac/mlin2/cluster_analysis/');
        clear matlabbatch;

        if exist('TabDat', 'var'); clear TabDat; end
        spmmat_path = fullfile(comparison_folder, 'SPM.mat');
        if ~exist(spmmat_path, 'file')
            error('SPM.mat file not found in %s', comparison_folder);
        end
   
        fwe_config = 'none';
        if fwep
            fwe_config = 'FWE'; pth = 0.05;
        else
            pth = 0.001;
        end

        % Create results to populate TabDat
        matlabbatch{1}.spm.stats.results.spmmat = {spmmat_path};
        matlabbatch{1}.spm.stats.results.conspec.titlestr = '';
        matlabbatch{1}.spm.stats.results.conspec.contrasts = 1;
        matlabbatch{1}.spm.stats.results.conspec.threshdesc = fwe_config;
        matlabbatch{1}.spm.stats.results.conspec.thresh = pth;
        matlabbatch{1}.spm.stats.results.conspec.extent = 20;
        matlabbatch{1}.spm.stats.results.conspec.conjunction = 1;
        matlabbatch{1}.spm.stats.results.conspec.mask.none = 1;
        matlabbatch{1}.spm.stats.results.units = 1;

        spm_jobman('run', matlabbatch);

        % Extract threshold from TabDat
        if fwep
            fwep_value = TabDat.ftr{5, 2}(1);
            disp(['FWE threshold for ', sprintf('%s_ttest_%s_s%d_vs_rest', type, column_name, cluster_num1), ': ', num2str(fwep_value)]);
        else
            fwep_value = TabDat.ftr{1, 2}(1);
            disp(['Unadjusted threshold for ', sprintf('%s_ttest_%s_s%d_vs_rest', type, column_name, cluster_num1), ': ', num2str(fwep_value)]);
        end

        % Save TabDat for record
        save(fullfile(comparison_folder, 'TabDat_data.mat'), 'TabDat');

        % Threshold NIfTI (keep pos >= thr; neg <= -thr)
        nii = spm_vol(fullfile(comparison_folder, 'spmT_0001.nii'));
        [data, ~] = spm_read_vols(nii);
        threshold = fwep_value;

        positive_data = data;
        positive_data(positive_data < threshold) = 0;

        negative_data = data;
        negative_data(negative_data > -threshold) = 0;
        negative_data(negative_data >= 0) = 0;

        combined_data = positive_data + negative_data;

        nii_combined = nii;
        nii_combined.fname = fullfile(comparison_folder, 'combined_contrast.nii');
        spm_write_vol(nii_combined, combined_data);

        cd('/home/mac/mlin2/cluster_analysis/');
    end % end unique_clusters loop
end % end i loop

%% Visualization
% Constants
surface = 'ICBM';
view = 'lm';
surface_path = '/home/mac/mlin2/Downloads/BrainNetViewer_20191031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv';
cfig = 'brainnet config/tmap_afni.mat';
folder_path = '/home/mac/mlin2/cluster_analysis/comparison 1vr/';
disp('generating images');

% Find all combined_contrast.nii files in subtype subfolders
nii_files = dir(fullfile(folder_path, '*subtype*', '*combined_contrast.nii'));

for i = 1:length(nii_files)
    [~, target, ~] = fileparts(nii_files(i).name);
    fpath = fullfile(nii_files(i).folder, [target, '.nii']);
    save_path = fullfile(nii_files(i).folder, sprintf('%s_%s_%s.jpg', surface, view, target));
    BrainNet_MapCfg(surface_path, fpath, save_path, cfig);
end
