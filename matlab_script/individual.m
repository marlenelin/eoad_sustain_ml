% Specify the scans folder
scansFolder = 'scans/';

% Get all NIfTI files in the scans folder
niiFiles = dir(fullfile(scansFolder, '*.nii'));

% Initialize list to hold files to smooth
filesToSmooth = {};

% Check for unsmoothed files
for i = 1:length(niiFiles)
    % Get file name and path
    [~, target, ext] = fileparts(niiFiles(i).name);
    
    % Skip already smoothed files (starting with 's')
    if startsWith(target, 's')
        continue;
    end
    
    % Check if smoothed file exists
    smoothedFile = fullfile(scansFolder, ['s', target, ext]);
   % if ~isfile(smoothedFile)
        % Add to the list for smoothing
    filesToSmooth{end+1} = fullfile(scansFolder, niiFiles(i).name);
   % end
end

% Smoothing using SPM
if ~isempty(filesToSmooth)
    % Initialize SPM batch
    matlabbatch = {};
    matlabbatch{1}.spm.spatial.smooth.data = filesToSmooth';
    matlabbatch{1}.spm.spatial.smooth.fwhm = [5.3 5.3 5.3]; % Smoothing kernel
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';

    % Run SPM batch
    spm('defaults', 'FMRI');
    spm_jobman('run', matlabbatch);
end

% Generate images from smoothed files
smoothedFiles = dir(fullfile(scansFolder, 's*.nii'));

% Constants for BrainNet
surface = 'ICBMs';
view = 'lm_jet';
surface_path = '/home/mac/mlin2/Downloads/BrainNetViewer_20191031/Data/SurfTemplate/BrainMesh_ICBM152_smoothed.nv'; % Path to the surface file
%cfig = fullfile('/home/mac/mlin2/cluster_analysis/brainnet config', 'lateral_medial_0pt5_3pt5_jet.mat'); % Configuration file
cfig = fullfile('/home/mac/mlin2/cluster_analysis/brainnet config','lm_0pt5_3pt5_jet.mat');
% Iterate over smoothed files
for i = 1:length(smoothedFiles)
    % Get the target (file name without the extension)
    [~, target, ~] = fileparts(smoothedFiles(i).name);
    
    % Construct the paths
    fpath = fullfile(scansFolder, smoothedFiles(i).name);  % Path to the NIfTI file
    save_path = fullfile(scansFolder, sprintf('%s_%s_%s.jpg', surface, view, target));  % Path for the output image
    
    % Check if the output image already exists
    if isfile(save_path)
        fprintf('Skipping %s: Image already exists.\n', smoothedFiles(i).name);
        continue;
    end
    
    % Print out current progress
    fprintf('Processing image %d out of %d: %s\n', i, length(smoothedFiles), smoothedFiles(i).name);
    
    % Call BrainNet_MapCfg with the updated configuration file
    BrainNet_MapCfg(surface_path, fpath, save_path, cfig);
end
