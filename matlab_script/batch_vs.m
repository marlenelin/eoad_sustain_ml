%% Batch Run Script
% remember to change the run_vs_script.m variabless
% Define models
models = {
   % 'fname ~ subtypetwo*day_to_baseline + subtypethree*day_to_baseline + CDR_SB_baseline + Gender_baseline + Age_baseline + Yrs_of_Education_baseline + Centiloids_MRI_Based_Composite_baseline + (1 + day_to_baseline|subj)', 'full_s1ref';
   % 'fname ~ subtypeone*day_to_baseline + subtypethree*day_to_baseline + CDR_SB_baseline + Gender_baseline + Age_baseline + Yrs_of_Education_baseline + Centiloids_MRI_Based_Composite_baseline + (1 + day_to_baseline|subj)', 'full_s2ref';
    %'fname ~ subtypeone*day_to_baseline + subtypetwo*day_to_baseline + CDR_SB_baseline + Gender_baseline + Age_baseline + Yrs_of_Education_baseline + Centiloids_MRI_Based_Composite_baseline + (1 + day_to_baseline|subj)', 'full_s3ref';
    %'fname ~ subtypetwo*day_to_baseline + subtypethree*day_to_baseline + CDR_SB_baseline + Gender_baseline + Age_baseline + Yrs_of_Education_baseline + Centiloids_MRI_Based_Composite_baseline + Stage_baseline + (1 + day_to_baseline|subj)', 'full_s1ref_stageadj'
 %'fname ~ subtypeone*day_to_baseline + subtypethree*day_to_baseline + CDR_SB_baseline + Gender_baseline + Age_baseline + Yrs_of_Education_baseline + Centiloids_MRI_Based_Composite_baseline + Stage_baseline + (1 + day_to_baseline|subj)', 'full_s2ref_stageadj'
%'fname ~ subtypeone*day_to_baseline + subtypetwo*day_to_baseline + CDR_SB_baseline + Gender_baseline + Age_baseline + Yrs_of_Education_baseline + Centiloids_MRI_Based_Composite_baseline + Stage_baseline + (1 + day_to_baseline|subj)', 'full_s3ref_stageadj'
%'mname ~ subtypeone*day_to_baseline + subtypetwo*day_to_baseline + cdr_sb_baseline + sex_baseline + age_at_pet_baseline + years_education_baseline + Centiloid_baseline + TIV_in_mL_baseline + (1 + day_to_baseline|subj)', 'full_s3ref_mri';
%'mname ~ subtypeone*day_to_baseline + subtypethree*day_to_baseline + cdr_sb_baseline + sex_baseline + age_at_pet_baseline + years_education_baseline + Centiloid_baseline + TIV_in_mL_baseline + (1 + day_to_baseline|subj)', 'full_s2ref_mri';
%'mname ~ subtypethree*day_to_baseline + subtypetwo*day_to_baseline + cdr_sb_baseline + sex_baseline + age_at_pet_baseline + years_education_baseline + Centiloid_baseline + TIV_in_mL_baseline + (1 + day_to_baseline|subj)', 'full_s1ref_mri';
%'fbb_file ~ subtypeone*day_to_baseline + subtypetwo*day_to_baseline + cdr_sb_baseline + sex_baseline + age_at_pet_baseline + years_education_baseline + Stage_baseline + (1 + day_to_baseline|subj)', 'stage_s3ref_amy';
%'fbb_file ~ subtypeone*day_to_baseline + subtypethree*day_to_baseline + cdr_sb_baseline + sex_baseline + age_at_pet_baseline + years_education_baseline + Stage_baseline + (1 + day_to_baseline|subj)', 'stage_s2ref_amy';
%'fbb_file ~ subtypethree*day_to_baseline + subtypetwo*day_to_baseline + cdr_sb_baseline + sex_baseline + age_at_pet_baseline + years_education_baseline + Stage_baseline + (1 + day_to_baseline|subj)', 'stage_s1ref_amy';
%'fbb_file ~ subtypeone*day_to_baseline + subtypetwo*day_to_baseline + Stage_baseline + (1 + day_to_baseline|subj)', 'no_s3ref_amy';
%'fbb_file ~ subtypeone*day_to_baseline + subtypethree*day_to_baseline + Stage_baseline + (1 + day_to_baseline|subj)', 'no_s2ref_amy';
%'fbb_file ~ subtypethree*day_to_baseline + subtypetwo*day_to_baseline + Stage_baseline + (1 + day_to_baseline|subj)', 'no_s1ref_amy';
'mname ~ subtypeone*day_to_baseline + subtypetwo*day_to_baseline + cdr_sb_baseline + sex_baseline + age_at_pet_baseline + years_education_baseline + Stage_baseline + TIV_in_mL_baseline + Centiloid_baseline +(1 + day_to_baseline|subj)', 'stage_s3ref_mri';
'mname ~ subtypeone*day_to_baseline + subtypethree*day_to_baseline + cdr_sb_baseline + sex_baseline + age_at_pet_baseline + years_education_baseline + Stage_baseline + TIV_in_mL_baseline + Centiloid_baseline + (1 + day_to_baseline|subj)', 'stage_s2ref_mri';
'mname ~ subtypethree*day_to_baseline + subtypetwo*day_to_baseline + cdr_sb_baseline + sex_baseline + age_at_pet_baseline + years_education_baseline + Stage_baseline + TIV_in_mL_baseline + Centiloid_baseline + (1 + day_to_baseline|subj)', 'stage_s1ref_mri';
'mname ~ subtypeone*day_to_baseline + subtypetwo*day_to_baseline + Stage_baseline + TIV_in_mL_baseline + Centiloid_baseline + (1 + day_to_baseline|subj)', 'no_s3ref_mri';
'mname ~ subtypeone*day_to_baseline + subtypethree*day_to_baseline + Stage_baseline + TIV_in_mL_baseline + Centiloid_baseline + (1 + day_to_baseline|subj)', 'no_s2ref_mri';
'mname ~ subtypethree*day_to_baseline + subtypetwo*day_to_baseline + Stage_baseline + TIV_in_mL_baseline + Centiloid_baseline + (1 + day_to_baseline|subj)', 'no_s1ref_mri';


    };

% Loop through models
for i = 1:size(models, 1)
    stringModel = models{i, 1};
    modelName = models{i, 2};
    
    fprintf('Running model: %s\n', modelName);
    run_vs_script(stringModel, modelName);  % Call the external script
    fprintf('Completed model: %s\n\n', modelName);
end
