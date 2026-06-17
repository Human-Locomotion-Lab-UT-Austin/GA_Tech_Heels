%% Batch analysis across participants
% Loops analyze_participant.m over every participant listed in your
% parameter spreadsheet and collects the 6 summary outputs into one table.

%% --- EDIT THESE PATHS ---
data_dir   = "/Users/andrewthornton/Documents/WCB2026/Data";
param_file = fullfile(data_dir, "participant_params_heels_Day3.csv");
% -------------------------

%% Load the three master data files (each contains every participant
%% as a top-level field, e.g. MECH_Data.SS11, MECH_Data.SS12, ...)
MECH_Data        = importdata(fullfile(data_dir, "MECH_Data.mat"));
JointMomentData  = importdata(fullfile(data_dir, "JointMomentData.mat"));
JointAngleData   = importdata(fullfile(data_dir, "JointAngleData.mat"));

%% Load subject-specific parameters
% Expected columns: ParticipantID, lin_k, rest_AtMA, CSA, Lo
% ParticipantID must exactly match the struct field name in the .mat
% files (e.g. "SS11").
param_table = readtable(param_file, 'TextType', 'string');

n_participants = height(param_table);

%% Loop through participants
results_struct = repmat(struct( ...
    'ParticipantID', "", ...
    'E', NaN, ...
    'mean_peak_Fmtu_EMA', NaN, ...
    'mean_lin_strain_impulse', NaN, ...
    'mean_peak_lin_strain', NaN, ...
    'mean_lin_strain', NaN, ...
    'mean_peak_Fr', NaN), n_participants, 1);

for k = 1:n_participants
    pid = param_table.ParticipantID(k);

    try
        r = analyze_participant(pid, ...
                param_table.lin_k(k), ...
                param_table.rest_AtMA(k), param_table.CSA(k), param_table.Lo(k), ...
                MECH_Data, JointMomentData, JointAngleData);

        results_struct(k).ParticipantID          = pid;
        results_struct(k).E                       = r.E;
        results_struct(k).mean_peak_Fmtu_EMA      = r.mean_peak_Fmtu_EMA;
        results_struct(k).mean_lin_strain_impulse = r.mean_lin_strain_impulse;
        results_struct(k).mean_peak_lin_strain    = r.mean_peak_lin_strain;
        results_struct(k).mean_lin_strain         = r.mean_lin_strain;
        results_struct(k).mean_peak_Fr             = r.mean_peak_Fr;

    catch ME
        warning('run_batch_analysis:participantFailed', ...
                'Participant %s failed: %s', pid, ME.message);
        results_struct(k).ParticipantID = pid; % numeric fields stay NaN
    end
end

%% Assemble into a table and save
results_table = struct2table(results_struct);
disp(results_table);

writetable(results_table, fullfile(data_dir, 'batch_results_heels_Day3.csv'));
fprintf('Saved results for %d participants to %s\n', n_participants, ...
    fullfile(data_dir, 'batch_results_heels_Day3.csv'));
