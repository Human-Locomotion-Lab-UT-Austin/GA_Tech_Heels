%% Batch analysis across participants, visits and footwear conditions
% Loops analyze_participant.m over every participant in the parameter file for
% each visit (Pre = DayA_01, Post = DayA_03) and footwear condition (HEEL,
% FLAT), and saves one results table per visit x condition:
%
%   data/processed/batch_results_{heels,flats}_{Day1,Day3}.csv
%
% Run from anywhere: paths are set relative to this file's location in the
% repository. The raw .mat files must be in data/raw (see README.md).

%% Paths (relative to the repository root)
repo_dir      = fileparts(fileparts(mfilename('fullpath')));
raw_dir       = fullfile(repo_dir, "data", "raw");
params_dir    = fullfile(repo_dir, "data", "params");
processed_dir = fullfile(repo_dir, "data", "processed");
if ~exist(processed_dir, 'dir'), mkdir(processed_dir); end

%% Load the three master data files (each contains every participant
%% as a top-level field, e.g. MECH_Data.SS11, MECH_Data.SE15, ...)
MECH_Data        = importdata(fullfile(raw_dir, "MECH_Data.mat"));
JointMomentData  = importdata(fullfile(raw_dir, "JointMomentData.mat"));
JointAngleData   = importdata(fullfile(raw_dir, "JointAngleData.mat"));

%% Visits and conditions
% day: label in the data structs; params: parameter file for that visit;
% tag: suffix of the output file name
visits = struct( ...
    'day',    {'DayA_01', 'DayA_03'}, ...
    'params', {"participant_params_Day1.csv", "participant_params_Day3.csv"}, ...
    'tag',    {"Day1", "Day3"});
conditions = struct('cond', {'HEEL', 'FLAT'}, 'tag', {"heels", "flats"});

result_fields = {'E', 'mean_peak_Fmtu_EMA', 'mean_pushoff_GRF_EMA', 'mean_lin_strain_impulse', ...
    'mean_peak_lin_strain', 'mean_lin_strain', 'mean_peak_Fr', ...
    'mean_peak_Fmtu', 'mean_Fmtu', 'mean_peak_ank_mom', 'mean_ank_mom'};

for v = 1:numel(visits)
    % Expected columns: ParticipantID, lin_k, rest_AtMA, CSA, Lo.
    % ParticipantID must exactly match the struct field name in the .mat files.
    param_table = readtable(fullfile(params_dir, visits(v).params), 'TextType', 'string');
    n_participants = height(param_table);

    for c = 1:numel(conditions)
        results_table = table(param_table.ParticipantID, 'VariableNames', {'ParticipantID'});
        for f = 1:numel(result_fields)
            results_table.(result_fields{f}) = nan(n_participants, 1);
        end

        for k = 1:n_participants
            pid = param_table.ParticipantID(k);
            try
                r = analyze_participant(pid, visits(v).day, conditions(c).cond, ...
                        param_table.lin_k(k), param_table.rest_AtMA(k), ...
                        param_table.CSA(k), param_table.Lo(k), ...
                        MECH_Data, JointMomentData, JointAngleData);
                for f = 1:numel(result_fields)
                    results_table.(result_fields{f})(k) = r.(result_fields{f});
                end
            catch ME
                % Missing or unusable data: warn and leave this row as NaN
                warning('run_batch_analysis:participantFailed', ...
                        'Participant %s (%s %s) failed: %s', pid, ...
                        visits(v).day, conditions(c).cond, ME.message);
            end
        end

        out_file = fullfile(processed_dir, sprintf('batch_results_%s_%s.csv', ...
                            conditions(c).tag, visits(v).tag));
        writetable(results_table, out_file);
        fprintf('Saved results for %d participants to %s\n', n_participants, out_file);
    end
end
