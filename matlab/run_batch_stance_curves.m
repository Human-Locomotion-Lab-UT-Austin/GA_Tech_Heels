%% Batch export of mean stance-phase curves for every participant, visit and condition
% Runs export_stance_curves.m for both visits (Pre = DayA_01 with the Day1
% parameters, Post = DayA_03 with the Day3 parameters) and both footwear
% conditions, and writes one long-format CSV (AT strain, resultant GRF, AT
% force and EMA over stance) that WCB2026_Stats.R averages across participants:
%
%   data/processed/stance_curves.csv
%
% Run from anywhere: paths are set relative to this file's location in the
% repository. The raw .mat files must be in data/raw (see README.md).

%% Paths (relative to the repository root)
repo_dir      = fileparts(fileparts(mfilename('fullpath')));
raw_dir       = fullfile(repo_dir, "data", "raw");
params_dir    = fullfile(repo_dir, "data", "params");
processed_dir = fullfile(repo_dir, "data", "processed");
if ~exist(processed_dir, 'dir'), mkdir(processed_dir); end
out_file = fullfile(processed_dir, "stance_curves.csv");

%% Load the master data files
MECH_Data        = importdata(fullfile(raw_dir, "MECH_Data.mat"));
JointMomentData  = importdata(fullfile(raw_dir, "JointMomentData.mat"));
JointAngleData   = importdata(fullfile(raw_dir, "JointAngleData.mat"));

% Visit label in the data structs, matching parameter file, and name used in R
visits = struct( ...
    'day',        {'DayA_01', 'DayA_03'}, ...
    'param_file', {"participant_params_Day1.csv", "participant_params_Day3.csv"}, ...
    'label',      {"Pre", "Post"});
conditions = ["HEEL", "FLAT"];

rows = {};
for v = 1:numel(visits)
    params = readtable(fullfile(params_dir, visits(v).param_file), 'TextType', 'string');
    for c = conditions
        for k = 1:height(params)
            pid = params.ParticipantID(k);
            try
                curve = export_stance_curves(pid, visits(v).day, c, ...
                    params.lin_k(k), params.rest_AtMA(k), params.Lo(k), ...
                    MECH_Data, JointMomentData, JointAngleData);
                rows{end+1} = table( ...
                    repmat(pid, 101, 1), repmat(visits(v).label, 101, 1), repmat(c, 101, 1), ...
                    curve.stance_pct, curve.strain_mean, curve.strain_sd, curve.grf_bw, curve.fmtu_bw, ...
                    curve.r_internal, curve.R_external, curve.ema_median, ...
                    repmat(curve.n_stances, 101, 1), ...
                    'VariableNames', {'ParticipantID', 'Visit', 'Condition', 'stance_pct', ...
                                      'strain_mean', 'strain_sd', 'grf_bw', 'fmtu_bw', 'r_internal', ...
                                      'R_external', 'ema_median', 'n_stances'}); %#ok<SAGROW>
            catch ME
                warning('run_batch_stance_curves:failed', '%s %s %s failed: %s', ...
                        pid, visits(v).label, c, ME.message);
            end
        end
    end
end

curves = vertcat(rows{:});
writetable(curves, out_file);
fprintf('Saved %d curves to %s\n', height(curves) / 101, out_file);
