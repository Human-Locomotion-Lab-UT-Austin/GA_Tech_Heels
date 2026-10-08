function curve = export_stance_curves(pid, day, cond, lin_k, rest_AtMA, Lo, ...
                                      MECH_Data, JointMomentData, JointAngleData)
% EXPORT_STANCE_CURVES  Mean stance-phase time series for one participant,
% visit and footwear condition: Achilles tendon linear strain (%), resultant
% GRF and AT force (body weights), and the parts of effective mechanical
% advantage (EMA).
%
%   curve = export_stance_curves(pid, day, cond, lin_k, rest_AtMA, Lo, ...
%                                MECH_Data, JointMomentData, JointAngleData)
%
% Uses the same processing as analyze_participant.m and Force_vs_EMA.m
% (filters, stance detection, body mass from all strides, moment arm scaling,
% strain, EMA), so the curves match the summary values used in the
% statistics. Every complete stance in the trial is time-normalized to 0-100%
% of stance (0 = heel strike, 100 = toe-off, 101 points), then the stances are
% averaged.
%
% INPUTS
%   pid       - participant field name, e.g. "SS11"
%   day       - visit label in the data structs: 'DayA_01' (Pre) or 'DayA_03' (Post)
%   cond      - 'HEEL' or 'FLAT'
%   lin_k     - linear tendon stiffness for that visit (N/mm)
%   rest_AtMA - resting Achilles tendon moment arm (m)
%   Lo        - tendon slack length for that visit (mm)
%   MECH_Data, JointMomentData, JointAngleData - the three loaded .mat structs
%
% OUTPUT
%   curve - struct of 101x1 vectors (mean across stances unless noted):
%     stance_pct     0-100% of stance
%     strain_mean    AT linear strain (%), and strain_sd (SD across stances)
%     grf_bw         resultant GRF / body weight
%     fmtu_bw        AT (muscle-tendon unit) force / body weight
%     r_internal     AT (internal) moment arm (m)
%     R_external     GRF (external) moment arm = ankle moment / resultant GRF (m)
%     ema_median     median across stances of r_internal / R_external
%   and the scalar n_stances.
%   EMA is also returned as its parts because it is undefined where the ankle
%   moment crosses zero (R_external -> 0), so a plain mean of EMA is unstable.

pid  = char(pid);
day  = char(day);
cond = char(cond);
trial = 'S13';

fs = 1000;        % Hz, rate of the moment and angle data
fn = fs/2;        % Nyquist frequency
fc = 10;          % GRF filter cutoff (Hz)
fc_mom = 10;      % moment filter cutoff (Hz)
threshold = 10;   % N, stance/swing threshold

%% Ground reaction forces
Fy_full = -(MECH_Data.(pid).(day).(cond).(trial).allData.Force_Fy1); % N
Fz_full = -(MECH_Data.(pid).(day).(cond).(trial).allData.Force_Fz1); % N
norm_ank_mom = (JointMomentData.(pid).(day).(cond).AnkMom'); % Nm/kg, 1000 Hz

% Some GRF recordings were sampled at a multiple of 1000 Hz (SE15 Post: 2000 Hz,
% 30000 samples for the same 15 s as the 15000-sample moment data). Filter at
% the native rate, then keep every grf_ratio-th sample so the GRF lines up
% sample-for-sample with the 1000 Hz moment and angle data.
grf_ratio = round(length(Fz_full) / length(norm_ank_mom));
n = min(15000, floor(length(Fz_full) / grf_ratio));

[b, a] = butter(4, fc / (fs * grf_ratio / 2));
filt_Fz = filtfilt(b, a, Fz_full(1:n * grf_ratio));
filt_Fy = filtfilt(b, a, Fy_full(1:n * grf_ratio));
filt_Fz = filt_Fz(1:grf_ratio:end);
filt_Fy = filt_Fy(1:grf_ratio:end);

%% Heel strike / toe off detection
heelstrike_index = [];
toeoff_index = [];
for i = 2:n
    if filt_Fz(i-1) < threshold && filt_Fz(i) > threshold
        heelstrike_index(end+1) = i;
    elseif filt_Fz(i-1) > threshold && filt_Fz(i) < threshold
        toeoff_index(end+1) = i;
    end
end

if isempty(heelstrike_index) || isempty(toeoff_index)
    error('export_stance_curves:noSteps', 'No stance phases detected for %s %s %s.', pid, day, cond);
end

% Remove partial stance phases at the start and end of the trial
if toeoff_index(1) < heelstrike_index(1)
    toeoff_index = toeoff_index(2:end);
end
if toeoff_index(end) < heelstrike_index(end)
    heelstrike_index = heelstrike_index(1:end-1);
end

num_strides = min(length(heelstrike_index), length(toeoff_index));
heelstrike_index = heelstrike_index(1:num_strides);
toeoff_index = toeoff_index(1:num_strides);

if num_strides < 2
    error('export_stance_curves:tooFewSteps', 'Fewer than 2 stances for %s %s %s.', pid, day, cond);
end

%% Body weight / mass from the mean vertical GRF over all strides
body_weight_stride = zeros(1, num_strides-1);
for i = 1:(num_strides-1)
    body_weight_stride(i) = mean(filt_Fz(heelstrike_index(i):heelstrike_index(i+1))) * 2; % N
end
body_weight = mean(body_weight_stride); % N
mass = body_weight / 9.81; % kg

filt_Fr = sqrt(filt_Fz.^2 + filt_Fy.^2); % N, resultant GRF

%% Ankle moment and moment arms
[b, a] = butter(4, fc_mom/fn);
filt_norm_ank_mom = filtfilt(b, a, norm_ank_mom(1:n));
filt_ank_mom = filt_norm_ank_mom * mass; % Nm

ankle_angles = (JointAngleData.(pid).(day).(cond).AnkAng');
ankle_angles = ankle_angles(1:n);

maganaris_ankle_angles = [-15, 0, 15, 30];     % deg
maganaris_moment_arms  = [4.3, 4.7, 5.2, 5.6]; % cm, resting values
poly_coeffs = polyfit(maganaris_ankle_angles, maganaris_moment_arms, 2);
scaling_factor = rest_AtMA / polyval(poly_coeffs, 0);
internal_moment_arms = polyval(poly_coeffs * scaling_factor, ankle_angles); % m

external_moment_arm = filt_ank_mom ./ filt_Fr; % m, used during stance only

%% MTU force and strain
Fmtu = (filt_norm_ank_mom ./ internal_moment_arms) * mass; % N
lin_strain = (Fmtu / lin_k) / Lo * 100; % percent

%% Time-normalize every stance to 0-100% and average across stances
stance_pct = (0:100)';
strain_s = nan(101, num_strides);
grf_s    = nan(101, num_strides);
fmtu_s   = nan(101, num_strides);
r_s      = nan(101, num_strides);
R_s      = nan(101, num_strides);
for i = 1:num_strides
    idx = heelstrike_index(i):toeoff_index(i);
    pct = (idx - idx(1)) / (idx(end) - idx(1)) * 100;
    strain_s(:, i) = interp1(pct, lin_strain(idx), stance_pct);
    grf_s(:, i)    = interp1(pct, filt_Fr(idx) / body_weight, stance_pct);
    fmtu_s(:, i)   = interp1(pct, Fmtu(idx) / body_weight, stance_pct);
    r_s(:, i)      = interp1(pct, internal_moment_arms(idx), stance_pct);
    R_s(:, i)      = interp1(pct, external_moment_arm(idx), stance_pct);
end

curve.stance_pct  = stance_pct;
curve.strain_mean = mean(strain_s, 2, 'omitnan');
curve.strain_sd   = std(strain_s, 0, 2, 'omitnan');
curve.grf_bw      = mean(grf_s, 2, 'omitnan');
curve.fmtu_bw     = mean(fmtu_s, 2, 'omitnan');
curve.r_internal  = mean(r_s, 2, 'omitnan');
curve.R_external  = mean(R_s, 2, 'omitnan');
curve.ema_median  = median(r_s ./ R_s, 2, 'omitnan');
curve.n_stances   = num_strides;

end
