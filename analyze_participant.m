function results = analyze_participant(pid, lin_k, rest_AtMA, CSA, Lo, ...
                                         MECH_Data, JointMomentData, JointAngleData)
% ANALYZE_PARTICIPANT  Run the GRF / Achilles tendon pipeline for one participant.
%
%   results = analyze_participant(pid, lin_k, rest_AtMA, CSA, Lo, ...
%                                  MECH_Data, JointMomentData, JointAngleData, speed)
%
% INPUTS
%   pid     - participant field name as it appears in the data structs, e.g. "SS11"
%   lin_k   - linear tendon stiffness (N/mm)
%   rest_AtMA - resting Achilles tendon moment arm (m)
%   CSA     - tendon CSA (mm^2)
%   Lo      - tendon slack length (mm)
%   MECH_Data, JointMomentData, JointAngleData - the three loaded .mat structs
%             (each one contains ALL participants as top-level fields)
%
% OUTPUT
%   results - struct with fields:
%       ParticipantID, E, mean_peak_Fmtu_EMA, mean_lin_strain_impulse,
%       mean_peak_lin_strain, mean_lin_strain, mean_peak_Fr,
%       mass, num_strides (diagnostics)

speed = 1.3; % m/s

pid = char(pid); % allow string or char input for dynamic field access

% --- Trial labels: identical for every participant ---
day   = 'DayA_03';
cond  = 'FLAT';
trial = 'S13';

fs = 1000;        % Hz, GRF sampling frequency
fn = fs/2;        % Nyquist frequency
fc = 10;          % GRF filter cutoff (Hz)
fc_mom = 10;      % moment filter cutoff (Hz)
threshold = 10;   % N, stance/swing threshold

%% Ground reaction forces
Fy_full = -(MECH_Data.(pid).(day).(cond).(trial).allData.Force_Fy1); % N
Fz_full = -(MECH_Data.(pid).(day).(cond).(trial).allData.Force_Fz1); % N

% Some GRF recordings were sampled at a multiple of 1000 Hz (SE15 DayA_03:
% 2000 Hz, 30000 samples for the same 15 s as the 15000-sample moment data).
% Filter at the native rate, then keep every grf_ratio-th sample so the GRF
% lines up sample-for-sample with the 1000 Hz moment and angle data.
n_mom = length(JointMomentData.(pid).(day).(cond).AnkMom);
grf_ratio = round(length(Fz_full) / n_mom);
n = min(15000, floor(length(Fz_full) / grf_ratio));
time = ((1:n)/fs)'; % s

[b, a] = butter(4, fc / (fs * grf_ratio / 2));
filt_Fy = filtfilt(b, a, Fy_full(1:n * grf_ratio));
filt_Fz = filtfilt(b, a, Fz_full(1:n * grf_ratio));
filt_Fy = filt_Fy(1:grf_ratio:end);
filt_Fz = filt_Fz(1:grf_ratio:end);

%% Heel strike / toe off detection
heelstrike_index = [];
toeoff_index = [];
for i = 2:length(time)
    if filt_Fz(i-1) < threshold && filt_Fz(i) > threshold
        heelstrike_index(end+1) = i; 
    elseif filt_Fz(i-1) > threshold && filt_Fz(i) < threshold
        toeoff_index(end+1) = i; 
    end
end

if isempty(heelstrike_index) || isempty(toeoff_index)
    error('analyze_participant:noSteps', 'No stance phases detected for %s.', pid);
end

% Remove partial stance phases
if toeoff_index(1) < heelstrike_index(1)
    toeoff_index = toeoff_index(2:end);
end
if toeoff_index(end) < heelstrike_index(end)
    heelstrike_index = heelstrike_index(1:end-1);
end

num_strides = min(length(heelstrike_index), length(toeoff_index));
heelstrike_index = heelstrike_index(1:num_strides);
toeoff_index = toeoff_index(1:num_strides);

%% Body weight / mass
body_weight_stride = zeros(1, num_strides-1);
for i = 1:(num_strides-1)
    body_weight_stride(i) = mean(filt_Fz(heelstrike_index(i):heelstrike_index(i+1)))*2; % N
end
body_weight = mean(body_weight_stride); % N
mass = body_weight/9.81; % kg

filt_Fr = sqrt(filt_Fz.^2 + filt_Fy.^2); % N, resultant GRF, not normalized

%% Mean peak resultant GRF
peak_Fr = nan(1, num_strides);
for i = 1:num_strides
    idx1 = heelstrike_index(i);
    idx2 = toeoff_index(i);
    segment = filt_Fr(idx1:idx2);
    pks = findpeaks(segment, "MinPeakHeight", max(segment)/1.5);
    if ~isempty(pks)
        peak_Fr(i) = max(pks);
    end
end
mean_peak_Fr = mean(peak_Fr, 'omitnan');

%% Ankle moment
norm_ank_mom = (JointMomentData.(pid).(day).(cond).AnkMom'); % Nm/kg
ank_mom = norm_ank_mom * mass;

[b, a] = butter(4, fc_mom/fn);
filt_norm_ank_mom = filtfilt(b, a, norm_ank_mom);
filt_ank_mom = filtfilt(b, a, ank_mom);

peak_ank_mom = nan(1, num_strides);
peak_ank_mom_index= nan(1, num_strides);
for i = 1:num_strides
    idx1 = heelstrike_index(i);
    idx2 = toeoff_index(i);
    segment = filt_norm_ank_mom(idx1:idx2);
    [pks, locs] = findpeaks(segment, "MinPeakHeight", max(segment)/1.5);
    if ~isempty(pks)
        [peak_ank_mom(i), m] = max(pks);
        peak_ank_mom_index(i) = locs(m) + idx1 - 1;
    end
end
mean_peak_ank_mom = mean(peak_ank_mom);
mean_ank_mom = mean(nonzeros(norm_ank_mom));

%% External (GRF) moment arm, computed during stance only
external_moment_arm = zeros(n, 1);
for j = 1:num_strides
    idx1 = heelstrike_index(j);
    idx2 = toeoff_index(j);
    external_moment_arm(idx1:idx2) = filt_ank_mom(idx1:idx2) ./ filt_Fr(idx1:idx2); % m
end

%% Internal (Achilles tendon) moment arm, scaled to subject
ankle_angles = (JointAngleData.(pid).(day).(cond).AnkAng');

maganaris_ankle_angles = [-15, 0, 15, 30];     % deg
maganaris_moment_arms  = [4.3, 4.7, 5.2, 5.6]; % cm, resting values

poly_coeffs = polyfit(maganaris_ankle_angles, maganaris_moment_arms, 2);
moment_arm_neutral_generic = polyval(poly_coeffs, 0);
scaling_factor = rest_AtMA / moment_arm_neutral_generic;
subject_poly_coeffs = poly_coeffs * scaling_factor;

internal_moment_arms = polyval(subject_poly_coeffs, ankle_angles); % m

%% Muscle-tendon unit force (Fmtu)
Fmtu = (filt_norm_ank_mom ./ internal_moment_arms) * mass; % N

peak_Fmtu = nan(1, num_strides);
peak_Fmtu_index = nan(1, num_strides);
for i = 1:num_strides
    idx1 = heelstrike_index(i);
    idx2 = toeoff_index(i);
    segment = Fmtu(idx1:idx2);
    [pks, locs] = findpeaks(segment, "MinPeakHeight", max(segment)/1.5);
    if ~isempty(pks)
        [peak_Fmtu(i), m] = max(pks);
        peak_Fmtu_index(i) = locs(m) + idx1 - 1;
    end
end
mean_peak_Fmtu = mean(peak_Fmtu);
mean_Fmtu = mean(nonzeros(Fmtu));


%% Strain (linear stiffness)
lin_displacement = Fmtu / lin_k; % mm
lin_strain = lin_displacement / Lo * 100; % percent

peak_lin_strain = nan(1, num_strides);
for i = 1:num_strides
    idx1 = heelstrike_index(i);
    idx2 = toeoff_index(i);
    segment = lin_strain(idx1:idx2);
    pks = findpeaks(segment, "MinPeakHeight", max(segment)/1.5);
    if ~isempty(pks)
        peak_lin_strain(i) = max(pks);
    end
end
mean_peak_lin_strain = mean(peak_lin_strain, 'omitnan');
mean_lin_strain = mean(nonzeros(lin_strain));

lin_strain_impulse = nan(1, num_strides);
for i = 1:num_strides
    idx1 = heelstrike_index(i);
    idx2 = toeoff_index(i);
    lin_strain_impulse(i) = trapz(time(idx1:idx2), lin_strain(idx1:idx2));
end
mean_lin_strain_impulse = mean(lin_strain_impulse, 'omitnan');

%% Effective mechanical advantage (EMA)
EMA = internal_moment_arms ./ external_moment_arm;

peak_Fmtu_EMA = nan(1, num_strides);
for i = 1:num_strides
    if ~isnan(peak_Fmtu_index(i))
        peak_Fmtu_EMA(i) = EMA(peak_Fmtu_index(i));
    end
end
mean_peak_Fmtu_EMA = mean(peak_Fmtu_EMA, 'omitnan');

%% Young's modulus at peak Fmtu
stress = Fmtu / CSA; % N/mm^2

peak_Fmtu_E = nan(1, num_strides);
for i = 1:num_strides
    if ~isnan(peak_Fmtu_index(i))
        idx = peak_Fmtu_index(i);
        peak_Fmtu_E(i) = stress(idx) / lin_strain(idx);
    end
end
E = mean(peak_Fmtu_E, 'omitnan');

%% Package results
results.ParticipantID = pid;
results.E = E;
results.mean_peak_Fmtu_EMA = mean_peak_Fmtu_EMA;
results.mean_lin_strain_impulse = mean_lin_strain_impulse;
results.mean_peak_lin_strain = mean_peak_lin_strain;
results.mean_lin_strain = mean_lin_strain;
results.mean_peak_Fr = mean_peak_Fr;
results.mass = mass;
results.num_strides = num_strides;
results.mean_peak_Fmtu = mean_peak_Fmtu;
results.mean_Fmtu = mean_Fmtu;
results.mean_peak_ank_mom = mean_peak_ank_mom;
results.mean_ank_mom = mean_ank_mom;


end
