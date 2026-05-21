% GRF_Analysis
%close all
%clear
%clc

%% Importing data
% [filename, pathname] = uigetfile('Please select a file');
%fname = fullfile(pathname, filename);
% load(fname);
MECH_Data = importdata("C:\Users\Kids\Downloads\Data (1)\Data\MECH_Data.mat");


% Input Stiffness Values
lin_k = 186.7119; % N/mm
fun_k = 356.5763; % N/mm

% Input CSA Values
CSA = 55.44444; % mm^2

% Input slack lenth values
Lo = 160; % mm

% Selecting GRF data
Fy = -(MECH_Data.AA17.DayA_01.FLAT.S13.allData.Force_Fy1(1:15000)); % N
Fz = -(MECH_Data.AA17.DayA_01.FLAT.S13.allData.Force_Fz1(1:15000)); % N

% Loading GRF sampling frequency
fs = 1000;

speed = 1.3; % m/s

time = ((1:length(Fz))/fs)'; % s

%% Filtering data

fn = fs/2; % Nyquist Frequency (half of sample rate, Hz)
fc = 20; % Cutoff frequency (Hz)

% Creating 4th order Butterworth filter
[b, a] = butter(4, fc/fn);

% filtering data
filt_Fy = filtfilt(b, a, Fy);
filt_Fz = filtfilt(b, a, Fz);

%% Identifying Toe off and Heel Strike

% Setting threshold for identify stance vs swing phase
threshold = 100;

% Initializing arrays
heelstrike_index = [];
toeoff_index = [];

% Iterating through data to identify toe off and heel strike based off
% threshold passing
for i = 2:length(time)
    if filt_Fz(i-1) < threshold && filt_Fz(i) > threshold
        heelstrike_index = [heelstrike_index, i];

    elseif filt_Fz(i-1) > threshold && filt_Fz(i) < threshold
        toeoff_index = [toeoff_index, i];
    end
end

% Eliminating partial stance phases
%if toeoff_index(1) < heelstrike_index(1)
  %  toeoff_index = toeoff_index(2:end);
%end

%if toeoff_index(end) < heelstrike_index(end)
%    heelstrike_index = heelstrike_index(1:end - 1);
%end

% Calculating number of steps per leg
num_steps = length(heelstrike_index);

% Plotting toe offs and heel strikes
 

sgtitle('Plotting Heel Strikes and Toe Offs')

%% Calculating Bodyweight

% calculating bodyweight for each step
body_weight_step = [];
for i = 1:(num_steps - 1)
    body_weight_step(i) = mean(filt_Fz(heelstrike_index(i):heelstrike_index(i + 1)))*2; % N
end

body_weight = mean(body_weight_step);

% converting to kg
mass = body_weight/9.81; % kg

% Normalizing GRFs by body weight
filt_Fy_norm = filt_Fy/body_weight; % Nm/kg
filt_Fz_norm = filt_Fz/body_weight; % Nm/kg

% Calculating Resultant GRF
filt_Fr = sqrt((filt_Fz.^2)+(filt_Fy.^2));

%% Calculating Output Variables

% Calculating contact time, stride time, duty factor, stride length
contact_time = [];
stride_time = [];
stride_freq = [];
stride_df = [];
stride_length = [];
for i = 1:(num_steps)
    contact_time(i) = time(toeoff_index(i)) - time(heelstrike_index(i)); % s
end

for i = 1:(num_steps - 1)
    stride_time(i) = time(heelstrike_index(i + 1)) - time(heelstrike_index(i)); % s

    stride_freq(i) = 1/stride_time(i); % Hz

    % duty factor = ratio of contact time to stride time??
    stride_df(i) = contact_time(i) / stride_time(i);

    stride_length(i) = stride_time(i)*speed; % m
end

% Finding impact and active peaks
[peak_Fz, peak_Fz_index] = findpeaks(filt_Fz_norm(heelstrike_index(1):toeoff_index(end)), ...
    "MinPeakDistance", mean(contact_time)*fs/3, ...
    "MinPeakHeight", max(filt_Fz_norm)/2);
impact_peak_Fz = peak_Fz(1:2:end);
impact_peak_Fz_index = peak_Fz_index(1:2:end) + heelstrike_index(1);
active_peak_Fz = peak_Fz(2:2:end);
active_peak_Fz_index = peak_Fz_index(2:2:end) + heelstrike_index(1);

% Plotting impact and active peaks
plot(time, filt_Fz_norm)
hold on
plot(time(impact_peak_Fz_index), impact_peak_Fz, 'rx')
plot(time(active_peak_Fz_index), active_peak_Fz, 'ro')
hold off
xlabel('Time (s)')
ylabel('Normalized Vertical GRF')
legend('F_{z}', 'Impact Peak', 'Active Peak')
title('Right')

% Finding Peak Braking and Propulsive forces
[propel_peak, propel_peak_index] = findpeaks(filt_Fy_norm(heelstrike_index(1):toeoff_index(end)), ...
    "MinPeakDistance", mean(contact_time)*fs/2, ...
    "MinPeakHeight", max(filt_Fy_norm)/1.5);

[brake_peak, brake_peak_index] = findpeaks(-filt_Fy_norm(heelstrike_index(1):toeoff_index(end)), ...
    "MinPeakDistance", mean(contact_time)*fs/2, ...
    "MinPeakHeight", max(-filt_Fy_norm)/1.5);

brake_peak = -brake_peak;

propel_peak_index = propel_peak_index + heelstrike_index(1);
brake_peak_index = brake_peak_index + heelstrike_index(1);

% Plotting propulsive and braking peaks
plot(time, filt_Fy_norm)
hold on
plot(time(propel_peak_index), propel_peak, 'rx')
plot(time(brake_peak_index), brake_peak, 'ro')
hold off
xlabel('Time (s)')
ylabel('Normalized Horizontal GRF')
legend('F_{y}', 'Propel Peak', 'Brake Peak')
title('Right')

% Calculating Impulses
stance_index = [];
vertical_impulse = [];
brake_index = [];
brake_impulse = [];
propel_index = [];
propel_impulse = [];
for i = 1:(num_steps)
    stance_index = heelstrike_index(i):toeoff_index(i);

    vertical_impulse(i) = trapz(time(stance_index), filt_Fz_norm(stance_index)); %Ns/kg

    brake_index = find(filt_Fy_norm(stance_index) < 0);

    brake_impulse(i) = trapz(time(stance_index(brake_index)), filt_Fy_norm(stance_index(brake_index))); %Ns/kg

    propel_index = find(filt_Fy_norm(stance_index) > 0);

    propel_impulse(i) = trapz(time(stance_index(propel_index)), filt_Fy_norm(stance_index(propel_index))); %Ns/kg
end

% Calculating Mean Peak GRF
peak_Fz = [];
peak_Fz_index = [];
for i = 1:num_steps
    idx1 = heelstrike_index(i);
    idx2 = toeoff_index(i);

    segment = filt_Fz(idx1:idx2);

    % Find all peaks in the segment
    [pks, locs] = findpeaks(segment, "MinPeakHeight", max(segment)/1.5);

    if isempty(pks)
        % No peak found
        peak_Fz(i) = NaN;
        peak_Fz_index(i) = NaN;
    else
        % Take the largest peak
        [peak_Fz(i), m] = max(pks);
        peak_Fz_index(i) = locs(m) + idx1 - 1;
    end
end
mean_peak_Fz = mean(peak_Fz);


peak_Fy = [];
peak_Fy_index = [];
for i = 1:num_steps
    idx1 = heelstrike_index(i);
    idx2 = toeoff_index(i);

    segment = filt_Fy(idx1:idx2);

    % Find all peaks in the segment
    [pks, locs] = findpeaks(segment, "MinPeakHeight", max(segment)/1.5);

    if isempty(pks)
        % No peak found
        peak_Fy(i) = NaN;
        peak_Fy_index(i) = NaN;
    else
        % Take the largest peak
        [peak_Fy(i), m] = max(pks);
        peak_Fy_index(i) = locs(m) + idx1 - 1;
    end
end
mean_peak_Fy = mean(peak_Fy);


peak_Fr = [];
peak_Fr_index = [];
for i = 1:num_steps
    idx1 = heelstrike_index(i);
    idx2 = toeoff_index(i);

    segment = filt_Fr(idx1:idx2);

    % Find all peaks in the segment
    [pks, locs] = findpeaks(segment, "MinPeakHeight", max(segment)/1.5);

    if isempty(pks)
        % No peak found
        peak_Fr(i) = NaN;
        peak_Fr_index(i) = NaN;
    else
        % Take the largest peak
        [peak_Fr(i), m] = max(pks);
        peak_Fr_index(i) = locs(m) + idx1 - 1;
    end
end
mean_peak_Fr = mean(peak_Fr);

%% Moment Calculations
% Importing data

% [filename, pathname] = uigetfile('Please select a file');
%fname = fullfile(pathname, filename);
% load(fname);
JointMomentData = importdata("C:\Users\Kids\Downloads\Data (1)\Data\JointMomentData.mat");

% Selecting moment data
norm_ank_mom = (JointMomentData.AA17.DayA_01.FLAT.AnkMom'); % Nm/kg


% filtering data
fc_mom = 6; % Cutoff frequency (Hz)

% Creating 4th order Butterworth filter
[b, a] = butter(4, fc_mom/fn);
filt_norm_ank_mom = filtfilt(b, a, norm_ank_mom);

% Plotting toe offs and heel strikes
figure

plot(time, filt_norm_ank_mom, heelstrike_index./fs, filt_norm_ank_mom(heelstrike_index), 'rx', toeoff_index./fs, filt_norm_ank_mom(toeoff_index), 'ro')
xlabel('Time (s)')
ylabel('Moment (Nm/kg)')
legend('Ankle Moment', 'Heel Strike', 'Toe off')
title('Right')

% Finding peak moments
peak_norm_ank_mom = [];
peak_norm_ank_mom_index = [];
for i = 1:num_steps
    idx1 = heelstrike_index(i);
    idx2 = toeoff_index(i);

    segment = filt_norm_ank_mom(idx1:idx2);

    % Find all peaks in the segment
    [pks, locs] = findpeaks(segment, "MinPeakHeight", max(segment)/1.5);

    if isempty(pks)
        % No peak found
        peak_norm_ank_mom(i) = NaN;
        peak_norm_ank_mom_index(i) = NaN;
    else
        % Take the largest peak
        [peak_norm_ank_mom(i), m] = max(pks);
        peak_norm_ank_mom_index(i) = locs(m) + idx1 - 1;
    end
end

% Plotting peak moments
plot(time, filt_norm_ank_mom)
hold on
plot(time(peak_norm_ank_mom_index), peak_norm_ank_mom, 'rx')
hold off
xlabel('Time (s)')
ylabel('Ankle Moment (Nm)')
legend('Ankle Moment', 'Peak')
title('Right')

% Calculating mean moment
step_mean_norm_ank_mom = [];
for i = 1:(num_steps)
 step_mean_norm_ank_mom(i) = mean(filt_norm_ank_mom(heelstrike_index(i):toeoff_index(i))); % Nm/kg
end
mean_norm_ank_mom = mean(step_mean_norm_ank_mom);

% Calculating moment impulse
moment_impulse = [];
for i = 1:(num_steps)
    moment_impulse(i) = trapz(time(heelstrike_index(i):toeoff_index(i)), filt_norm_ank_mom(heelstrike_index(i):toeoff_index(i))); % Ns/kg
end


%% GRF Moment Arm
GRF_MA = zeros(15000,1); % m

% Calculate through whole stance phases
for i = 1:(length(time))
    for j = 1:(num_steps)
        if i >= heelstrike_index(j) && i <= toeoff_index(j) 
            GRF_MA(i) = filt_norm_ank_mom(i) * mass ./ filt_Fr(i); % m 
        end
    end
end
mean_GRF_MA = mean(GRF_MA); % m



%% Calculate Dynamic Moment Arms
% Importing data

% [filename, pathname] = uigetfile('Please select a file');
% fname = fullfile(pathname, filename);
% load(fname);
JointAngleData = importdata("C:\Users\Kids\Downloads\Data (1)\Data\JointAngleData.mat");

% Selecting resting MA data
rest_AtMA = 100 * (JointAngleData.AA17.rest_AtMA'); % cm

% Selecting ankle angle data
ankle_angles = (JointAngleData.AA17.DayA_01.FLAT.AnkAng');

% Maganaris et al. (2000) data
maganaris_ankle_angles = [-15, 0, 15, 30]; % degrees
maganaris_moment_arms = [4.3, 4.7, 5.2, 5.6]; % cm, resting values

poly_coeffs = polyfit(maganaris_ankle_angles, maganaris_moment_arms, 2);
moment_arm_neutral_generic = polyval(poly_coeffs, 0);
subject_moment_arm_neutral = rest_AtMA; % cm, subject-specific
scaling_factor = subject_moment_arm_neutral / moment_arm_neutral_generic;
subject_poly_coeffs = poly_coeffs * scaling_factor;

% Calculate moment arms
moment_arms = [];
for i = 1:(length(time))
    moment_arms(i) = polyval(subject_poly_coeffs, ankle_angles(i)) / 100; % m
end

%% Calculating Fmtu
Fmtu = [];
for i = 1:(length(time))
    Fmtu(i) = filt_norm_ank_mom(i) * mass / moment_arms(i); % N
end

% Average Fmtu
step_mean_Fmtu = [];
for i = 1:(num_steps)
    step_mean_Fmtu(i) = mean(Fmtu(heelstrike_index(i):toeoff_index(i))); % N
end
mean_Fmtu = mean(step_mean_Fmtu);

% Peak Fmtu
peak_Fmtu_index = [];
peak_Fmtu = [];
for i = 1:num_steps
    idx1 = heelstrike_index(i);
    idx2 = toeoff_index(i);

    segment = Fmtu(idx1:idx2);

    % Find all peaks in the segment
    [pks, locs] = findpeaks(segment, "MinPeakHeight", max(segment)/1.5);

    if isempty(pks)
        % No peak found
        peak_Fmtu(i) = NaN;
        peak_Fmtu_index(i) = NaN;
    else
        % Take the largest peak
        [peak_Fmtu(i), m] = max(pks);
        peak_Fmtu_index(i) = locs(m) + idx1 - 1;
    end
end

% Plotting peak Fmtu
plot(time, Fmtu)
hold on
plot(time(peak_Fmtu_index), peak_Fmtu, 'rx')
hold off
xlabel('Time (s)')
ylabel('Plantar Flexor Fmtu(N)')
legend('Fmtu (N/kg)', 'Peak')
title('Right')

%% Calculating Strain

% Calculating Displacement using Linear Stiffness
lin_displacement_magnitude = [];
for i = 1:(length(time))
    lin_displacement_magnitude(i) = abs(Fmtu(i) / lin_k); % mm; 
end

% Calculating Displacement using Functional Stiffness
fun_displacement_magnitude = [];
for i = 1:(length(time))
    fun_displacement_magnitude(i) = abs(Fmtu(i) / fun_k); % mm; 
end

% Calculating strain (Normalizing displacement)
lin_strain = [];
for i = 1:(length(time))
    lin_strain(i) = lin_displacement_magnitude(i) / Lo * 100; % percent
end

fun_strain = [];
for i = 1:(length(time))
    fun_strain(i) = fun_displacement_magnitude(i) / Lo * 100; % percent
end

% Finding Peak Strain (Linear stiffness)
peak_lin_strain_index = [];
peak_lin_strain = [];
for i = 1:num_steps
    idx1 = heelstrike_index(i);
    idx2 = toeoff_index(i);

    segment = lin_strain(idx1:idx2);

    % Find all peaks in the segment
    [pks, locs] = findpeaks(segment, "MinPeakHeight", max(segment)/1.5);

    if isempty(pks)
        % No peak found
        peak_lin_strain(i) = NaN;
        peak_lin_strain_index(i) = NaN;
    else
        % Take the largest peak
        [peak_lin_strain(i), m] = max(pks);
        peak_lin_strain_index(i) = locs(m) + idx1 - 1;
    end
end

% Plotting peak strain
plot(time, lin_strain)
hold on
plot(time(peak_lin_strain_index), peak_lin_strain, 'rx')
hold off
xlabel('Time (s)')
ylabel('Strain (%)')
legend('Strain', 'Peak')
title('Right')

avg_peak_lin_strain = mean(peak_lin_strain);

% Finding Peak Strain (Functional stiffness)
peak_fun_strain = [];
peak_fun_strain_index = [];
for i = 1:num_steps
    idx1 = heelstrike_index(i);
    idx2 = toeoff_index(i);

    segment = fun_strain(idx1:idx2);

    % Find all peaks in the segment
    [pks, locs] = findpeaks(segment, "MinPeakHeight", max(segment)/1.5);

    if isempty(pks)
        % No peak found
        peak_fun_strain(i) = NaN;
        peak_fun_strain_index(i) = NaN;
    else
        % Take the largest peak
        [peak_fun_strain(i), m] = max(pks);
        peak_fun_strain_index(i) = locs(m) + idx1 - 1;
    end
end

% Plotting peak strain
plot(time, fun_strain)
hold on
plot(time(peak_fun_strain_index), peak_fun_strain, 'rx')
hold off
xlabel('Time (s)')
ylabel('Strain (%)')
legend('Strain', 'Peak')
title('Right')

avg_peak_fun_strain = mean(peak_fun_strain);

% Finding mean strain
step_mean_lin_strain = [];
for i = 1:(num_steps)
 step_mean_lin_strain(i) = mean(lin_strain(heelstrike_index(i):toeoff_index(i)));
end
mean_lin_strain = mean(step_mean_lin_strain);

step_mean_fun_strain = [];
for i = 1:(num_steps)
 step_mean_fun_strain(i) = mean(fun_strain(heelstrike_index(i):toeoff_index(i)));
end
mean_fun_strain = mean(step_mean_fun_strain);

% Calculating strain impulse
lin_strain_impulse = [];
for i = 1:(num_steps)
    lin_strain_impulse(i) = trapz(time(heelstrike_index(i):toeoff_index(i)), lin_strain(heelstrike_index(i):toeoff_index(i)));
end
mean_lin_strain_impulse = mean(lin_strain_impulse);

fun_strain_impulse = [];
for i = 1:(num_steps)
    fun_strain_impulse(i) = trapz(time(heelstrike_index(i):toeoff_index(i)), fun_strain(heelstrike_index(i):toeoff_index(i)));
end
mean_fun_strain_impulse = mean(fun_strain_impulse);

%% EMA
EMA = [];
for i=1:length(time)
   EMA(i) = moment_arms(i) / GRF_MA(i);
end

peak_mom_EMA = [];
for i = 1:(num_steps)
    peak_mom_EMA(i) = EMA(peak_norm_ank_mom_index(i));
end
avg_peak_EMA = mean(peak_mom_EMA);

%% Young's Modulus
stress = [];
for i = 1:length(time)
    stress(i) = Fmtu(i) / CSA; % N/mm^2
end

for i = 1:(num_steps)
    peak_mom_E(i) = stress(peak_norm_ank_mom_index(i)) / lin_strain(peak_norm_ank_mom_index(i));
end
E = mean(peak_mom_E); % N/mm^2