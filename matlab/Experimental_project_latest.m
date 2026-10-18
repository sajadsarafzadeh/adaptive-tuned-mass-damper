%% Experimental project

clear all
close all
clc

fsamp = 2000;
lb = 10;                            % Lower bound in Hz
ub = 500;                           % Upper bound in Hz
Tseg = 10;                          % Window length
olap = 0.6;                         % Overlap 60 %
Nseg = Tseg * fsamp;                % Samples per window
Win = hann(Nseg);                   % Hanning window function
step = round(Nseg * (1 - olap));    % Shift between windows
sensor_sens = 0.0102;               % Sensitivity in V/[m/s^2]
actuator_sens = 0.0224;             % Sensitivity in V/N

m_clamp = 0.16;                     % Mass of the clamping system [kg]

%% Current 0 A
data0 = load('Selected_data/Current0/merged_data0.mat');
acc0   = data0.Dati(:,1)/sensor_sens;       % Acceleration [m/s^2]
force0 = data0.Dati(:,2)/actuator_sens;     % Force [N]

nSamples = length(acc0);
nWindows = floor((nSamples - Nseg)/step) + 1;

% Initialize using the correct input signal (acceleration) for frequency vector
[~, f] = autocross(acc0(1:Nseg), acc0(1:Nseg), fsamp, Nseg, Win);
f = f(:);
nFreq = length(f);

G_AA_all = zeros(nFreq, nWindows);
G_AF_all = zeros(nFreq, nWindows);
G_FF_all = zeros(nFreq, nWindows);

for iWin = 1:nWindows
    idx_start = (iWin-1)*step + 1;
    idx_end   = idx_start + Nseg - 1;
    acc_win   = acc0(idx_start:idx_end);
    force_win = force0(idx_start:idx_end);
    
    % Correct H1 Spectral Densities: Input = Acc, Output = Force
    G_AA_all(:, iWin) = autocross(acc_win, acc_win, fsamp, Nseg, Win);
    G_AF_all(:, iWin) = autocross(acc_win, force_win, fsamp, Nseg, Win);
    G_FF_all(:, iWin) = autocross(force_win, force_win, fsamp, Nseg, Win);
end
%/MATLAB Drive/Experimental dynamics_1/Selected_data/Current6.0/test@3861706643.mat
% Average spectra
G_AA_0_avg = mean(G_AA_all, 2);
G_AF_0_avg = mean(G_AF_all, 2);
G_FF_0_avg = mean(G_FF_all, 2);

gamma2_0_avg = (abs(G_AF_0_avg).^2) ./ (G_AA_0_avg .* G_FF_0_avg);

% Compute measured H1 Apparent Mass [N / (m/s^2)]
H_app_0_avg = G_AF_0_avg ./ G_AA_0_avg;

% Select active frequency bounds
idx_freq    = (f >= lb) & (f <= ub);
f_plot      = f(idx_freq);
H_app_0_plot = H_app_0_avg(idx_freq);
gamma2_0_plot = gamma2_0_avg(idx_freq);

% Apply physical clamp mass cancellation
H_app_0_true = H_app_0_plot + m_clamp;

% Convert Apparent Mass to Dynamic Stiffness (N/m) 
omega = 2 * pi * f_plot;
G_SY_0_plot = H_app_0_true .* (-omega.^2);

%% Current 0.5 A
data05 = load('Selected_data/Current0.5/test@3861698019.mat');
acc05   = data05.Dati(:,1)/sensor_sens;     % Acceleration [m/s^2]
force05 = data05.Dati(:,2)/actuator_sens;   % Force [N]
nSamples = length(acc05);
nWindows = floor((nSamples - Nseg)/step) + 1;
G_AA_all = zeros(nFreq, nWindows);
G_AF_all = zeros(nFreq, nWindows);
G_FF_all = zeros(nFreq, nWindows);
for iWin = 1:nWindows
    idx_start = (iWin-1)*step + 1;
    idx_end   = idx_start + Nseg - 1;
    acc_win   = acc05(idx_start:idx_end);
    force_win = force05(idx_start:idx_end);
    
    % Correct H1 Spectral Densities: Input = Acc, Output = Force
    G_AA_all(:, iWin) = autocross(acc_win, acc_win, fsamp, Nseg, Win);
    G_AF_all(:, iWin) = autocross(acc_win, force_win, fsamp, Nseg, Win);
    G_FF_all(:, iWin) = autocross(force_win, force_win, fsamp, Nseg, Win);
end
% Average spectra
G_AA_05_avg = mean(G_AA_all, 2);
G_AF_05_avg = mean(G_AF_all, 2);
G_FF_05_avg = mean(G_FF_all, 2);
gamma2_05_avg = (abs(G_AF_05_avg).^2) ./ (G_AA_05_avg .* G_FF_05_avg);

% Compute measured H1 Apparent Mass [N / (m/s^2)]
H_app_05_avg = G_AF_05_avg ./ G_AA_05_avg;

% Select active frequency bounds
H_app_05_plot = H_app_05_avg(idx_freq);
gamma2_05_plot = gamma2_05_avg(idx_freq);

% Apply physical clamp mass cancellation
H_app_05_true = H_app_05_plot + m_clamp;

% Convert Apparent Mass to Dynamic Stiffness (N/m) 
G_SY_05_plot = H_app_05_true .* (-omega.^2);

%% Current 1 A

data1 = load('Selected_data/Current1.0/merged_data1.mat');
acc1   = data1.Dati(:,1)/sensor_sens;       % Acceleration [m/s^2]
force1 = data1.Dati(:,2)/actuator_sens;     % Force [N]
nSamples = length(acc1);
nWindows = floor((nSamples - Nseg)/step) + 1;
% Initialize using the acceleration signal for frequency vector
[~, f] = autocross(acc1(1:Nseg), acc1(1:Nseg), fsamp, Nseg, Win);
f = f(:);
nFreq = length(f);
G_AA_all = zeros(nFreq, nWindows);
G_AF_all = zeros(nFreq, nWindows);
G_FF_all = zeros(nFreq, nWindows);
for iWin = 1:nWindows
    idx_start = (iWin-1)*step + 1;
    idx_end   = idx_start + Nseg - 1;
    acc_win   = acc1(idx_start:idx_end);
    force_win = force1(idx_start:idx_end);
    
    % Correct H1 Spectral Densities: Input = Acc, Output = Force
    G_AA_all(:, iWin) = autocross(acc_win, acc_win, fsamp, Nseg, Win);
    G_AF_all(:, iWin) = autocross(acc_win, force_win, fsamp, Nseg, Win);
    G_FF_all(:, iWin) = autocross(force_win, force_win, fsamp, Nseg, Win);
end
% Average spectra
G_AA_avg = mean(G_AA_all, 2);
G_AF_avg = mean(G_AF_all, 2);
G_FF_avg = mean(G_FF_all, 2);
gamma2_1_avg = (abs(G_AF_avg).^2) ./ (G_AA_avg .* G_FF_avg);
% Compute measured H1 Apparent Mass [N / (m/s^2)]
H_app_1_avg = G_AF_avg ./ G_AA_avg;

% Select active frequency bounds
idx_freq    = (f >= lb) & (f <= ub);
f_plot      = f(idx_freq);
H_app_1_plot = H_app_1_avg(idx_freq);
gamma2_1_plot = gamma2_1_avg(idx_freq);

% Apply physical clamp mass cancellation
H_app_1_true = H_app_1_plot + m_clamp;

% Convert Apparent Mass to Dynamic Stiffness (N/m)
omega = 2 * pi * f_plot; 
G_SY_1_plot = H_app_1_true .* (-omega.^2);

% Plot final averaged FRF

% figure
% semilogy(f_plot, abs(H1_1_plot),'r', 'LineWidth', 1.5)
% grid on
% xlabel('Frequency [Hz]')
% ylabel('|H_1(f)| [m/s^2/N]')
% title('Averaged FRF - H1 Estimator for 1 A')
% xlim([lb ub])

%% Current 2 A

data2 = load('Selected_data/Current2.0/test@3861701622.mat');
acc2   = data2.Dati(:,1)/sensor_sens;       % Acceleration [m/s^2]
force2 = data2.Dati(:,2)/actuator_sens;     % Force [N]

nSamples = length(acc2);
nWindows = floor((nSamples - Nseg)/step) + 1;

% Initialize using the acceleration signal for frequency vector
[~, f] = autocross(acc2(1:Nseg), acc2(1:Nseg), fsamp, Nseg, Win);
f = f(:);
nFreq = length(f);

G_AA_all = zeros(nFreq, nWindows);
G_AF_all = zeros(nFreq, nWindows);
G_FF_all = zeros(nFreq, nWindows);

for iWin = 1:nWindows
    idx_start = (iWin-1)*step + 1;
    idx_end   = idx_start + Nseg - 1;
    acc_win   = acc2(idx_start:idx_end);
    force_win = force2(idx_start:idx_end);
    
    % Correct H1 Spectral Densities: Input = Acc, Output = Force
    G_AA_all(:, iWin) = autocross(acc_win, acc_win, fsamp, Nseg, Win);
    G_AF_all(:, iWin) = autocross(acc_win, force_win, fsamp, Nseg, Win);
    G_FF_all(:, iWin) = autocross(force_win, force_win, fsamp, Nseg, Win);
end

% Average spectra
G_AA_avg = mean(G_AA_all, 2);
G_AF_avg = mean(G_AF_all, 2);
G_FF_avg = mean(G_FF_all, 2);
gamma2_2_avg = (abs(G_AF_avg).^2) ./ (G_AA_avg .* G_FF_avg);


% Compute measured H1 Apparent Mass [N / (m/s^2)]
H_app_2_avg = G_AF_avg ./ G_AA_avg;

% Select active frequency bounds
idx_freq    = (f >= lb) & (f <= ub);
f_plot      = f(idx_freq);
H_app_2_plot = H_app_2_avg(idx_freq);
gamma2_2_plot = gamma2_2_avg(idx_freq);

% Apply physical clamp mass cancellation
H_app_2_true = H_app_2_plot + m_clamp;

% Convert Apparent Mass to Dynamic Stiffness (N/m) 
omega = 2 * pi * f_plot; 
G_SY_2_plot = H_app_2_true .* (-omega.^2);

% Plot final averaged FRF

% figure
% semilogy(f_plot, abs(H1_2_plot), 'LineWidth', 1.5)
% grid on
% xlabel('Frequency [Hz]')
% ylabel('|H_1(f)| [m/s^2/N]')
% title('Averaged FRF - H1 Estimator for 2 A')
% xlim([lb ub])

%% Current 3 A

data3 = load('Selected_data/Current3.0/test@3861703028.mat');
acc3   = data3.Dati(:,1)/sensor_sens;       % Acceleration [m/s^2]
force3 = data3.Dati(:,2)/actuator_sens;     % Force [N]

nSamples = length(acc3);
nWindows = floor((nSamples - Nseg)/step) + 1;

% Initialize using the acceleration signal for frequency vector
[~, f] = autocross(acc3(1:Nseg), acc3(1:Nseg), fsamp, Nseg, Win);
f = f(:);
nFreq = length(f);

G_AA_all = zeros(nFreq, nWindows);
G_AF_all = zeros(nFreq, nWindows);
G_FF_all = zeros(nFreq, nWindows);

for iWin = 1:nWindows
    idx_start = (iWin-1)*step + 1;
    idx_end   = idx_start + Nseg - 1;
    acc_win   = acc3(idx_start:idx_end);
    force_win = force3(idx_start:idx_end);
    
    % Correct H1 Spectral Densities: Input = Acc, Output = Force
    G_AA_all(:, iWin) = autocross(acc_win, acc_win, fsamp, Nseg, Win);
    G_AF_all(:, iWin) = autocross(acc_win, force_win, fsamp, Nseg, Win);
    G_FF_all(:, iWin) = autocross(force_win, force_win, fsamp, Nseg, Win);
end

% Average spectra
G_AA_avg = mean(G_AA_all, 2);
G_AF_avg = mean(G_AF_all, 2);
G_FF_avg = mean(G_FF_all, 2);
gamma2_3_avg = (abs(G_AF_avg).^2) ./ (G_AA_avg .* G_FF_avg);

% Compute measured H1 Apparent Mass [N / (m/s^2)]
H_app_3_avg = G_AF_avg ./ G_AA_avg;

% Select active frequency bounds
idx_freq    = (f >= lb) & (f <= ub);
f_plot      = f(idx_freq);
H_app_3_plot = H_app_3_avg(idx_freq);
gamma2_3_plot = gamma2_3_avg(idx_freq);

% Apply physical clamp mass cancellation
H_app_3_true = H_app_3_plot + m_clamp;

% Convert Apparent Mass to Dynamic Stiffness (N/m) 
omega = 2 * pi * f_plot; 
G_SY_3_plot = H_app_3_true .* (-omega.^2);

% Plot final averaged FRF

% figure
% semilogy(f_plot, abs(H1_3_plot), 'LineWidth', 1.5)
% grid on
% xlabel('Frequency [Hz]')
% ylabel('|H_1(f)| [m/s^2/N]')
% title('Averaged FRF - H1 Estimator for 3 A')
% xlim([lb ub])

%% Current 4 A

data4 = load('Selected_data/Current4.0/test@3861704381.mat');
acc4   = data4.Dati(:,1)/sensor_sens;       % Acceleration [m/s^2]
force4 = data4.Dati(:,2)/actuator_sens;     % Force [N]

nSamples = length(acc4);
nWindows = floor((nSamples - Nseg)/step) + 1;

% Initialize using the acceleration signal for frequency vector
[~, f] = autocross(acc4(1:Nseg), acc4(1:Nseg), fsamp, Nseg, Win);
f = f(:);
nFreq = length(f);

G_AA_all = zeros(nFreq, nWindows);
G_AF_all = zeros(nFreq, nWindows);
G_FF_all = zeros(nFreq, nWindows);

for iWin = 1:nWindows
    idx_start = (iWin-1)*step + 1;
    idx_end   = idx_start + Nseg - 1;
    acc_win   = acc4(idx_start:idx_end);
    force_win = force4(idx_start:idx_end);
    
    % Correct H1 Spectral Densities: Input = Acc, Output = Force
    G_AA_all(:, iWin) = autocross(acc_win, acc_win, fsamp, Nseg, Win);
    G_AF_all(:, iWin) = autocross(acc_win, force_win, fsamp, Nseg, Win);
    G_FF_all(:, iWin) = autocross(force_win, force_win, fsamp, Nseg, Win);
end

% Average spectra
G_AA_avg = mean(G_AA_all, 2);
G_AF_avg = mean(G_AF_all, 2);
G_FF_avg = mean(G_FF_all, 2);
gamma2_4_avg = (abs(G_AF_avg).^2) ./ (G_AA_avg .* G_FF_avg);

% Compute measured H1 Apparent Mass [N / (m/s^2)]
H_app_4_avg = G_AF_avg ./ G_AA_avg;

% Select active frequency bounds
idx_freq    = (f >= lb) & (f <= ub);
f_plot      = f(idx_freq);
H_app_4_plot = H_app_4_avg(idx_freq);
gamma2_4_plot = gamma2_4_avg(idx_freq);

% Apply physical clamp mass cancellation
H_app_4_true = H_app_4_plot + m_clamp;

% Convert Apparent Mass to Dynamic Stiffness (N/m) 
omega = 2 * pi * f_plot; 
G_SY_4_plot = H_app_4_true .* (-omega.^2);

% Plot final averaged FRF

% figure
% semilogy(f_plot, abs(H1_4_plot), 'LineWidth', 1.5)
% grid on
% xlabel('Frequency [Hz]')
% ylabel('|H_1(f)| [m/s^2/N]')
% title('Averaged FRF - H1 Estimator for 4 A')
% xlim([lb ub])

%% Current 5 A

data5 = load('Selected_data/Current5.0/test@3861705458.mat');
acc5   = data5.Dati(:,1)/sensor_sens;       % Acceleration [m/s^2]
force5 = data5.Dati(:,2)/actuator_sens;     % Force [N]

nSamples = length(acc5);
nWindows = floor((nSamples - Nseg)/step) + 1;

% Initialize using the acceleration signal for frequency vector
[~, f] = autocross(acc5(1:Nseg), acc5(1:Nseg), fsamp, Nseg, Win);
f = f(:);
nFreq = length(f);

G_AA_all = zeros(nFreq, nWindows);
G_AF_all = zeros(nFreq, nWindows);
G_FF_all = zeros(nFreq, nWindows);

for iWin = 1:nWindows
    idx_start = (iWin-1)*step + 1;
    idx_end   = idx_start + Nseg - 1;
    acc_win   = acc5(idx_start:idx_end);
    force_win = force5(idx_start:idx_end);
    
    % Correct H1 Spectral Densities: Input = Acc, Output = Force
    G_AA_all(:, iWin) = autocross(acc_win, acc_win, fsamp, Nseg, Win);
    G_AF_all(:, iWin) = autocross(acc_win, force_win, fsamp, Nseg, Win);
    G_FF_all(:, iWin) = autocross(force_win, force_win, fsamp, Nseg, Win);
end

% Average spectra
G_AA_avg = mean(G_AA_all, 2);
G_AF_avg = mean(G_AF_all, 2);
G_FF_avg = mean(G_FF_all, 2);
gamma2_5_avg = (abs(G_AF_avg).^2) ./ (G_AA_avg .* G_FF_avg);

% Compute measured H1 Apparent Mass [N / (m/s^2)]
H_app_5_avg = G_AF_avg ./ G_AA_avg;

% Select active frequency bounds
idx_freq    = (f >= lb) & (f <= ub);
f_plot      = f(idx_freq);
H_app_5_plot = H_app_5_avg(idx_freq);
gamma2_5_plot = gamma2_5_avg(idx_freq);

% Apply physical clamp mass cancellation
H_app_5_true = H_app_5_plot + m_clamp;

% Convert Apparent Mass to Dynamic Stiffness (N/m) 
omega = 2 * pi * f_plot; 
G_SY_5_plot = H_app_5_true .* (-omega.^2);

% Plot final averaged FRF

% figure
% semilogy(f_plot, abs(H1_5_plot), 'LineWidth', 1.5)
% grid on
% xlabel('Frequency [Hz]')
% ylabel('|H_1(f)| [m/s^2/N]')
% title('Averaged FRF - H1 Estimator for 5 A')
% xlim([lb ub])

%% Current 6 A

data6 = load('Selected_data/Current6.0/test@3861706643.mat');
acc6   = data6.Dati(:,1)/sensor_sens;       % Acceleration [m/s^2]
force6 = data6.Dati(:,2)/actuator_sens;     % Force [N]

nSamples = length(acc6);
nWindows = floor((nSamples - Nseg)/step) + 1;

% Initialize using the acceleration signal for frequency vector
[~, f] = autocross(acc6(1:Nseg), acc6(1:Nseg), fsamp, Nseg, Win);
f = f(:);
nFreq = length(f);

G_AA_all = zeros(nFreq, nWindows);
G_AF_all = zeros(nFreq, nWindows);
G_FF_all = zeros(nFreq, nWindows);

for iWin = 1:nWindows
    idx_start = (iWin-1)*step + 1;
    idx_end   = idx_start + Nseg - 1;
    acc_win   = acc6(idx_start:idx_end);
    force_win = force6(idx_start:idx_end);
    
    % Correct H1 Spectral Densities: Input = Acc, Output = Force
    G_AA_all(:, iWin) = autocross(acc_win, acc_win, fsamp, Nseg, Win);
    G_AF_all(:, iWin) = autocross(acc_win, force_win, fsamp, Nseg, Win);
    G_FF_all(:, iWin) = autocross(force_win, force_win, fsamp, Nseg, Win);
end

% Average spectra
G_AA_avg = mean(G_AA_all, 2);
G_AF_avg = mean(G_AF_all, 2);
G_FF_avg = mean(G_FF_all, 2);
gamma2_6_avg = (abs(G_AF_avg).^2) ./ (G_AA_avg .* G_FF_avg);

% Compute measured H1 Apparent Mass [N / (m/s^2)]
H_app_6_avg = G_AF_avg ./ G_AA_avg;

% Select active frequency bounds
idx_freq    = (f >= lb) & (f <= ub);
f_plot      = f(idx_freq);
H_app_6_plot = H_app_6_avg(idx_freq);
gamma2_6_plot = gamma2_6_avg(idx_freq);

% Apply physical clamp mass cancellation
H_app_6_true = H_app_6_plot + m_clamp;

% Convert Apparent Mass to Dynamic Stiffness (N/m) 
omega = 2 * pi * f_plot; 
G_SY_6_plot = H_app_6_true .* (-omega.^2);

% Plot final averaged FRF

% figure
% semilogy(f_plot, abs(H1_6_plot), 'LineWidth', 1.5)
% grid on
% xlabel('Frequency [Hz]')
% ylabel('|H_1(f)| [m/s^2/N]')
% title('Averaged FRF - H1 Estimator for 6 A')
% xlim([lb ub])

%% Current 6.5 A

data65 = load('Selected_data/Current6.5/test@3861707544.mat');
acc65   = data65.Dati(:,1)/sensor_sens;     % Acceleration [m/s^2]
force65 = data65.Dati(:,2)/actuator_sens;   % Force [N]

nSamples = length(acc65);
nWindows = floor((nSamples - Nseg)/step) + 1;

% Initialize using the acceleration signal for frequency vector
[~, f] = autocross(acc65(1:Nseg), acc65(1:Nseg), fsamp, Nseg, Win);
f = f(:);
nFreq = length(f);

G_AA_all = zeros(nFreq, nWindows);
G_AF_all = zeros(nFreq, nWindows);
G_FF_all = zeros(nFreq, nWindows);

for iWin = 1:nWindows
    idx_start = (iWin-1)*step + 1;
    idx_end   = idx_start + Nseg - 1;
    acc_win   = acc65(idx_start:idx_end);
    force_win = force65(idx_start:idx_end);
    
    % Correct H1 Spectral Densities: Input = Acc, Output = Force
    G_AA_all(:, iWin) = autocross(acc_win, acc_win, fsamp, Nseg, Win);
    G_AF_all(:, iWin) = autocross(acc_win, force_win, fsamp, Nseg, Win);
    G_FF_all(:, iWin) = autocross(force_win, force_win, fsamp, Nseg, Win);
end

% Average spectra
G_AA_avg = mean(G_AA_all, 2);
G_AF_avg = mean(G_AF_all, 2);
G_FF_avg = mean(G_FF_all, 2);
gamma2_65_avg = (abs(G_AF_avg).^2) ./ (G_AA_avg .* G_FF_avg);

% Compute measured H1 Apparent Mass [N / (m/s^2)]
H_app_65_avg = G_AF_avg ./ G_AA_avg;

% Select active frequency bounds
idx_freq    = (f >= lb) & (f <= ub);
f_plot      = f(idx_freq);
H_app_65_plot = H_app_65_avg(idx_freq);
gamma2_65_plot = gamma2_65_avg(idx_freq);

% Apply physical clamp mass cancellation
H_app_65_true = H_app_65_plot + m_clamp;

% Convert Apparent Mass to Dynamic Stiffness (N/m) 
omega = 2 * pi * f_plot; 
G_SY_65_plot = H_app_65_true .* (-omega.^2);

% Plot final averaged FRF

% figure
% semilogy(f_plot, abs(H1_65_plot), 'LineWidth', 1.5)
% grid on
% xlabel('Frequency [Hz]')
% ylabel('|H_1(f)| [m/s^2/N]')
% title('Averaged FRF - H1 Estimator for 6.5 A')
% xlim([lb ub])

%% Current 7 A

data7 = load('Selected_data/Current7.0/test@3861708838.mat');
acc7   = data7.Dati(:,1)/sensor_sens;       % Acceleration [m/s^2]
force7 = data7.Dati(:,2)/actuator_sens;     % Force [N]

nSamples = length(acc7);
nWindows = floor((nSamples - Nseg)/step) + 1;

% Initialize using the acceleration signal for frequency vector
[~, f] = autocross(acc7(1:Nseg), acc7(1:Nseg), fsamp, Nseg, Win);
f = f(:);
nFreq = length(f);

G_AA_all = zeros(nFreq, nWindows);
G_AF_all = zeros(nFreq, nWindows);
G_FF_all = zeros(nFreq, nWindows);

for iWin = 1:nWindows
    idx_start = (iWin-1)*step + 1;
    idx_end   = idx_start + Nseg - 1;
    acc_win   = acc7(idx_start:idx_end);
    force_win = force7(idx_start:idx_end);
    
    % Correct H1 Spectral Densities: Input = Acc, Output = Force
    G_AA_all(:, iWin) = autocross(acc_win, acc_win, fsamp, Nseg, Win);
    G_AF_all(:, iWin) = autocross(acc_win, force_win, fsamp, Nseg, Win);
    G_FF_all(:, iWin) = autocross(force_win, force_win, fsamp, Nseg, Win);
end

% Average spectra
G_AA_avg = mean(G_AA_all, 2);
G_AF_avg = mean(G_AF_all, 2);
G_FF_avg = mean(G_FF_all, 2);
gamma2_7_avg = (abs(G_AF_avg).^2) ./ (G_AA_avg .* G_FF_avg);

% Compute measured H1 Apparent Mass [N / (m/s^2)]
H_app_7_avg = G_AF_avg ./ G_AA_avg;

% Select active frequency bounds
idx_freq    = (f >= lb) & (f <= ub);
f_plot      = f(idx_freq);
H_app_7_plot = H_app_7_avg(idx_freq);
gamma2_7_plot = gamma2_7_avg(idx_freq);

% Apply physical clamp mass cancellation
H_app_7_true = H_app_7_plot + m_clamp;

% Convert Apparent Mass to Dynamic Stiffness (N/m) 
omega = 2 * pi * f_plot; 
G_SY_7_plot = H_app_7_true .* (-omega.^2);

% Plot final averaged FRF

% figure
% semilogy(f_plot, abs(H1_7_plot), 'LineWidth', 1.5)
% grid on
% xlabel('Frequency [Hz]')
% ylabel('|H_1(f)| [m/s^2/N]')
% title('Averaged FRF - H1 Estimator for 7 A')
% xlim([lb ub])

%% Current 7.5 A

data75 = load('Selected_data/Current7.5/test@3861709727.mat');
acc75   = data75.Dati(:,1)/sensor_sens;     % Acceleration [m/s^2]
force75 = data75.Dati(:,2)/actuator_sens;   % Force [N]

nSamples = length(acc75);
nWindows = floor((nSamples - Nseg)/step) + 1;

% Initialize using the acceleration signal for frequency vector
[~, f] = autocross(acc75(1:Nseg), acc75(1:Nseg), fsamp, Nseg, Win);
f = f(:);
nFreq = length(f);

G_AA_all = zeros(nFreq, nWindows);
G_AF_all = zeros(nFreq, nWindows);
G_FF_all = zeros(nFreq, nWindows);

for iWin = 1:nWindows
    idx_start = (iWin-1)*step + 1;
    idx_end   = idx_start + Nseg - 1;
    acc_win   = acc75(idx_start:idx_end);
    force_win = force75(idx_start:idx_end);
    
    % Correct H1 Spectral Densities: Input = Acc, Output = Force
    G_AA_all(:, iWin) = autocross(acc_win, acc_win, fsamp, Nseg, Win);
    G_AF_all(:, iWin) = autocross(acc_win, force_win, fsamp, Nseg, Win);
    G_FF_all(:, iWin) = autocross(force_win, force_win, fsamp, Nseg, Win);
end

% Average spectra
G_AA_avg = mean(G_AA_all, 2);
G_AF_avg = mean(G_AF_all, 2);
G_FF_avg = mean(G_FF_all, 2);
gamma2_75_avg = (abs(G_AF_avg).^2) ./ (G_AA_avg .* G_FF_avg);

% Compute measured H1 Apparent Mass [N / (m/s^2)]
H_app_75_avg = G_AF_avg ./ G_AA_avg;

% Select active frequency bounds
idx_freq    = (f >= lb) & (f <= ub);
f_plot      = f(idx_freq);
H_app_75_plot = H_app_75_avg(idx_freq);
gamma2_75_plot = gamma2_75_avg(idx_freq);

% Apply physical clamp mass cancellation
H_app_75_true = H_app_75_plot + m_clamp;

% Convert Apparent Mass to Dynamic Stiffness (N/m) 
omega = 2 * pi * f_plot; 
G_SY_75_plot = H_app_75_true .* (-omega.^2);

% Plot final averaged FRF

% figure
% semilogy(f_plot, abs(H1_75_plot), 'LineWidth', 1.5)
% grid on
% xlabel('Frequency [Hz]')
% ylabel('|H_1(f)| [m/s^2/N]')
% title('Averaged FRF - H1 Estimator for 7.5 A')
% xlim([lb ub])

%% Current 8.5 A

data85 = load('Selected_data/Current8.5/test@3861710602.mat');
acc85   = data85.Dati(:,1)/sensor_sens;     % Acceleration [m/s^2]
force85 = data85.Dati(:,2)/actuator_sens;   % Force [N]

nSamples = length(acc85);
nWindows = floor((nSamples - Nseg)/step) + 1;

% Initialize using the acceleration signal for frequency vector
[~, f] = autocross(acc85(1:Nseg), acc85(1:Nseg), fsamp, Nseg, Win);
f = f(:);
nFreq = length(f);

G_AA_all = zeros(nFreq, nWindows);
G_AF_all = zeros(nFreq, nWindows);
G_FF_all = zeros(nFreq, nWindows);

for iWin = 1:nWindows
    idx_start = (iWin-1)*step + 1;
    idx_end   = idx_start + Nseg - 1;
    acc_win   = acc85(idx_start:idx_end);
    force_win = force85(idx_start:idx_end);
    
    % Correct H1 Spectral Densities: Input = Acc, Output = Force
    G_AA_all(:, iWin) = autocross(acc_win, acc_win, fsamp, Nseg, Win);
    G_AF_all(:, iWin) = autocross(acc_win, force_win, fsamp, Nseg, Win);
    G_FF_all(:, iWin) = autocross(force_win, force_win, fsamp, Nseg, Win);
end

% Average spectra
G_AA_avg = mean(G_AA_all, 2);
G_AF_avg = mean(G_AF_all, 2);
G_FF_avg = mean(G_FF_all, 2);
gamma2_85_avg = (abs(G_AF_avg).^2) ./ (G_AA_avg .* G_FF_avg);

% Compute measured H1 Apparent Mass [N / (m/s^2)]
H_app_85_avg = G_AF_avg ./ G_AA_avg;

% Select active frequency bounds
idx_freq    = (f >= lb) & (f <= ub);
f_plot      = f(idx_freq);
H_app_85_plot = H_app_85_avg(idx_freq);
gamma2_85_plot = gamma2_85_avg(idx_freq);

% Apply physical clamp mass cancellation
H_app_85_true = H_app_85_plot + m_clamp;

% Convert Apparent Mass to Dynamic Stiffness (N/m) 
omega = 2 * pi * f_plot; 
G_SY_85_plot = H_app_85_true .* (-omega.^2);
% Plot final averaged FRF

% figure
% semilogy(f_plot, abs(H1_85_plot), 'LineWidth', 1.5)
% grid on
% xlabel('Frequency [Hz]')
% ylabel('|H_1(f)| [m/s^2/N]')
% title('Averaged FRF - H1 Estimator for 8.5 A')
% xlim([lb ub])

%% Current 9.5 A

data95 = load('Selected_data/Current9.5/test@3861711600.mat');
acc95   = data95.Dati(:,1)/sensor_sens;     % Acceleration [m/s^2]
force95 = data95.Dati(:,2)/actuator_sens;   % Force [N]

nSamples = length(acc95);
nWindows = floor((nSamples - Nseg)/step) + 1;

% Initialize using the acceleration signal for frequency vector
[~, f] = autocross(acc95(1:Nseg), acc95(1:Nseg), fsamp, Nseg, Win);
f = f(:);
nFreq = length(f);

G_AA_all = zeros(nFreq, nWindows);
G_AF_all = zeros(nFreq, nWindows);
G_FF_all = zeros(nFreq, nWindows);

for iWin = 1:nWindows
    idx_start = (iWin-1)*step + 1;
    idx_end   = idx_start + Nseg - 1;
    acc_win   = acc95(idx_start:idx_end);
    force_win = force95(idx_start:idx_end);
    
    % Correct H1 Spectral Densities: Input = Acc, Output = Force
    G_AA_all(:, iWin) = autocross(acc_win, acc_win, fsamp, Nseg, Win);
    G_AF_all(:, iWin) = autocross(acc_win, force_win, fsamp, Nseg, Win);
    G_FF_all(:, iWin) = autocross(force_win, force_win, fsamp, Nseg, Win);
end

% Average spectra
G_AA_avg = mean(G_AA_all, 2);
G_AF_avg = mean(G_AF_all, 2);
G_FF_avg = mean(G_FF_all, 2);
gamma2_95_avg = (abs(G_AF_avg).^2) ./ (G_AA_avg .* G_FF_avg);

% Compute measured H1 Apparent Mass [N / (m/s^2)]
H_app_95_avg = G_AF_avg ./ G_AA_avg;

% Select active frequency bounds
idx_freq    = (f >= lb) & (f <= ub);
f_plot      = f(idx_freq);
H_app_95_plot = H_app_95_avg(idx_freq);
gamma2_95_plot = gamma2_95_avg(idx_freq);

% Apply physical clamp mass cancellation
H_app_95_true = H_app_95_plot + m_clamp;

% Convert Apparent Mass to Dynamic Stiffness (N/m) 
omega = 2 * pi * f_plot; 
G_SY_95_plot = H_app_95_true .* (-omega.^2);

% Plot final averaged FRF

% figure
% semilogy(f_plot, abs(H1_95_plot), 'LineWidth', 1.5)
% grid on
% xlabel('Frequency [Hz]')
% ylabel('|H_1(f)| [m/s^2/N]')
% title('Averaged FRF - H1 Estimator for 9.5 A')
% xlim([lb ub])

%% Coherence evaluation

figure;
hold on; grid on; box on;
plot(f_plot, gamma2_0_plot,  'k-',  'LineWidth',1);
plot(f_plot, gamma2_05_plot, 'Color',[0.850 0.325 0.098], 'LineWidth',1.2);
plot(f_plot, gamma2_1_plot,  'Color',[0.929 0.694 0.125], 'LineWidth',1.2);
plot(f_plot, gamma2_2_plot,  'Color',[0.494 0.184 0.556], 'LineWidth',1.2);
plot(f_plot, gamma2_3_plot,  'Color',[0.466 0.674 0.188], 'LineWidth',1.2);
plot(f_plot, gamma2_4_plot,  'Color',[0.301 0.745 0.933], 'LineWidth',1.2);
plot(f_plot, gamma2_5_plot,  'Color',[0.635 0.078 0.184], 'LineWidth',1.2);
plot(f_plot, gamma2_6_plot,  'Color',[0.000 0.447 0.741], 'LineWidth',1.2);
plot(f_plot, gamma2_65_plot, 'Color',[0.750 0.750 0.000], 'LineWidth',1.2);
plot(f_plot, gamma2_7_plot,  'Color',[0.000 0.500 0.500], 'LineWidth',1.2);
plot(f_plot, gamma2_75_plot, 'Color',[0.500 0.000 0.500], 'LineWidth',1.2);
plot(f_plot, gamma2_85_plot, 'Color',[0.800 0.400 0.000], 'LineWidth',1.2);
plot(f_plot, gamma2_95_plot, 'Color',[0.200 0.200 0.200], 'LineWidth',1.2);
legend({'0 A','0.5 A','1 A','2 A','3 A','4 A','5 A','6 A','6.5 A','7 A','7.5 A','8.5 A','9.5 A'}, 'NumColumns',2, 'Location','southwest');
ylim([0.7 1.05]); % Cap the axis just above 1 to see the top clearly
xlabel('Frequency [Hz]');
ylabel('Coherence \gamma^2');
title('Coherence Function Comparison');

%% Final comparison plot for acceleration

figure;
semilogy(f_plot, abs(H_app_0_true), 'k-', 'LineWidth',2, 'Marker','o', 'MarkerIndices',1:200:length(f_plot), 'MarkerSize',2, 'MarkerFaceColor','k')
hold on; grid on;
semilogy(f_plot, abs(H_app_05_true), 'LineStyle','-',  'Color',[0.902, 0.098, 0.294], 'LineWidth',1.5)
semilogy(f_plot, abs(H_app_1_true),  'LineStyle','--', 'Color',[0.235, 0.706, 0.294], 'LineWidth',1.5)
semilogy(f_plot, abs(H_app_2_true),  'LineStyle','-',  'Color',[0.000, 0.510, 0.784], 'LineWidth',1.5)
semilogy(f_plot, abs(H_app_3_true),  'LineStyle','-.', 'Color',[1.000, 0.882, 0.098], 'LineWidth',1.5)
semilogy(f_plot, abs(H_app_4_true),  'LineStyle','-',  'Color',[0.569, 0.118, 0.706], 'LineWidth',1.5)
semilogy(f_plot, abs(H_app_5_true),  'LineStyle','-',  'Color',[0.961, 0.510, 0.188], 'LineWidth',1.5)
semilogy(f_plot, abs(H_app_6_true),  'LineStyle','--', 'Color',[0.275, 0.941, 0.941], 'LineWidth',1.5)
semilogy(f_plot, abs(H_app_65_true), 'LineStyle','-',  'Color',[0.941, 0.196, 0.902], 'LineWidth',1.5)
semilogy(f_plot, abs(H_app_7_true),  'LineStyle','-',  'Color',[0.824, 0.961, 0.235], 'LineWidth',1.5)
semilogy(f_plot, abs(H_app_75_true), 'LineStyle','-',  'Color',[0.980, 0.745, 0.745], 'LineWidth',1.5)
semilogy(f_plot, abs(H_app_85_true), 'LineStyle','--', 'Color',[0.667, 0.431, 0.157], 'LineWidth',1.5)
semilogy(f_plot, abs(H_app_95_true), 'LineStyle','-',  'Color',[0.000, 0.502, 0.502], 'LineWidth',1.5)

xlabel('Frequency [Hz]')
ylabel('|H_{app}(f)| [N / (m/s^2)]')
title('Apparent Mass FRF Comparison for Different Currents (Force / Acceleration)')
xlim([lb ub])
legend('0.00 A','0.50 A','1.00 A','2.00 A','3.00 A','4.00 A','5.00 A','6.00 A','6.50 A','7.00 A','7.50 A','8.50 A','9.50 A','Location','northeast')

%% Transformation in displacement domain

figure;
semilogy(f_plot, abs(G_SY_0_plot), 'k-', 'LineWidth',2, 'Marker','o', 'MarkerIndices',1:200:length(f_plot), 'MarkerSize',2, 'MarkerFaceColor','k')
hold on; grid on;
semilogy(f_plot, abs(G_SY_05_plot), 'LineStyle','-',  'Color',[0.902, 0.098, 0.294], 'LineWidth',1.5)
semilogy(f_plot, abs(G_SY_1_plot),  'LineStyle','--', 'Color',[0.235, 0.706, 0.294], 'LineWidth',1.5)
semilogy(f_plot, abs(G_SY_2_plot),  'LineStyle','-',  'Color',[0.000, 0.510, 0.784], 'LineWidth',1.5)
semilogy(f_plot, abs(G_SY_3_plot),  'LineStyle','-.', 'Color',[1.000, 0.882, 0.098], 'LineWidth',1.5)
semilogy(f_plot, abs(G_SY_4_plot),  'LineStyle','-',  'Color',[0.569, 0.118, 0.706], 'LineWidth',1.5)
semilogy(f_plot, abs(G_SY_5_plot),  'LineStyle','-',  'Color',[0.961, 0.510, 0.188], 'LineWidth',1.5)
semilogy(f_plot, abs(G_SY_6_plot),  'LineStyle','--', 'Color',[0.275, 0.941, 0.941], 'LineWidth',1.5)
semilogy(f_plot, abs(G_SY_65_plot), 'LineStyle','-',  'Color',[0.941, 0.196, 0.902], 'LineWidth',1.5)
semilogy(f_plot, abs(G_SY_7_plot),  'LineStyle','-',  'Color',[0.824, 0.961, 0.235], 'LineWidth',1.5)
semilogy(f_plot, abs(G_SY_75_plot), 'LineStyle','-',  'Color',[0.980, 0.745, 0.745], 'LineWidth',1.5)
semilogy(f_plot, abs(G_SY_85_plot), 'LineStyle','--', 'Color',[0.667, 0.431, 0.157], 'LineWidth',1.5)
semilogy(f_plot, abs(G_SY_95_plot), 'LineStyle','-',  'Color',[0.000, 0.502, 0.502], 'LineWidth',1.5)

xlabel('Frequency [Hz]')
ylabel('|G_{SY}(f)| [N/m]')
title('True Dynamic Stiffness Comparison for Different Currents')
xlim([lb ub])
legend('0.00 A','0.50 A','1.00 A','2.00 A','3.00 A','4.00 A','5.00 A','6.00 A','6.50 A','7.00 A','7.50 A','8.50 A','9.50 A','Location','northeast')

%% Modal extraction - Peak picking

currents = [0 0.5 1 2 3 4 5 6 6.5 7 7.5 8.5 9.5];

% Update to use the corrected, mass-cancelled Apparent Mass variables
H_all = {H_app_0_true, H_app_05_true, H_app_1_true, H_app_2_true, H_app_3_true, ...
         H_app_4_true, H_app_5_true, H_app_6_true, H_app_65_true, H_app_7_true, ...
         H_app_75_true, H_app_85_true, H_app_95_true};

fn1 = zeros(size(currents));
fn2 = zeros(size(currents));
h1 = zeros(size(currents));
h2 = zeros(size(currents));

% Choose approximate frequency ranges for mode 1 and mode 2
% (Adjust these tightly around your specific physical peak zones)
mode1_range = [15 45];    
mode2_range = [150 350]; 

for k = 1:length(currents)
    Hmag = abs(H_all{k});
    
    % ----- Mode 1 -----
    idx1 = f_plot >= mode1_range(1) & f_plot <= mode1_range(2); 
    f1_range = f_plot(idx1);
    H1_range = Hmag(idx1);
    
    [peak1, loc1] = max(H1_range);
    fn1(k) = f1_range(loc1);
    
    half1 = peak1 / sqrt(2);
    left_idx = find(H1_range(1:loc1) <= half1, 1, 'last');
    right_idx = loc1 - 1 + find(H1_range(loc1:end) <= half1, 1, 'first');
    
    if ~isempty(left_idx) && ~isempty(right_idx)
        f_left = f1_range(left_idx);
        f_right = f1_range(right_idx);
        h1(k) = (f_right - f_left) / (2 * fn1(k));
    else
        h1(k) = NaN;
    end
    
    % ----- Mode 2 -----
    idx2 = f_plot >= mode2_range(1) & f_plot <= mode2_range(2);
    f2_range = f_plot(idx2);
    H2_range = Hmag(idx2);
    
    [peak2, loc2] = max(H2_range);
    fn2(k) = f2_range(loc2);
    
    half2 = peak2 / sqrt(2);
    left_idx = find(H2_range(1:loc2) <= half2, 1, 'last');
    right_idx = loc2 - 1 + find(H2_range(loc2:end) <= half2, 1, 'first');
    
    if ~isempty(left_idx) && ~isempty(right_idx)
        f_left = f2_range(left_idx);
        f_right = f2_range(right_idx);
        h2(k) = (f_right - f_left) / (2 * fn2(k));
    else
        h2(k) = NaN;
    end
end

%% Modal extraction by SDOF curve fitting using lsqnonlin

mode1_fit_range = [15 45];     
mode2_fit_range = [150 350];    

fit_ranges = [mode1_fit_range;
              mode2_fit_range];

fn_fit   = zeros(length(currents), 2);
h_fit = zeros(length(currents), 2);
p_all_fit = cell(length(currents), 2);

opts = optimoptions('lsqnonlin', ...
    'Display','off', ...
    'MaxFunctionEvaluations',5000, ...
    'MaxIterations',1000);

for k = 1:length(currents)

    H_current = H_all{k};

    for m = 1:2

        % Broad search range for current mode
        search_range = fit_ranges(m,:);

        idx_search = f_plot >= search_range(1) & f_plot <= search_range(2);
        f_search = f_plot(idx_search);
        H_search = H_current(idx_search);

        % Peak detection inside broad range
        [~, imax_search] = max(abs(H_search));
        fn_peak = f_search(imax_search);

        % Fitting range: +/- 10% around detected peak
        fit_range = [0.9*fn_peak, 1.1*fn_peak];

        % Limit fitting range inside original search interval
        fit_range(1) = max(fit_range(1), search_range(1));
        fit_range(2) = min(fit_range(2), search_range(2));

        idx_fit = f_plot >= fit_range(1) & f_plot <= fit_range(2);
        f_fit = f_plot(idx_fit);
        H_fit = H_current(idx_fit);
        w = 2*pi*f_fit;

        % Initial guess
        fn0 = fn_peak;
        wn0 = 2*pi*fn0;
        zeta0 = 0.01;

        A0 = max(abs(H_fit)) * (2*zeta0*wn0^2);
        C0 = mean(H_fit);

        % p = [wn, zeta, real(A), imag(A), real(C), imag(C)]
        p0 = [wn0, zeta0, real(A0), imag(A0), real(C0), imag(C0)];

        modelFun = @(p,w) ...
            (p(5) + 1i*p(6)) + ...
            (p(3) + 1i*p(4)) ./ ...
            (p(1)^2 - w.^2 + 1i*2*p(2)*p(1).*w);

        resFun = @(p) [real(modelFun(p,w) - H_fit);
                       imag(modelFun(p,w) - H_fit)];

        lb_fit = [2*pi*fit_range(1), 0.0001, -Inf, -Inf, -Inf, -Inf];
        ub_fit = [2*pi*fit_range(2), 0.30,    Inf,  Inf,  Inf,  Inf];

        p_est = lsqnonlin(resFun, p0, lb_fit, ub_fit, opts);

        p_all_fit{k,m} = p_est;
        fn_fit(k,m)   = p_est(1)/(2*pi);
        h_fit(k,m) = p_est(2);

    end
end

fn1_fit = fn_fit(:,1);
fn2_fit = fn_fit(:,2);
h1_fit = h_fit(:,1);
h2_fit = h_fit(:,2);

results_table = table(currents(:), fn1_fit, h1_fit, fn2_fit, h2_fit, ...
    'VariableNames', {'Current_A','fn1_Hz','zeta1','fn2_Hz','zeta2'});

disp(results_table)

%% Plotting Extracted Parameters vs Current

% Natural Frequencies vs Current ( mode 1)
figure;
plot(currents, fn1_fit, 'o-', 'Color', [0.000, 0.510, 0.784], 'LineWidth', 2, 'MarkerFaceColor', [0.000, 0.510, 0.784]);
grid on; box on;
xlabel('Current [A]');
ylabel('Natural Frequency [Hz]');
title('Tracked Natural Frequencies vs Control Current for mode 1');


% Natural Frequencies vs Current ( mode 2)
figure;
plot(currents, fn2_fit, 's-', 'Color', [0.902, 0.098, 0.294], 'LineWidth', 2, 'MarkerFaceColor', [0.902, 0.098, 0.294]);
grid on; box on;
xlabel('Current [A]');
ylabel('Natural Frequency [Hz]');
title('Tracked Natural Frequencies vs Control Current for mode 2');


% Plot Damping Ratio Mode 1 vs Current
figure;
plot(currents, h1_fit * 100, 'o-', 'Color', [0.235, 0.706, 0.294], 'LineWidth', 2);
hold on
plot(currents, h1 * 100, 'k-o');
grid on; box on;

legend('Minimization','Peak peaking')
xlabel('Current [A]');
ylabel('Damping Ratio \zeta [%]');
title('Mode 1 Damping Ratio Evolution vs Current');

% Plot Damping Ratio Mode 2 vs Current
figure;
plot(currents, h2_fit * 100,'o-', 'Color', [0.235, 0.706, 0.294], 'LineWidth', 2);
hold on
plot(currents, h2 * 100, 'k-o');
legend('Minimization','Peak peaking')

grid on; box on;
xlabel('Current [A]');
ylabel('Damping Ratio \zeta [%]');
title('Mode 2 Damping Ratio Evolution vs Current');

%% Comparison in a single figure
% Original vs Peak Picking vs Minimization
% Mode 2 only

for k = 1:length(currents)
    H_original = H_all{k};
    
    % Create a single wide figure for the current step
    figure('Name', ['FRF Reconstruction - Current = ', num2str(currents(k)), ' A']); 
    
    % ================= MODE 1 Subplot =================
    subplot(1, 2, 1);
    % Retrieve saved optimized vector for Mode 1 (m=1)
    p_opt1 = p_all_fit{k, 1}; 
    fn_opt1 = p_opt1(1)/(2*pi);

    % Plot range: +/- 10% around fitted eigenfrequency
    plot_range1 = [0.9*fn_opt1, 1.1*fn_opt1];

    idx_mode1 = f_plot >= plot_range1(1) & f_plot <= plot_range1(2);
    f_rec1 = f_plot(idx_mode1);
    w_rec1 = 2*pi*f_rec1;
    H_original_mode1 = H_original(idx_mode1);
    
    % Reconstruct Mode 1 with residuals
    H_rec_mode1_min = (p_opt1(5) + 1i*p_opt1(6)) + ...
                      (p_opt1(3) + 1i*p_opt1(4)) ./ ...
                      (p_opt1(1)^2 - w_rec1.^2 + 1i*2*p_opt1(2)*p_opt1(1).*w_rec1);
                      
    semilogy(f_rec1, abs(H_original_mode1), 'k-', 'LineWidth', 2);
    hold on; grid on; box on;
    semilogy(f_rec1, abs(H_rec_mode1_min), 'r--', 'LineWidth', 2);
    xline(fn_opt1, 'b:', 'LineWidth', 1.5);
    
    xlabel('Frequency [Hz]');
    ylabel('|H_{app}(f)| [N / (m/s^2)]');
    title(sprintf('Mode 1 - %.1f A', currents(k)));
    legend('Experimental', 'Reconstructed Model', 'f_{n,1}', 'Location','best');
    xlim(plot_range1);
    
    % ================= MODE 2 Subplot =================
    subplot(1, 2, 2);
    
    % Retrieve saved optimized vector for Mode 2 (m=2)
    p_opt2 = p_all_fit{k, 2}; 
    fn_opt2 = p_opt2(1)/(2*pi);

    plot_range2 = [0.9*fn_opt2, 1.1*fn_opt2];

    idx_mode2 = f_plot >= plot_range2(1) & f_plot <= plot_range2(2);
    f_rec2 = f_plot(idx_mode2);
    w_rec2 = 2*pi*f_rec2;
    H_original_mode2 = H_original(idx_mode2);
    
    % Reconstruct Mode 2 with residuals
    H_rec_mode2_min = (p_opt2(5) + 1i*p_opt2(6)) + ...
                      (p_opt2(3) + 1i*p_opt2(4)) ./ ...
                      (p_opt2(1)^2 - w_rec2.^2 + 1i*2*p_opt2(2)*p_opt2(1).*w_rec2);
                      
    semilogy(f_rec2, abs(H_original_mode2), 'k-', 'LineWidth', 2);
    hold on; grid on; box on;
    semilogy(f_rec2, abs(H_rec_mode2_min), 'r--', 'LineWidth', 2);
    xline(fn_opt2, 'b:', 'LineWidth', 1.5);
    
    xlabel('Frequency [Hz]');
    ylabel('|H_{app}(f)| [N / (m/s^2)]');
    title(sprintf('Mode 2 - %.1f A', currents(k)));
    legend('Experimental', 'Reconstructed Model', 'f_{n,2}', 'Location','best');
    xlim(plot_range2);
    
    % Main overarching figure title
    sgtitle(sprintf('FRF SDOF Minimization Reconstruction - Current = %.1f A', currents(k)), 'FontSize', 14, 'FontWeight', 'bold');
end

%% MAC between Experimental and Reconstructed FRF (to validate fitting quality)

% Preallocate MAC matrices: rows = currents, columns = modes (1 and 2)
MAC_values = zeros(length(currents), 2);

for k = 1:length(currents)
    
    H_current = H_all{k};
    
    for m = 1:2
        
        % Get the fitted parameters for this mode
        p_opt = p_all_fit{k, m};
        fn_opt = p_opt(1) / (2*pi);
        
        % Define the frequency band for comparison (same as fitting range used)
        % I use +/- 10% around the detected natural frequency
        band_limits = [0.9 * fn_opt, 1.1 * fn_opt];
        
        % Get indices of frequencies strictly inside this band
        idx_band = (f_plot >= band_limits(1)) & (f_plot <= band_limits(2));
        
        % Extract experimental and reconstructed complex vectors
        H_exp_band = H_current(idx_band);
        
        % Reconstruct the model for these exact frequencies
        w_band = 2 * pi * f_plot(idx_band);
        H_rec_band = (p_opt(5) + 1i*p_opt(6)) + ...
                     (p_opt(3) + 1i*p_opt(4)) ./ ...
                     (p_opt(1)^2 - w_band.^2 + 1i*2*p_opt(2)*p_opt(1).*w_band);
        
        % Compute MAC (ensuring vectors are column vectors)
        num = abs( H_exp_band' * H_rec_band )^2;
        den = ( H_exp_band' * H_exp_band ) * ( H_rec_band' * H_rec_band );
        MAC_values(k, m) = num / den;
        
    end
end

% Display results in a table
MAC_table = table(currents(:), MAC_values(:,1), MAC_values(:,2), ...
    'VariableNames', {'Current_A', 'MAC_Mode1', 'MAC_Mode2'});
disp('MAC between Experimental and Reconstructed FRF (1 = perfect fit):');
disp(MAC_table);

%% Optional: Plot MAC vs Current
figure;
subplot(1,2,1);
plot(currents, MAC_values(:,1), 'o-', 'LineWidth', 2, 'Color', [0, 0.5, 0.8]);
grid on; box on;
xlabel('Current [A]'); ylabel('MAC');
title('Mode 1 - Fit Quality'); ylim([0.95 1.01]);

subplot(1,2,2);
plot(currents, MAC_values(:,2), 's-', 'LineWidth', 2, 'Color', [0.9, 0.1, 0.3]);
grid on; box on;
xlabel('Current [A]'); ylabel('MAC');
title('Mode 2 - Fit Quality'); ylim([0.95 1.01]);

sgtitle('MAC between Experimental and Reconstructed FRFs', 'FontWeight', 'bold');


%% Curve fitting for controller design (Piecewise Approach) 

[min_f1, min_idx1] = min(fn1_fit);
P_inv1_M = polyfit(fn1_fit(1:min_idx1), currents(1:min_idx1), 2);
P_inv1_A = polyfit(fn1_fit(min_idx1:end), currents(min_idx1:end), 3);

[min_f2, min_idx2] = min(fn2_fit);
P_inv2_M = polyfit(fn2_fit(1:min_idx2), currents(1:min_idx2), 2);
P_inv2_A = polyfit(fn2_fit(min_idx2:end), currents(min_idx2:end), 3);

max_f1_M = fn1_fit(1); 
max_f2_M = fn2_fit(1);

P_forward2 = polyfit(currents, fn2_fit, 4);
P_forward1 = polyfit(currents, fn1_fit, 4);


%% 

plotNyquistCircle(H_app_0_true,...
                  f_plot,...
                  [14 60],...
                  'Mode 1 - Current = 0 A');
plotNyquistCircle(H_app_0_true,...
                  f_plot,...
                  [200 450],...
                  'Mode 2 - Current = 0 A');

plotNyquistCircle(H_app_05_true,...
                  f_plot,...
                  [14 60],...
                  'Mode 1 - Current = 0.5 A');
plotNyquistCircle(H_app_05_true,...
                  f_plot,...
                  [200 450],...
                  'Mode 2 - Current = 0.5 A');

plotNyquistCircle(H_app_1_true,...
                  f_plot,...
                  [14 60],...
                  'Mode 1 - Current = 1 A');
plotNyquistCircle(H_app_1_true,...
                  f_plot,...
                  [200 450],...
                  'Mode 2 - Current = 1 A');
plotNyquistCircle(H_app_2_true,...
                  f_plot,...
                  [14 60],...
                  'Mode 1 - Current = 2 A');
plotNyquistCircle(H_app_2_true,...
                  f_plot,...
                  [200 450],...
                  'Mode 2 - Current = 2 A');
plotNyquistCircle(H_app_3_true,...
                  f_plot,...
                  [14 60],...
                  'Mode 1 - Current = 3 A');
plotNyquistCircle(H_app_3_true,...
                  f_plot,...
                  [200 450],...
                  'Mode 2 - Current = 3 A');
plotNyquistCircle(H_app_4_true,...
                  f_plot,...
                  [14 60],...
                  'Mode 1 - Current = 4 A');
plotNyquistCircle(H_app_4_true,...
                  f_plot,...
                  [200 450],...
                  'Mode 2 - Current = 4 A');
plotNyquistCircle(H_app_5_true,...
                  f_plot,...
                  [14 60],...
                  'Mode 1 - Current = 5 A');
plotNyquistCircle(H_app_5_true,...
                  f_plot,...
                  [200 450],...
                  'Mode 2 - Current = 5 A');
plotNyquistCircle(H_app_6_true,...
                  f_plot,...
                  [14 60],...
                  'Mode 1 - Current = 6 A');
plotNyquistCircle(H_app_6_true,...
                  f_plot,...
                  [200 450],...
                  'Mode 2 - Current = 6 A');
plotNyquistCircle(H_app_65_true,...
                  f_plot,...
                  [14 60],...
                  'Mode 1 - Current = 6.5 A');
plotNyquistCircle(H_app_65_true,...
                  f_plot,...
                  [200 450],...
                  'Mode 2 - Current = 6.5 A');
plotNyquistCircle(H_app_7_true,...
                  f_plot,...
                  [14 60],...
                  'Mode 1 - Current = 7 A');
plotNyquistCircle(H_app_7_true,...
                  f_plot,...
                  [200 450],...
                  'Mode 2 - Current = 7 A');
plotNyquistCircle(H_app_75_true,...
                  f_plot,...
                  [14 60],...
                  'Mode 1 - Current = 7.5 A');
plotNyquistCircle(H_app_75_true,...
                  f_plot,...
                  [200 450],...
                  'Mode 2 - Current = 7.5 A');
plotNyquistCircle(H_app_85_true,...
                  f_plot,...
                  [14 60],...
                  'Mode 1 - Current = 8.5 A');
plotNyquistCircle(H_app_85_true,...
                  f_plot,...
                  [200 450],...
                  'Mode 2 - Current = 8.5 A');
plotNyquistCircle(H_app_95_true,...
                  f_plot,...
                  [14 60],...
                  'Mode 1 - Current = 9.5 A');
plotNyquistCircle(H_app_95_true,...
                  f_plot,...
                  [200 450],...
                  'Mode 2 - Current = 9.5 A');

