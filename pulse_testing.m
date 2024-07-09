% Getting started guide for MERIT
% A basic guide to:
%  - loading and visualising the sample data;
%  - processing signals using the MERIT functions;
%  - imaging with this toolbox.

%% Loading sample data
% Details of the breast phantoms used to collect the sample data
% are given in "Microwave Breast Imaging: experimental
% tumour phantoms for the evaluation of new breast cancer diagnosis
% systems", 2018 Biomed. Phys. Eng. Express 4 025036.
% The antenna locations, frequency points and scattered signals
% are given in the /data folder:
%   antenna_locations.csv: the antenna locations in metres;
%   frequencies.csv: the frequency points in Hertz;
%   channel_names.csv: the descriptions of the channels in the scattered data;
%   B0_P3_p000.csv: homogeneous breast phantom with an 11 mm diameter
%     tumour located at (15, 0, 35) mm.
%   B0_P5_p000.csv: homogeneous breast phantom with an 20 mm diameter
%     tumour located at (15, 0, 30) mm.
% For both phantoms, a second scan rotated by 36 degrees from the first
% was acquired for artefact removal:
% B0_P3_p036.csv and B0_P5_p036.csv respectively.

frequencies = [-2e9:4.5e6:2e9]';
antenna_locations = dlmread('data/antenna_locations.csv');
channel_names = dlmread('data/channel_names.csv');
times = [0:10e-12:10e-9]';
timesT = [0:10e-12:0.75e-9]';

scan1 = dlmread('data/B0_P3_p000.csv');
scan2 = dlmread('data/B0_P3_p036.csv');

%% Gaussian Pulse generation
% Parameters
Fs = 1e12;            % Sampling frequency (samples per second)
T = 1/Fs;             % Sampling period (seconds per sample)
L = 1000;             % Length of signal
t = (-L/2:L/2)*T;   % Time vector
% Gaussian pulse parameters
sigma = 5e-11;         % Standard deviation of the Gaussian (seconds)

% Generate the Gaussian pulse
gausspulseT = exp(-(t/sigma).^2);
gausspulse = gausspulseT';
gausspulsefreq = merit.process.td2fd(gausspulse,times,frequencies);

gausspulsefreq_magnitude = abs(gausspulsefreq);
%% Plot the acquired scans.
figure(1)
subplot(2, 1, 1);
plot(times, gausspulse);
xlabel('Time (s)');
ylabel('value');
legend('gaussian pulse time domian');

subplot(2, 1, 2);
plot(frequencies, gausspulsefreq_magnitude);
xlabel('Frequency (Hz)');
ylabel('vaule');
legend('gaussian pulse frequency domian');


%% Perform rotation subtraction
signals = scan1-scan2;

%% Plot artefact removed: channel 1
%subplot(2, 1, 1);
%plot(times, [data_channel1, signals(:,1)]);


%xlabel('time');
%ylabel('data');
%title(sprintf('Channel (%d, %d) Magnitude—Artefact removed', channel_names(1, :)));
%subplot(2, 1, 2);

%plot(frequencies, [channel1_phase, rotated_channel1_phase]);
%xlabel('Frequency (Hz)');
%ylabel('Phase (rad)');
%legend('Original Scan', 'Rotated Scan', 'Artefact removed');
%title(sprintf('Channel (%d, %d) Phase—Artefact removed', channel_names(1, :)));