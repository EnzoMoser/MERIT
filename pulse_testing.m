%Gaussian Pulse generator and testing
% uses the MERIT functions

%% Time and Frequency arrays


frequencies = [0:4.5e6:4.5e9]';
times = [0:10e-12:10e-9]';



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

%use time to frequency shifter to get gaussian pulse in frequency domain.
gausspulsefreq = merit.process.td2fd(gausspulse,times,frequencies);

gausspulsefreq_magnitude = abs(gausspulsefreq);

%modulate gausssian pulse with a 3 GHz sine wave
sinewave_mod = sin(2*pi*3e+9*times);
gausspulsemod = gausspulse .* sinewave_mod;

gausspulsemodfreq = merit.process.td2fd(gausspulsemod,times,frequencies);
gausspulsemodfreq_magnitude = abs(gausspulsemodfreq);
%% Plot gaussian pulse
figure(1)
subplot(2, 1, 1);
plot(times, gausspulse);
xlabel('Time (s)');
ylabel('value');
title(sprintf('Gaussian Pulse Time Domain'));

subplot(2, 1, 2);
plot(frequencies, gausspulsefreq_magnitude);
xlabel('Frequency (Hz)');
ylabel('vaule');
title(sprintf('Gaussian Pulse Frequency Domain'));


%% Plot modulated gaussian pulse
figure(2)
subplot(2, 1, 1);
plot(times, gausspulsemod);
xlabel('Time (s)');
ylabel('value');
title(sprintf('Modulated Gaussian Pulse Time Domain'));

subplot(2, 1, 2);
plot(frequencies, gausspulsemodfreq_magnitude);
xlabel('Frequency (Hz)');
ylabel('vaule');
title(sprintf('Modulated Gaussian Pulse Frequency Domain'));
