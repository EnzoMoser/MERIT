% Parameters
Fs = 1e12;            % Sampling frequency (samples per second)
T = 1/Fs;             % Sampling period (seconds per sample)
L = 1001;             % Length of signal
t = (0:L)*T;   % Time vector

% Gaussian pulse parameters
sigma = 5e-11;         % Standard deviation of the Gaussian (seconds)

% Generate the Gaussian pulse
gauspulse = exp(-(t/sigma).^2);

% Plot the Gaussian pulse in the time domain
figure(10);
plot(t, gauspulse);
title('Gaussian Pulse Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;
% Plot the power spectrum
