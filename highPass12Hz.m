function Hd = highPass12Hz
%HIGHPASS12HZ Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 23.2 and Signal Processing Toolbox 23.2.
% Generated on: 31-Mar-2025 16:20:45

% Equiripple Highpass filter designed using the FIRPM function.

% All frequency values are in Hz.
Fs = 1000;  % Sampling Frequency

Fstop = 11;              % Stopband Frequency
Fpass = 12;              % Passband Frequency
Dstop = 0.0001;          % Stopband Attenuation
Dpass = 0.057501127785;  % Passband Ripple
dens  = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fstop, Fpass]/(Fs/2), [0 1], [Dstop, Dpass]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);

% [EOF]
