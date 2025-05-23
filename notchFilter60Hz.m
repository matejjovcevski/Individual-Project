function Hd = notchFilter60Hz
%NOTCHFILTER60HZ Returns a discrete-time filter object.

% MATLAB Code
% Generated by MATLAB(R) 23.2 and Signal Processing Toolbox 23.2.
% Generated on: 24-Feb-2025 21:08:32

% Equiripple Bandstop filter designed using the FIRPM function.

% All frequency values are in Hz.
Fs = 1000;  % Sampling Frequency

Fpass1 = 56;              % First Passband Frequency
Fstop1 = 59.5;            % First Stopband Frequency
Fstop2 = 60.5;            % Second Stopband Frequency
Fpass2 = 64;              % Second Passband Frequency
Dpass1 = 0.057501127785;  % First Passband Ripple
Dstop  = 0.0001;          % Stopband Attenuation
Dpass2 = 0.057501127785;  % Second Passband Ripple
dens   = 20;              % Density Factor

% Calculate the order from the parameters using FIRPMORD.
[N, Fo, Ao, W] = firpmord([Fpass1 Fstop1 Fstop2 Fpass2]/(Fs/2), [1 0 ...
                          1], [Dpass1 Dstop Dpass2]);

% Calculate the coefficients using the FIRPM function.
b  = firpm(N, Fo, Ao, W, {dens});
Hd = dfilt.dffir(b);

% [EOF]
