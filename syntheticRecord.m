function syntheticRecord(subStr, SNRStr)

% Initialise some strings
recordName = strcat('testSynthDB/sub0',subStr,'_snr',SNRStr,'dB_l1_c0_');
fStr = strcat(recordName,'fecg1');
mStr = strcat(recordName,'mecg');
n1Str = strcat(recordName,'noise1');
n2Str = strcat(recordName,'noise2');

% 0. Read Data
[fecg1DIR, ~, ~] = rdsamp(fStr,1);
[fecg2DIR, ~, ~] = rdsamp(fStr,2);
[fecg3DIR, ~, ~] = rdsamp(fStr,3);
[fecg4DIR, ~, ~] = rdsamp(fStr,4);

[mecg1, ~, ~] = rdsamp(mStr,1);
[mecg2, ~, ~] = rdsamp(mStr,2);
[mecg3, ~, ~] = rdsamp(mStr,3);
[mecg4, ~, ~] = rdsamp(mStr,4);

[noise11, ~, ~] = rdsamp(n1Str,1);
[noise12, ~, ~] = rdsamp(n1Str,2);
[noise13, ~, ~] = rdsamp(n1Str,3);
[noise14, ~, ~] = rdsamp(n1Str,4);

[noise21, ~, ~] = rdsamp(n2Str,1);
[noise22, ~, ~] = rdsamp(n2Str,2);
[noise23, ~, ~] = rdsamp(n2Str,3);
[noise24, fs, ~] = rdsamp(n2Str,4);

[time, m, atrtimed, annotd] = rdann(fStr,'qrs',[],1000000);
directImpulse = time;


sig1 = fecg1DIR+mecg1+noise11+noise21;
sig2 = fecg2DIR+mecg2+noise12+noise22;
sig3 = fecg3DIR+mecg3+noise13+noise23;
sig4 = fecg4DIR+mecg4+noise14+noise24;

X_raw = [sig1 sig2 sig3 sig4];
[ndt, ns] = size(X_raw);
arrayRows = 7; %SEVEN METHODS 
arrayCols = 2*ns + 3; % format: fecg 1-4, ica 1-4, + max PCA and max ICA
outputMat = zeros([arrayRows arrayCols]);
% Preprocessing
order = 633;
[X_bp, ~] = bandpassFilter(X_raw,BandPass2HzSamp250Hz,order,0); % Bandpass filtering, delay is not used throughout since signal is centred
X_im=FecgImpArtCanc(X_bp,fs); % Impulse artefact cancellation

%START LOOP
count=1;
templateSize = 20;

% !!! 1. cancelMecg !!!
%---------------------------------------------------------------------------%
% mECG cancellation
[fecg1, ~, ~] = cancelMecg(X_im(:,1), fs, templateSize,0,0);
[fecg2, ~, ~] = cancelMecg(X_im(:,2), fs, templateSize,0,0);
[fecg3, ~, ~] = cancelMecg(X_im(:,3), fs, templateSize,0,0);
[fecg4, ~, delayValue] = cancelMecg(X_im(:,4), fs, templateSize,0,0);

% ICA/Post-processing
X_fecg = [fecg1, fecg2, fecg3, fecg4];
X_fecg = transpose(X_fecg);
r = ns;
[X_ica, ~, ~, ~] = fastICA(X_fecg,r);
X_ica = transpose(X_ica);
X_fecg = transpose(X_fecg);
% 4. F1 Calculation
outputMat(count,1) = 1;
for j=1:ns
    [f1, ~] = compareDirect(X_fecg(:,j), directImpulse, delayValue, fs); %PCA F1s
    outputMat(count,j+1) = f1;
    [f1, ~] = compareDirect(X_ica(:,j), directImpulse, delayValue, fs); %ICA F1s
    outputMat(count,j+ns+1) = f1;
end
outputMat(count,2*ns+2) = max(outputMat(count,2:5));
outputMat(count,2*ns+3) = max(outputMat(count,6:9));
count=count+1;
%---------------------------------------------------------------------------%

% !!! 2. cancelMecgAverage !!!
%---------------------------------------------------------------------------%
% mECG cancellation

[fecg1, ~, ~] = cancelMecgAverage(X_im(:,1), fs, templateSize,0,0);
[fecg2, ~, ~] = cancelMecgAverage(X_im(:,2), fs, templateSize,0,0);
[fecg3, ~, ~] = cancelMecgAverage(X_im(:,3), fs, templateSize,0,0);
[fecg4, ~, delayValue] = cancelMecgAverage(X_im(:,4), fs, templateSize,0,0);

% ICA/Post-processing
X_fecg = [fecg1, fecg2, fecg3, fecg4];
X_fecg = transpose(X_fecg);
r = ns;
[X_ica, ~, ~, ~] = fastICA(X_fecg,r);
X_ica = transpose(X_ica);
X_fecg = transpose(X_fecg);
% F1 Calculation
outputMat(count,1) = 2;
for j=1:ns
    [f1, ~] = compareDirect(X_fecg(:,j), directImpulse, delayValue, fs); %PCA F1s
    outputMat(count,j+1) = f1;
    [f1, ~] = compareDirect(X_ica(:,j), directImpulse, delayValue, fs); %ICA F1s
    outputMat(count,j+ns+1) = f1;
end
outputMat(count,2*ns+2) = max(outputMat(count,2:5));
outputMat(count,2*ns+3) = max(outputMat(count,6:9));
count=count+1;
%---------------------------------------------------------------------------%

% !!! 3. cancelMecgCausal!!!
%---------------------------------------------------------------------------%
% mECG cancellation

[fecg1, ~, ~] = cancelMecgCausal(X_im(:,1), fs, templateSize,0,0);
[fecg2, ~, ~] = cancelMecgCausal(X_im(:,2), fs, templateSize,0,0);
[fecg3, ~, ~] = cancelMecgCausal(X_im(:,3), fs, templateSize,0,0);
[fecg4, ~, delayValue] = cancelMecgCausal(X_im(:,4), fs, templateSize,0,0);

% ICA/Post-processing
X_fecg = [fecg1, fecg2, fecg3, fecg4];
X_fecg = transpose(X_fecg);
r = ns;
[X_ica, ~, ~, ~] = fastICA(X_fecg,r);
X_ica = transpose(X_ica);
X_fecg = transpose(X_fecg);
% F1 Calculation
outputMat(count,1) = 3;
for j=1:ns
    [f1, ~] = compareDirect(X_fecg(:,j), directImpulse, delayValue, fs); %PCA F1s
    outputMat(count,j+1) = f1;
    [f1, ~] = compareDirect(X_ica(:,j), directImpulse, delayValue, fs); %ICA F1s
    outputMat(count,j+ns+1) = f1;
end
outputMat(count,2*ns+2) = max(outputMat(count,2:5));
outputMat(count,2*ns+3) = max(outputMat(count,6:9));
count=count+1;
%---------------------------------------------------------------------------%

% !!! 4. cancelMecgCausalSVD!!!
%---------------------------------------------------------------------------%
% mECG cancellation
[fecg1, ~, ~] = cancelMecgCausalSVD(X_im(:,1), fs, templateSize,0,0);
[fecg2, ~, ~] = cancelMecgCausalSVD(X_im(:,2), fs, templateSize,0,0);
[fecg3, ~, ~] = cancelMecgCausalSVD(X_im(:,3), fs, templateSize,0,0);
[fecg4, ~, delayValue] = cancelMecgCausalSVD(X_im(:,4), fs, templateSize,0,0);

% ICA/Post-processing
X_fecg = [fecg1, fecg2, fecg3, fecg4];
X_fecg = transpose(X_fecg);
r = ns;
[X_ica, ~, ~, ~] = fastICA(X_fecg,r);
X_ica = transpose(X_ica);
X_fecg = transpose(X_fecg);
% F1 Calculation
outputMat(count,1) = 4;
for j=1:ns
    [f1, ~] = compareDirect(X_fecg(:,j), directImpulse, delayValue, fs); %PCA F1s
    outputMat(count,j+1) = f1;
    [f1, ~] = compareDirect(X_ica(:,j), directImpulse, delayValue, fs); %ICA F1s
    outputMat(count,j+ns+1) = f1;
end
outputMat(count,2*ns+2) = max(outputMat(count,2:5));
outputMat(count,2*ns+3) = max(outputMat(count,6:9));
count=count+1;
%---------------------------------------------------------------------------%

% !!! 5. cancelMecgCausalAverage!!!
%---------------------------------------------------------------------------%
% mECG cancellation
[fecg1, ~, ~] = cancelMecgCausalAverage(X_im(:,1), fs, templateSize,0,0);
[fecg2, ~, ~] = cancelMecgCausalAverage(X_im(:,2), fs, templateSize,0,0);
[fecg3, ~, ~] = cancelMecgCausalAverage(X_im(:,3), fs, templateSize,0,0);
[fecg4, ~, delayValue] = cancelMecgCausalAverage(X_im(:,4), fs, templateSize,0,0);

% ICA/Post-processing
X_fecg = [fecg1, fecg2, fecg3, fecg4];
X_fecg = transpose(X_fecg);
r = ns;
[X_ica, ~, ~, ~] = fastICA(X_fecg,r);
X_ica = transpose(X_ica);
X_fecg = transpose(X_fecg);
% F1 Calculation
outputMat(count,1) = 5;
for j=1:ns
    [f1, ~] = compareDirect(X_fecg(:,j), directImpulse, delayValue, fs); %PCA F1s
    outputMat(count,j+1) = f1;
    [f1, ~] = compareDirect(X_ica(:,j), directImpulse, delayValue, fs); %ICA F1s
    outputMat(count,j+ns+1) = f1;
end
outputMat(count,2*ns+2) = max(outputMat(count,2:5));
outputMat(count,2*ns+3) = max(outputMat(count,6:9));
count=count+1;
%---------------------------------------------------------------------------%

% !!! 6. cancelMecgSmooth!!!
%---------------------------------------------------------------------------%
% mECG cancellation
windowTime = 0.049;
[fecg1, ~, ~] = cancelMecgTrial(X_im(:,1), fs, templateSize,windowTime,0,0);
[fecg2, ~, ~] = cancelMecgTrial(X_im(:,2), fs, templateSize,windowTime,0,0);
[fecg3, ~, ~] = cancelMecgTrial(X_im(:,3), fs, templateSize,windowTime,0,0);
[fecg4, ~, delayValue] = cancelMecgTrial(X_im(:,4), fs, templateSize,windowTime,0,0);

% ICA/Post-processing
X_fecg = [fecg1, fecg2, fecg3, fecg4];
X_fecg = transpose(X_fecg);
r = ns;
[X_ica, ~, ~, ~] = fastICA(X_fecg,r);
X_ica = transpose(X_ica);
X_fecg = transpose(X_fecg);
% F1 Calculation
outputMat(count,1) = 6;
for j=1:ns
    [f1, ~] = compareDirect(X_fecg(:,j), directImpulse, delayValue, fs); %PCA F1s
    outputMat(count,j+1) = f1;
    [f1, ~] = compareDirect(X_ica(:,j), directImpulse, delayValue, fs); %ICA F1s
    outputMat(count,j+ns+1) = f1;
end
outputMat(count,2*ns+2) = max(outputMat(count,2:5));
outputMat(count,2*ns+3) = max(outputMat(count,6:9));
count=count+1;
%---------------------------------------------------------------------------%

% !!! 7. cancelMecgAvgSmooth!!!
%---------------------------------------------------------------------------%
% mECG cancellation
windowTime = 0.029;
[fecg1, ~, ~] = cancelMecgAvgSmooth(X_im(:,1), fs, templateSize,windowTime,0,0);
[fecg2, ~, ~] = cancelMecgAvgSmooth(X_im(:,2), fs, templateSize,windowTime,0,0);
[fecg3, ~, ~] = cancelMecgAvgSmooth(X_im(:,3), fs, templateSize,windowTime,0,0);
[fecg4, ~, delayValue] = cancelMecgAvgSmooth(X_im(:,4), fs, templateSize,windowTime,0,0);

% ICA/Post-processing
X_fecg = [fecg1, fecg2, fecg3, fecg4];
X_fecg = transpose(X_fecg);
r = ns;
[X_ica, ~, ~, ~] = fastICA(X_fecg,r);
X_ica = transpose(X_ica);
X_fecg = transpose(X_fecg);
% F1 Calculation
outputMat(count,1) = 7;
for j=1:ns
    [f1, ~] = compareDirect(X_fecg(:,j), directImpulse, delayValue, fs); %PCA F1s
    outputMat(count,j+1) = f1;
    [f1, ~] = compareDirect(X_ica(:,j), directImpulse, delayValue, fs); %ICA F1s
    outputMat(count,j+ns+1) = f1;
end
outputMat(count,2*ns+2) = max(outputMat(count,2:5));
outputMat(count,2*ns+3) = max(outputMat(count,6:9));
count=count+1;
%---------------------------------------------------------------------------%

outputFileName = strcat('testSynthDB/Results/','sub_',subStr,'_SNR_',SNRStr,'dB_','Noise','.xlsx');
text = {'Method', 'PCA_1', 'PCA_2', 'PCA_3', 'PCA_4', 'ICA_1', 'ICA_2', 'ICA_3', 'ICA_4', 'MAX_PCA', 'MAX_ICA'};
writecell(text,outputFileName,'WriteMode','overwritesheet');
writematrix(outputMat, outputFileName, 'WriteMode','append');

sig1 = fecg1DIR+mecg1;
sig2 = fecg2DIR+mecg2;
sig3 = fecg3DIR+mecg3;
sig4 = fecg4DIR+mecg4;

X_raw = [sig1 sig2 sig3 sig4];
[ndt, ns] = size(X_raw);
arrayRows = 7; %SEVEN METHODS 
arrayCols = 2*ns + 3; % format: fecg 1-4, ica 1-4, + max PCA and max ICA
outputMat = zeros([arrayRows arrayCols]);
% Preprocessing
order = 633;
[X_bp, ~] = bandpassFilter(X_raw,BandPass2HzSamp250Hz,order,0); % Bandpass filtering, delay is not used throughout since signal is centred
X_im=FecgImpArtCanc(X_bp,fs); % Impulse artefact cancellation

%START LOOP
count=1;
templateSize = 20;

% !!! 1. cancelMecg !!!
%---------------------------------------------------------------------------%
% mECG cancellation
[fecg1, ~, ~] = cancelMecg(X_im(:,1), fs, templateSize,0,0);
[fecg2, ~, ~] = cancelMecg(X_im(:,2), fs, templateSize,0,0);
[fecg3, ~, ~] = cancelMecg(X_im(:,3), fs, templateSize,0,0);
[fecg4, ~, delayValue] = cancelMecg(X_im(:,4), fs, templateSize,0,0);

% ICA/Post-processing
X_fecg = [fecg1, fecg2, fecg3, fecg4];
X_fecg = transpose(X_fecg);
r = ns;
[X_ica, ~, ~, ~] = fastICA(X_fecg,r);
X_ica = transpose(X_ica);
X_fecg = transpose(X_fecg);
% 4. F1 Calculation
outputMat(count,1) = 1;
for j=1:ns
    [f1, ~] = compareDirect(X_fecg(:,j), directImpulse, delayValue, fs); %PCA F1s
    outputMat(count,j+1) = f1;
    [f1, ~] = compareDirect(X_ica(:,j), directImpulse, delayValue, fs); %ICA F1s
    outputMat(count,j+ns+1) = f1;
end
outputMat(count,2*ns+2) = max(outputMat(count,2:5));
outputMat(count,2*ns+3) = max(outputMat(count,6:9));
count=count+1;
%---------------------------------------------------------------------------%

% !!! 2. cancelMecgAverage !!!
%---------------------------------------------------------------------------%
% mECG cancellation

[fecg1, ~, ~] = cancelMecgAverage(X_im(:,1), fs, templateSize,0,0);
[fecg2, ~, ~] = cancelMecgAverage(X_im(:,2), fs, templateSize,0,0);
[fecg3, ~, ~] = cancelMecgAverage(X_im(:,3), fs, templateSize,0,0);
[fecg4, ~, delayValue] = cancelMecgAverage(X_im(:,4), fs, templateSize,0,0);

% ICA/Post-processing
X_fecg = [fecg1, fecg2, fecg3, fecg4];
X_fecg = transpose(X_fecg);
r = ns;
[X_ica, ~, ~, ~] = fastICA(X_fecg,r);
X_ica = transpose(X_ica);
X_fecg = transpose(X_fecg);
% F1 Calculation
outputMat(count,1) = 2;
for j=1:ns
    [f1, ~] = compareDirect(X_fecg(:,j), directImpulse, delayValue, fs); %PCA F1s
    outputMat(count,j+1) = f1;
    [f1, ~] = compareDirect(X_ica(:,j), directImpulse, delayValue, fs); %ICA F1s
    outputMat(count,j+ns+1) = f1;
end
outputMat(count,2*ns+2) = max(outputMat(count,2:5));
outputMat(count,2*ns+3) = max(outputMat(count,6:9));
count=count+1;
%---------------------------------------------------------------------------%

% !!! 3. cancelMecgCausal!!!
%---------------------------------------------------------------------------%
% mECG cancellation

[fecg1, ~, ~] = cancelMecgCausal(X_im(:,1), fs, templateSize,0,0);
[fecg2, ~, ~] = cancelMecgCausal(X_im(:,2), fs, templateSize,0,0);
[fecg3, ~, ~] = cancelMecgCausal(X_im(:,3), fs, templateSize,0,0);
[fecg4, ~, delayValue] = cancelMecgCausal(X_im(:,4), fs, templateSize,0,0);

% ICA/Post-processing
X_fecg = [fecg1, fecg2, fecg3, fecg4];
X_fecg = transpose(X_fecg);
r = ns;
[X_ica, ~, ~, ~] = fastICA(X_fecg,r);
X_ica = transpose(X_ica);
X_fecg = transpose(X_fecg);
% F1 Calculation
outputMat(count,1) = 3;
for j=1:ns
    [f1, ~] = compareDirect(X_fecg(:,j), directImpulse, delayValue, fs); %PCA F1s
    outputMat(count,j+1) = f1;
    [f1, ~] = compareDirect(X_ica(:,j), directImpulse, delayValue, fs); %ICA F1s
    outputMat(count,j+ns+1) = f1;
end
outputMat(count,2*ns+2) = max(outputMat(count,2:5));
outputMat(count,2*ns+3) = max(outputMat(count,6:9));
count=count+1;
%---------------------------------------------------------------------------%

% !!! 4. cancelMecgCausalSVD!!!
%---------------------------------------------------------------------------%
% mECG cancellation
[fecg1, ~, ~] = cancelMecgCausalSVD(X_im(:,1), fs, templateSize,0,0);
[fecg2, ~, ~] = cancelMecgCausalSVD(X_im(:,2), fs, templateSize,0,0);
[fecg3, ~, ~] = cancelMecgCausalSVD(X_im(:,3), fs, templateSize,0,0);
[fecg4, ~, delayValue] = cancelMecgCausalSVD(X_im(:,4), fs, templateSize,0,0);

% ICA/Post-processing
X_fecg = [fecg1, fecg2, fecg3, fecg4];
X_fecg = transpose(X_fecg);
r = ns;
[X_ica, ~, ~, ~] = fastICA(X_fecg,r);
X_ica = transpose(X_ica);
X_fecg = transpose(X_fecg);
% F1 Calculation
outputMat(count,1) = 4;
for j=1:ns
    [f1, ~] = compareDirect(X_fecg(:,j), directImpulse, delayValue, fs); %PCA F1s
    outputMat(count,j+1) = f1;
    [f1, ~] = compareDirect(X_ica(:,j), directImpulse, delayValue, fs); %ICA F1s
    outputMat(count,j+ns+1) = f1;
end
outputMat(count,2*ns+2) = max(outputMat(count,2:5));
outputMat(count,2*ns+3) = max(outputMat(count,6:9));
count=count+1;
%---------------------------------------------------------------------------%

% !!! 5. cancelMecgCausalAverage!!!
%---------------------------------------------------------------------------%
% mECG cancellation
[fecg1, ~, ~] = cancelMecgCausalAverage(X_im(:,1), fs, templateSize,0,0);
[fecg2, ~, ~] = cancelMecgCausalAverage(X_im(:,2), fs, templateSize,0,0);
[fecg3, ~, ~] = cancelMecgCausalAverage(X_im(:,3), fs, templateSize,0,0);
[fecg4, ~, delayValue] = cancelMecgCausalAverage(X_im(:,4), fs, templateSize,0,0);

% ICA/Post-processing
X_fecg = [fecg1, fecg2, fecg3, fecg4];
X_fecg = transpose(X_fecg);
r = ns;
[X_ica, ~, ~, ~] = fastICA(X_fecg,r);
X_ica = transpose(X_ica);
X_fecg = transpose(X_fecg);
% F1 Calculation
outputMat(count,1) = 5;
for j=1:ns
    [f1, ~] = compareDirect(X_fecg(:,j), directImpulse, delayValue, fs); %PCA F1s
    outputMat(count,j+1) = f1;
    [f1, ~] = compareDirect(X_ica(:,j), directImpulse, delayValue, fs); %ICA F1s
    outputMat(count,j+ns+1) = f1;
end
outputMat(count,2*ns+2) = max(outputMat(count,2:5));
outputMat(count,2*ns+3) = max(outputMat(count,6:9));
count=count+1;
%---------------------------------------------------------------------------%

% !!! 6. cancelMecgSmooth!!!
%---------------------------------------------------------------------------%
% mECG cancellation
windowTime = 0.049;
[fecg1, ~, ~] = cancelMecgTrial(X_im(:,1), fs, templateSize,windowTime,0,0);
[fecg2, ~, ~] = cancelMecgTrial(X_im(:,2), fs, templateSize,windowTime,0,0);
[fecg3, ~, ~] = cancelMecgTrial(X_im(:,3), fs, templateSize,windowTime,0,0);
[fecg4, ~, delayValue] = cancelMecgTrial(X_im(:,4), fs, templateSize,windowTime,0,0);

% ICA/Post-processing
X_fecg = [fecg1, fecg2, fecg3, fecg4];
X_fecg = transpose(X_fecg);
r = ns;
[X_ica, ~, ~, ~] = fastICA(X_fecg,r);
X_ica = transpose(X_ica);
X_fecg = transpose(X_fecg);
% F1 Calculation
outputMat(count,1) = 6;
for j=1:ns
    [f1, ~] = compareDirect(X_fecg(:,j), directImpulse, delayValue, fs); %PCA F1s
    outputMat(count,j+1) = f1;
    [f1, ~] = compareDirect(X_ica(:,j), directImpulse, delayValue, fs); %ICA F1s
    outputMat(count,j+ns+1) = f1;
end
outputMat(count,2*ns+2) = max(outputMat(count,2:5));
outputMat(count,2*ns+3) = max(outputMat(count,6:9));
count=count+1;
%---------------------------------------------------------------------------%

% !!! 7. cancelMecgAvgSmooth!!!
%---------------------------------------------------------------------------%
% mECG cancellation
windowTime = 0.029;
[fecg1, ~, ~] = cancelMecgAvgSmooth(X_im(:,1), fs, templateSize,windowTime,0,0);
[fecg2, ~, ~] = cancelMecgAvgSmooth(X_im(:,2), fs, templateSize,windowTime,0,0);
[fecg3, ~, ~] = cancelMecgAvgSmooth(X_im(:,3), fs, templateSize,windowTime,0,0);
[fecg4, ~, delayValue] = cancelMecgAvgSmooth(X_im(:,4), fs, templateSize,windowTime,0,0);

% ICA/Post-processing
X_fecg = [fecg1, fecg2, fecg3, fecg4];
X_fecg = transpose(X_fecg);
r = ns;
[X_ica, ~, ~, ~] = fastICA(X_fecg,r);
X_ica = transpose(X_ica);
X_fecg = transpose(X_fecg);
% F1 Calculation
outputMat(count,1) = 7;
for j=1:ns
    [f1, ~] = compareDirect(X_fecg(:,j), directImpulse, delayValue, fs); %PCA F1s
    outputMat(count,j+1) = f1;
    [f1, ~] = compareDirect(X_ica(:,j), directImpulse, delayValue, fs); %ICA F1s
    outputMat(count,j+ns+1) = f1;
end
outputMat(count,2*ns+2) = max(outputMat(count,2:5));
outputMat(count,2*ns+3) = max(outputMat(count,6:9));
count=count+1;
%---------------------------------------------------------------------------%

outputFileName = strcat('testSynthDB/Results/','sub_',subStr,'_SNR_',SNRStr,'dB_','Clean','.xlsx');
text = {'Method', 'PCA_1', 'PCA_2', 'PCA_3', 'PCA_4', 'ICA_1', 'ICA_2', 'ICA_3', 'ICA_4', 'MAX_PCA', 'MAX_ICA'};
writecell(text,outputFileName,'WriteMode','overwritesheet');
writematrix(outputMat, outputFileName, 'WriteMode','append');

end