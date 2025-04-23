function templateSizeAverageVarRecord(recordName, lowerBound, upperBound)

%templateSizeVarRecord: Outputs the F1 performance characteristics of the
%record for variations in the QRS complex template size to be used in SVD.
%The output is a 9-d array of each fecg and ica F1 value for each
%templateSize value in the cancelMecg method.

%   Detailed explanation goes here

% 0. Read Data
% Read abdomenal signals
strRecordLocation = strcat('/mcode/databases/', recordName); % Records of set A must be in this directory
[sig1, ~, ~] = rdsamp(strRecordLocation, 1);
[sig2, ~, ~] = rdsamp(strRecordLocation, 2);
[sig3, ~, ~] = rdsamp(strRecordLocation, 3);
[sig4, fs, time] = rdsamp(strRecordLocation, 4);

% Read annotated fetal QRS locations
strFileLocation = strcat('setAfqrs\',recordName,'.fqrs.txt');
fileID = fopen(strFileLocation,'r');
directImpulse = fscanf(fileID,'%d');
fclose(fileID);

X_raw = [sig1 sig2 sig3 sig4];
[ndt, ns] = size(X_raw);
arrayRows = 1+(upperBound-lowerBound)/2; 
arrayCols = 2*ns + 1; % format: fecg 1-4, ica 1-4,
outputMat = zeros([arrayRows arrayCols]);
% 1. Preprocessing
order = 2530;
[X_bp, ~] = bandpassFilter(X_raw,BandPassCutoff3Hz,order,0); % Bandpass filtering, delay is not used throughout since signal is centred
X_nf = FecgNotchFilt(X_bp,fs);% Notch filtering
X_im=FecgImpArtCanc(X_nf,fs); % Impulse artefact cancellation


count=1;
% !!! Template size variation !!!
for i=lowerBound:2:upperBound
    templateSize = i;
    % 2. mECG cancellation
    [fecg1, ~, ~] = cancelMecgAverage(X_im(:,1), fs, templateSize,0,0);
    [fecg2, ~, ~] = cancelMecgAverage(X_im(:,2), fs, templateSize,0,0);
    [fecg3, ~, ~] = cancelMecgAverage(X_im(:,3), fs, templateSize,0,0);
    [fecg4, ~, delayValue] = cancelMecgAverage(X_im(:,4), fs, templateSize,0,0);

    % 3. ICA/Post-processing
    X_fecg = [fecg1, fecg2, fecg3, fecg4];
    X_fecg = transpose(X_fecg);
    r = ns;
    [X_ica, ~, ~, ~] = fastICA(X_fecg,r);
    X_ica = transpose(X_ica);
    X_fecg = transpose(X_fecg);
    % 4. F1 Calculation
    outputMat(count,1) = i;
    for j=1:ns
        [f1, ~] = compareDirect(X_fecg(:,j), directImpulse, delayValue, fs); %PCA F1s
        outputMat(count,j+1) = f1;
        [f1, ~] = compareDirect(X_ica(:,j), directImpulse, delayValue, fs); %ICA F1s
        outputMat(count,j+ns+1) = f1;
    end
    count=count+1;
end

outputFileName = strcat(recordName,'Average','.xlsx');
text = {'Template Size', 'PCA_1', 'PCA_2', 'PCA_3', 'PCA_4', 'ICA_1', 'ICA_2', 'ICA_3', 'ICA_4'};
writecell(text,outputFileName,'WriteMode','overwritesheet');
writematrix(outputMat, outputFileName, 'WriteMode','append');

end