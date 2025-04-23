function [fecgSignal, mecgSignal, delay] = cancelMecgHybrid(ecgSignal, fs, templateSize, windowTime, debugFlag, plotFlag)
% Syntax:       [fecgSignal, rPeakSignal] = cancelMecg(ecgSignal, templateSize)
%               
% Inputs:       ecgSignal: a one-dimensional ecg data array containing
%               both maternal and fetal components
%
%               fs: sampling rate of signal
%               
%               templateSize: number of detected maternal QRS-complexes
%               used to build the template before applying PCA.
%               Importantly, templateSize is expected to be an even
%               number.The function will result in an error if it is odd.
%               This is done so that the number of sub-signals either side
%               of an R-peak (before and after) is identical.
%
%               windowTime: length of interQRS averaging filter (in
%               seconds)
%
%               exclusionTime: length of QRS window to be excluded from
%               average filter
%
%               debugFlag: '0' by default, '1' simultaneously plots the
%               detected R-peaks above the corresponding ecg signal
%
%               plotFlag: '0' by default, '1' simultaneously plots the
%               original signal, the cancelled maternal and estimated fetal
%               ecg signals
%
%               
% Outputs:      fecgSignal: estimated fetal ecg signal after
%               cancelling the maternal components
%               
%               mecgSignal: estimated maternal signal using PCA
%
%               delay: increment at beginning of signal due to PCA
%               extension
%               
% Description:  Identical to cancelMecg, however incorporates an averaging
% filter to smoothen the estimated mECG.

if(nargin<5), debugFlag=0; plotFlag=0; end
if(nargin<6), plotFlag=0; end

if(mod(templateSize,2) ~= 0)
    error('Template size should be even.')
end

mili = fs/1000; %sub-signal duration is 601 ms in total
samplesBefore= ceil(200 * mili);
samplesAfter = ceil(400 * mili);
samplesPer = samplesBefore + samplesAfter + 1;

n = length(ecgSignal);

% Extend the signal to avoid any index errors due to location of R-peaks
% near ends of ecg input array

extendedSignal = zeros([1 n+samplesPer]);
for i=1:length(ecgSignal)
    extendedSignal(i+samplesBefore)=ecgSignal(i); 
end

%% Initialisations
n = length(extendedSignal);
[qrs_amp_raw,qrs_i_raw,delay] = pan_tompkin(ecgSignal,fs,0); % find R-peak locations
impulseTrain =  zeros([1 n]);
mecgSignal = zeros([1 n]);
Rpeaks = length(qrs_i_raw);
resultMatrix = zeros([Rpeaks samplesPer]);
resultMatrixAvg = zeros([Rpeaks samplesPer]);
dataMatrix = zeros([(templateSize+1) samplesPer]);
averageSignal = zeros([1 samplesPer]);

%% Smoothing parameters
QRinterval = 0.030; % QRS interval begins approx. 30 ms before beginning of R-wave
RSinterval = 0.090; % lasts in total ~100 ms so add 70 ms
criticalTime = windowTime; % apply linear filter as QRS complex begins (eliminates jitter).

smoothedMecgSignal = zeros([1 n]);
QRsamples = ceil(QRinterval*fs);
RSsamples = ceil(RSinterval*fs);
windowSamples = ceil(windowTime*fs);
windowSize = floor(windowSamples/2);
criticalSamples = windowSamples;
criticalSize = windowSize;
smoothingMarker = ones([1 n]); %value 2 if criticalLength linear filter is applied, 1 if window smoothing


for i=1:Rpeaks
   impulseTrain(qrs_i_raw(i)) = 1; % indicate locations of R-wave beginning
end

if(debugFlag) % Plot impulse train of detected R-peaks vs. raw signal
    plotTime =(0:length(ecgSignal)-1)*1/fs;
    plot(plotTime,50*impulseTrain(1:length(ecgSignal)),LineWidth=1, Color='red')
    hold on
    plot(plotTime,ecgSignal,Color='blue');
    hold off
    legend({'Impulse Train', 'Direct_1'});
end

if Rpeaks <= templateSize % No. of detected R-peaks is less than the no. of maternal templates used within the SVD data matrix
    error('Template size too large for detected maternal R-peaks OR maternal R-peaks could not be adequately detected')
end

subComplexMatrix = zeros([Rpeaks samplesPer]); %create matrix of sub-complexes (QRS sub-signals)
for i=1:Rpeaks
    for j=1:samplesPer
        subComplexMatrix(i,j) = extendedSignal(qrs_i_raw(i) + j);
    end
end

% Procedure is now divided into three parts.
%
% I. PCA for the first templateSize/2 sub-complexes. This is done by simply
% using the first templateSize+1 sub-complexes, since they are the first
% entries.

for i=1:templateSize+1
    for j=1:samplesPer
        dataMatrix(i,j) = subComplexMatrix(i,j);
    end
    averageSignal = averageSignal + subComplexMatrix(i,:);
end

% Apply PCA and estimate first templateSize/2 mECG sub-complexes using first 
% 3 principal components
[Zpca, U, mu, eigVecs] = PCA(dataMatrix,3);
Zreduced = U * Zpca + repmat(mu,1,samplesPer);

averageSignal = averageSignal/(templateSize+1);
for i=1:templateSize/2
    for j=1:samplesPer
        resultMatrix(i,j) = Zreduced(i,j);
    end
    resultMatrixAvg(i,:) = averageSignal;
end
averageSignal = zeros([1 samplesPer]);

% II. PCA for the (templateSize/2, Rpeaks-templateSize/2) sub-complexes
for k=templateSize/2 + 1:Rpeaks-templateSize/2

    for i=1:templateSize+1
        for j=1:samplesPer
            dataMatrix(i,j) = subComplexMatrix(k - templateSize/2 - 1 + i, j);
        end
        averageSignal = averageSignal + subComplexMatrix(k - templateSize/2 - 1 + i,:);
    end
    
    % PCA for constructed data matrix 
    [Zpca, U, mu, eigVecs] = PCA(dataMatrix,3);
    Zreduced = U * Zpca + repmat(mu,1,samplesPer);

    resultMatrix(k,:) = Zreduced(templateSize/2+1,:);
    
    averageSignal = averageSignal/(templateSize+1);
    resultMatrixAvg(k,:) = averageSignal;
    averageSignal = zeros([1 samplesPer]);
end

% III. PCA for the last templateSize/2 sub-complexes

for i=1:templateSize+1
    for j=1:samplesPer
        dataMatrix(i,j) = subComplexMatrix(Rpeaks-templateSize-1+i,j);
    end
    averageSignal = averageSignal + subComplexMatrix(Rpeaks-templateSize-1+i,:);
end

% Apply PCA and estimate first templateSize/2 mECG sub-complexes using first 
% 3 principal components
[Zpca, U, mu, eigVecs] = PCA(dataMatrix,3);
Zreduced = U * Zpca + repmat(mu,1,samplesPer);
averageSignal = averageSignal/(templateSize+1);

for i=Rpeaks-templateSize/2+1:Rpeaks
    resultMatrix(i,:) = Zreduced(templateSize + 1 - Rpeaks + i,:);
    resultMatrixAvg(i,:) = averageSignal;
end
%A=extendedSignal(qrs_i_raw(1):qrs_i_raw(1)+samplesPer)
%B=resultMatrixAvg(i,:)
countAVG=0;
countPCA=0;
for i=1:Rpeaks
    residualAVG = extendedSignal(qrs_i_raw(i)+1:qrs_i_raw(i)+samplesPer) - resultMatrixAvg(i,:);
    residualPCA = extendedSignal(qrs_i_raw(i)+1:qrs_i_raw(i)+samplesPer) - resultMatrix(i,:);
    %if(sampen(residualAVG, 2, 0.2,'chebychev') < sampen(residualPCA, 2, 0.2,'chebychev') && mean(residualPCA.^2)<mean(residualAVG.^2))
    if(mean(residualPCA.^2)<mean(residualAVG.^2))
        for j=1:samplesPer
            mecgSignal(qrs_i_raw(i)+j)= resultMatrixAvg(i,j);
        end
        countAVG=countAVG+1;   
    else
        for j=1:samplesPer
            mecgSignal(qrs_i_raw(i)+j)= resultMatrix(i,j);
        end
        countPCA=countPCA+1;
    end
    locs = qrs_i_raw(i)+samplesBefore;  % index of beginning of R wave
    if(i==1) % avoid edge case issues by not smoothing signal before first maternal beat
        smoothingMarker((locs+RSsamples):(locs+RSsamples+criticalSamples)) = 2;
        smoothingMarker((locs-QRsamples):(locs+RSsamples)) = 0;
    elseif(i==Rpeaks) % avoid edge case issues by not smoothing signal after last maternal beat
        smoothingMarker((locs-QRsamples-criticalSamples):(locs-QRsamples)) = 2;
        smoothingMarker((locs-QRsamples):(locs+RSsamples)) = 0;
    else
        smoothingMarker((locs-QRsamples-criticalSamples):(locs-QRsamples)) = 2;
        smoothingMarker((locs+RSsamples):(locs+RSsamples+criticalSamples)) = 2;
        smoothingMarker((locs-QRsamples):(locs+RSsamples)) = 0;
    end
    
end
countPCA = countPCA
countAVG = countAVG
% Do not smooth beginning and end of signal
smoothedMecgSignal = mecgSignal;
smoothingMarker(1:(qrs_i_raw(1)+samplesBefore)) = 0;
smoothingMarker((qrs_i_raw(Rpeaks)+samplesBefore):length(smoothingMarker)) = 0;
%smoothingMarker = smoothingMarker;

varCritCount=criticalSamples;
varCritSize =criticalSize;

for i=1:n
    if(smoothingMarker(i)==0)
        smoothedMecgSignal(i)=mecgSignal(i);
        varCritCount=criticalSamples;
        varCritSize =criticalSize;
    elseif(smoothingMarker(i)==2) %if equals to 2 
        varCritCount = varCritCount-1;
        if(varCritCount<2)
            continue;
        end
        sum = 0;
        if(rem(varCritCount,2)==0)
            for j=(-varCritSize):(varCritSize-1)
                sum = sum + mecgSignal(i+j);
            end
        else
            varCritSize = varCritSize-1;
            for j=-varCritSize:varCritSize
                sum = sum + mecgSignal(i+j);
            end
        end
        average = sum/(varCritCount);
        smoothedMecgSignal(i)=average;
    else
        sum = 0;
        for j=-windowSize:windowSize
            sum = sum + mecgSignal(i+j);
        end
        average = sum/(2*windowSize+1);
        smoothedMecgSignal(i)=average;
        varCritCount=criticalSamples;
        varCritSize =criticalSize;
    end
end    

fecgSignal = extendedSignal - mecgSignal;
%fecgSignal = extendedSignal - smoothedMecgSignal;
fecgSignal = transpose(fecgSignal);
delay = samplesBefore;

if(plotFlag)
    plotTime =(0:length(extendedSignal)-1)*1/fs;
    subplot(3,1,1)
    plot(plotTime,extendedSignal,Color='blue')
    title('Original ECG');
    subplot(3,1,2)
    plot(plotTime,mecgSignal,Color=	'green')
    hold on;
    plot(plotTime,smoothedMecgSignal,Color='red')
    legend('PCA maternal', 'Smoothed maternal')
    hold off;
    title('Maternal ECG');
    subplot(3,1,3)
    plot(plotTime,fecgSignal,Color='cyan')
    title('Fetal ECG')
end



