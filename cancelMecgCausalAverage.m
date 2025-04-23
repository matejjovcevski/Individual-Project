function [fecgSignal, mecgSignal, delay] = cancelMecgCausalAverage(ecgSignal, fs, templateSize, debugFlag, plotFlag)
% Syntax:       [fecgSignal, rPeakSignal] = cancelMecg(ecgSignal, templateSize)
%               
% Inputs:       ecgSignal: a one-dimensional ecg data array containing
%               both maternal and fetal components
%               
%               templateSize: number of detected maternal QRS-complexes
%               used to build the template before applying PCA.
%               Importantly, templateSize is expected to be an even
%               number.The function will result in an error if it is odd.
%               This is done so that the number of sub-signals either side
%               of an R-peak (before and after) is identical.
%               
%               fs: sampling rate of signal
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
% Description:  Uses the Pan and Tompkins R-peak detection algorithm to
% find the locations of the maternal R-peaks, and subsequently create
% sub-signals of duration 200 ms before and 300 ms after the R-peak. Thus,
% each sub-signal represents the maternal QRS complex for each heartbeat.
% PCA is then applied to a total of templateSize + 1 QRS complexes (except 
% for initial and final templateSize heartbeats) around a single R-peak, to 
% get a good estimate for the maternal QRS complex for that particular 
% R-peak. Only the first three principal components are chosen as they 
% contain 90 % of the signal (the maternal component is considered to be 
% much greater than the fetal). All of these estimated maternal sub-signals
% are then concatenated, with zeros padded between them, to fit the length
% of the initial ecgSignal input.

if(nargin<4), debugFlag=0; plotFlag=0; end
if(nargin<5), plotFlag=0; end

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

n = length(extendedSignal);
[qrs_amp_raw,qrs_i_raw,delay] = pan_tompkin(ecgSignal,fs,0); % find R-peak locations
impulseTrain =  zeros([1 n]);
Rpeaks = length(qrs_i_raw);
for i=1:Rpeaks
   impulseTrain(qrs_i_raw(i)) = 1; % indicate locations of R-peaks
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

resultMatrix = zeros([Rpeaks samplesPer]);
averageSignal = zeros([1 samplesPer]);

% Procedure is now divided into two parts.
%
% I. Avg for the first templateSize complexes. This is done by simply
% using the first templateSize sub-complexes, since they are the first
% entries.

for i=1:templateSize+1
    averageSignal = averageSignal + subComplexMatrix(i,:);
end
averageSignal = averageSignal/(templateSize+1);
for i=1:templateSize
    resultMatrix(i,:) = averageSignal;
end
averageSignal = zeros([1 samplesPer]);

% II. Avg for the rest of the sub-complexes
for k=(templateSize + 1):Rpeaks

    for i=1:templateSize+1
        averageSignal = averageSignal + subComplexMatrix(k - templateSize - 1 + i,:);
    end

    averageSignal = averageSignal/(templateSize+1);
    resultMatrix(k,:) = averageSignal;
    averageSignal = zeros([1 samplesPer]);

end

% Estimated maternal signal
mecgSignal = zeros([1 n]);
for i=1:Rpeaks
    for j=1:samplesPer
        mecgSignal(qrs_i_raw(i)+j)= resultMatrix(i,j);
    end
end

fecgSignal = extendedSignal - mecgSignal;
fecgSignal = transpose(fecgSignal);
delay = samplesBefore;

if(plotFlag)
    plotTime =(0:length(extendedSignal)-1)*1/fs;
    subplot(3,1,1)
    plot(plotTime,extendedSignal,Color='blue')
    title('Original ECG');
    subplot(3,1,2)
    plot(plotTime,mecgSignal,Color=	'green')
    title('Maternal ECG');
    subplot(3,1,3)
    plot(plotTime,fecgSignal,Color='red')
    title('Fetal ECG')
end



