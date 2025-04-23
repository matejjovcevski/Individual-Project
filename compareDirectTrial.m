function [f1, estimatedR] = compareDirect(estimfECG, directImpulse, delayValue, fs)
% Syntax:       [f1] = compare(estimfECG,directfECG)
%               
% Inputs:       estimfECG: one-dimensional ecg data array containing
%               the estimated fECG signal
%               
%               directfECG: fECG signal measured directly (reference)
%               
%               fs: sampling rate
%               
%               delayValue: amount of delay added to direct impulse signal
%               (due to extension in PCA process in cancelMecg function)
%               
% Outputs:      F1: performance metric for the ability to detect a
%               fetal heart beat within a 100 sampe (ms for 1 kHz fs) 
%               window, it is expressed as a value between 0 and 1.
%
%               estimatedR: array of indeces of R peaks for the estimated
%               fECG
%
%               directR: array of indeces of R peaks for the reference
%               fECG
%
%               
% Description:  Uses the Pan and Tompkins R-peak detection algorithm to
% find R-peaks of reference(scalp) measurement and estimated fECG, and then
% calculates the F1 value for the estimated fECG.
% 

%if length(estimfECG) ~= length(directfECG)
%    error("Array sizes do not match")
%end

n = size(estimfECG);
milis = floor(0.050 *fs); % 50 miliseconds interval before and after

dir_impulseTrain =  zeros([1 n]);
for i=1:length(directImpulse)
   dir_impulseTrain(directImpulse(i)+delayValue) = 1;
end

[~,estim_qrs_i_raw,~] = pan_tompkin(estimfECG,fs,0);
estim_impulseTrain =  zeros([1 n]);
for i=1:length(estim_qrs_i_raw)
   estim_impulseTrain(estim_qrs_i_raw(i)) = 1;
end

tp = 0; % True positive
fp = 0; % False positive
fn = 0; % False negative

for i=1:length(estim_impulseTrain)
    if estim_impulseTrain(i) == 1
        for j=0:2*milis
            if j==2*milis
                fp = fp + 1;
                break;
            end
            if(i-milis+j<=0 || i-milis+j>=length(dir_impulseTrain))
                continue;
            end
            s=i-milis+j;
            if dir_impulseTrain(i - milis + j) == 1    
                tp = tp + 1;
                dir_impulseTrain(i - milis + j) = 2; %mark corresponding point to avoid double counting
                break;
            end
        end
    end
end

for i=1:length(dir_impulseTrain)
    if dir_impulseTrain(i) == 1
        fn = fn+1; %count remaining points that weren't detected
    end
end

f1 = 2 * tp / (2*tp + fp + fn);
estimatedR = estim_impulseTrain;

end