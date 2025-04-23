function [X_raw,directSignal] = edfECGRead(inputString,plotFlag)
% Syntax:       [X_raw,directSignal] = edfECGRead(inputString,plotFlag)
%               
% Inputs:       inputString: location and file name of edf file containing 
%               ecg signal data 
%               
%               plotFlag: '0' by default, '1' plots the ecg signals as read
%               from the edf file
%               
% Outputs:      X_raw: matrix of unprocessed ecg signals, where each column
%               vector is an ecg signal
%               
%               directSignal: scalp measurement of associated recording
%               
% Description:  Store the unprocessed filters in the X_raw matrix, whose
% column vectors represent each signal

if(nargin<2), plotFlag=0; end

data = edfread(inputString);

raw1 = data.Abdomen_1{1};
raw2 = data.Abdomen_2{1};
raw3 = data.Abdomen_3{1};
raw4 = data.Abdomen_4{1};
direct1 = data.Direct_1{1};

for i = 2:60
    raw1 = cat(1,raw1,data.Abdomen_1{i});
    raw2 = cat(1,raw2,data.Abdomen_2{i});
    raw3 = cat(1,raw3,data.Abdomen_3{i});
    raw4 = cat(1,raw4,data.Abdomen_4{i});
    direct1 = cat(1,direct1,data.Direct_1{i});
end

X_raw = [raw1 raw2 raw3 raw4];
directSignal = direct1;

if(plotFlag)
    % Display raw data
    subplot(4,1,1);
    plot(raw1)
    title('Raw 1')
    
    subplot(4,1,2); 
    plot(raw2)
    title('Raw 2')  
    
    subplot(4,1,3); 
    plot(raw3)
    title('Raw 3')
    
    subplot(4,1,4); 
    plot(raw4)
    title('Raw 4')
end
end