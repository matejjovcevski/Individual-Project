function [smoothedSignal] = smoothMecg(inputSignal,fs,windowTime,exclusionTime,plotFlag)

% Applies an averagin filter of size windowSize to inputSignal

if(nargin<6); plotFlag=0; end
windowSize = ceil(windowTime*fs/2); % divide by two since averages values before and after current sample
exclusionWidth = ceil(exclusionTime*fs);
criticalLength = ceil(0.010 *fs/2); % use 11ms averaging filter when points near exclusionWidth are reached
smoothedSignal = inputSignal;
medianFilterLength = 11; % use 11 points when applying median filter to find max and min points of QRS complex


medianSignal = medfilt1(inputSignal,medianFilterLength);
[RpeakHI,indexHI] = max(medianSignal);
[RpeakLO,indexLO] = min(medianSignal);

zeroCrossing = floor((indexHI + indexLO)/2);
lowerBound = zeroCrossing - exclusionWidth - criticalLength;
upperBound = zeroCrossing + exclusionWidth + criticalLength;

if(lowerBound<1)
    lowerBound=1;
end
if(upperBound>length(inputSignal)-1)
    upperBound=length(inputSignal);
end

for i=1:lowerBound
    if(i<=(windowSize+1))
        sum = 0;
        for j=1:i
            sum = sum + inputSignal(j);
        end
        smoothedSignal(i)=sum/(2*windowSize+1);
    else
        sum = 0;
        for j=(i-windowSize):(i+windowSize)
            sum = sum + inputSignal(j);
        end
        smoothedSignal(i)=sum/(2*windowSize+1);
    end
end

for i=(lowerBound+1):(lowerBound+criticalLength)
    sum=0;
    for j=(i-criticalLength):(i+criticalLength)
        if(j<1 || j > length(inputSignal))
            continue;
        end
        sum = sum + inputSignal(j);
    end
    smoothedSignal(i) = sum/(2*criticalLength+1);
end

for i=upperBound:length(inputSignal)
    if(i>length(inputSignal)-windowSize-1)
        sum = 0;
        for j=i:length(inputSignal)
            sum = sum + inputSignal(j);
        end
        smoothedSignal(i)=sum/(2*windowSize+1);
    else
        sum = 0;
        for j=(i-windowSize):(i+windowSize)
            sum = sum + inputSignal(j);
        end
        smoothedSignal(i)=sum/(2*windowSize+1);
    end
end

for i=(upperBound+1):(upperBound+criticalLength)
    sum=0;
    for j=(i-criticalLength):(i+criticalLength)
        if(j<1 || j > length(inputSignal))
            continue;
        end
        sum = sum + inputSignal(j);
    end
    smoothedSignal(i) = sum/(2*criticalLength+1);
end

if(plotFlag)
    plotTime =(0:length(inputSignal)-1)*1/fs;
    subplot(2,1,1)
    plot(plotTime,inputSignal,Color='blue')
    title('Original ECG');
    subplot(2,1,2)
    plot(plotTime,smoothedSignal,Color=	'green')
    title('Smoothed ECG');
    hold on;
    plot(upperBound/fs,smoothedSignal(upperBound),Marker="*",Color='r');
    plot(lowerBound/fs,smoothedSignal(lowerBound),Marker="*",Color='r');
    plot(zeroCrossing/fs,smoothedSignal(zeroCrossing),Marker='+',Color='r');
    hold off;
end

end