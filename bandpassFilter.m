function [X_bp,delay] = bandpassFilter(X_in,filterType,order,plotFlag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if(nargin<4), plotFlag=0; end
zeroPadding = order;%added to beginning and end
[ndt, ns]=size(X_in);

% Delayed response
%{
for i=1:ns
    X_bp(:,i) = filter(filterType,X_in(:,i));
end
delay=0;
subplot(4,1,1);
plot(X_bp(:,1))
title('Imp 1')
    
subplot(4,1,2); 
plot(X_bp(:,2))
title('Imp 2')  

subplot(4,1,3); 
plot(X_bp(:,3))
title('Imp 3')
   
subplot(4,1,4); 
plot(X_bp(:,4))
title('Imp 4')
%}

% Accounting for phase
for i=1:ns
    X_temp(:,i) = zeros([1 ( 1 + zeroPadding*2 + length(X_in(:,i)) )]);
end

for i=1:ns
    for j=1:ndt
        X_temp(j+zeroPadding+1,i)=X_in(j,i);
    end
end

delay=zeroPadding;
for i=1:ns
    X_f(:,i) = filter(filterType,X_temp(:,i));
end

for i=1:ns
    for j=1:ndt
        X_bp(j,i)=X_f(j+3*ceil(delay/2) ,i);
    end
end    
if(plotFlag)
    % Display raw data
    subplot(4,1,1);
    plot(X_bp(:,1))
    title('Bp 1')
    
    subplot(4,1,2); 
    plot(X_bp(:,2))
    title('Bp 2')  
    
    subplot(4,1,3); 
    plot(X_bp(:,3))
    title('Bp 3')
    
    subplot(4,1,4); 
    plot(X_bp(:,4))
    title('Bp 4')
end
%}
end