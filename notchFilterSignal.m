function [X_nf, delay] = notchFilterSignal(X_in,filterType,order)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

zeroPadding = order; %added to beginning and end
ndt=length(X_in);

% Accounting for phase
X_temp = zeros( [1 (1 + zeroPadding*2 + ndt)] );
X_nf = zeros( [1 ndt] );

for i=1:ndt
    X_temp(i+zeroPadding+1)=X_in(i);
end
delay=zeroPadding;
X_f = filter(filterType,X_temp);

for i=1:ndt
    X_nf(i)=X_f(i+3*delay/2);
end

end