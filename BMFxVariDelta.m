function [x_f] = BMFxVariDelta(m,n,Delta,DT,pixel)
%BMFX Summary of this function goes here
%   Detailed explanation goes here
delta = Delta;
dt = DT;
x = randi(pixel,m,1);
r = zeros(m,n);

for k = 1:m
    r(k,:) = normrnd(0,sqrt(delta(k)*2*dt),[1,n]);    
end

x_f = x(:,1) + cumsum(r')';
end

