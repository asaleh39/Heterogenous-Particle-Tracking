a = 2; % Number of Particles

b = 300; % Steps

dt = 0.02;

pixel = 512;

delta = [0.503519500255143;0.663660403841073]; % Delta Constant


%% Load Sim Data
partdatanoshift = readtable('bmsimdata300-2noshift.xlsx');
endpart = partdatanoshift.particle(end);

h = 0;
k = find(partdatanoshift.particle == h);
x = partdatanoshift.x(k);
y = partdatanoshift.y(k);

for h = 1:endpart
    k = find(partdatanoshift.particle == h);
    
    x = [x,partdatanoshift.x(k)]';
    y = [y,partdatanoshift.y(k)]';
end


partdatashift = readtable('bmsimdata300-2shift.xlsx');
endpart = partdatashift.particle(end);

h = 0;
k = find(partdatashift.particle == h);
xgood = partdatashift.x(k)';
ygood = partdatashift.y(k)';

for h = 1:1
    k = find(partdatashift.particle == h);
    
    xgood = [xgood;partdatashift.x(k)'];
    ygood = [ygood;partdatashift.y(k)'];
end

%% Use x and y for shift, xgood ygood for shift
% Plots trajectory for particles without and with shift

for h = 1:a
    plot(x(h,1:b),y(h,1:b))
    axis([0,pixel,0,pixel])
    drawnow
    hold on
end
hold off
title('2D Brownian Motion')
xlabel('x')
ylabel('y')


figure

for h = 1:a
    plot(x(h,1:b),y(h,1:b),xgood(h,1:b),ygood(h,1:b))
    axis([0,pixel,0,pixel])
    drawnow
    hold on
end

hold off
title('2D Brownian Motion with Shift')
xlabel('x')
ylabel('y')



%% MSD
% Calculates the MSD | Needs to be checked
xdisp = 0;
ydisp = 0;
for j = 1:length(xgood)-1
    xdisp = 0;
    ydisp = 0;
    
    for k = 1:length(xgood)-j
        xdisp(k) = xgood(1,k+j)-xgood(1,k);
        ydisp(k) = ygood(1,k+j)-ygood(1,k);
        
    end
    msd(j) = mean(xdisp)^2 + mean(ydisp)^2;
end

% Please ignore the code below I was just testing different things and
% plotting them to see if any changes worked

% loglog(1:b-1,msd);
% axis([0 100 0 1000])
% hold on
% loglog(1:b,MSD_1_1D);
% axis([0 100 0 1000])
% hold off
% 
% xdisp = 0;
% ydisp = 0;
% 
% for k = 1:length(xgood)-1
%     xdisp(k) = xgood(1,k+1)-xgood(1,k);
%     ydisp(k) = ygood(1,k+1)-ygood(1,k);
% end
% xydispsqr = 0;
% for k = 1:length(xdisp)
%     xydispsqr(k) = xdisp(k)^2 + ydisp(k)^2;
% end
% 
% msd2 = [sqrt(mean(xdisp)^2 + mean(ydisp)^2)];
% 


%% Haar Wavelets
% Taken from Supplementary Information From 
% Disentangling Random Motion and Flow in a Complex Medium

%%%%%%%%A^n_k Haar%%%%%%%%%%

ak_k_less_eq_n = @(n,k) (1/(60*n^2*(n+1)^2))*(3*k^5-k^4*(10*n+5)...
    +k^3*(40*n^2+40*n-5)-k^2*(80*n^3+120*n^2+30*n-5)...
    +k*(60*n^4+120*n^3+80*n^2+20*n+2));

ak_n_k_2n = @(n,k) (1/(60*n^2*(n+1)^2))*(-k^5+k^4*(10*n+5)-k^3*(40*n^2+40*n+5)...
    +k^2*(80*n^3+120*n^2+30*n-5)-k*(80*n^4+160*n^3+60*n^2-20*n-6)...
    +(44*n^5+110*n^4+80*n^3+10*n^2-4*n));

ak_k_great_eq_2n =  @(n,k)((2*n+1)*(3*n^2+3*n+4))/(30*n*(n+1));

%%%%%%%%B^n_k Haar%%%%%%%%%%

bk_k_less_eq_n = @(n,k) 1/(n^2*(n+1)^2)*(-k^3+k^2*(2*n+1)-4*k*n*(n+1)+2*n*(n+1)...
    *(n^2+n+1));

bk_n_k_2n = @(n,k) 1/(3*n^2*(n+1)^2)*(k^3-3*k^2*(2*n+1)+2*k*(6*n^2+6*n+1)...
    +2*n^3*(3*n+4)-12*n^2*(n+1)-2*n);

bk_k_great_eq_2n = @(n,k) (6*n^2-2*n+2)/(3*n*(n+1));

% MSDn = delta(1)*dt*akn + 4 * 1^2 * bkn;
% 
% MSD_Rescaledn = MSDn/bkn;



% for k23 = 1:b
%     MSD_1_1D(k23) = 4*delta(1)*k23;
% end

% figure
% loglog(1:b-1,msd);
% axis([0 100 0 1000])

% nidx = 1:20;
akn = zeros(20,b-1); % b-1 is used to create arrays that are inline with the dimensions of the MSD which loses 1 data value
for n = 1:20
    for k = 1:b-1
        
        if k<=n
            akn(n,k) = ak_k_less_eq_n(n,k);
            
        elseif n<k && k<2*n
            akn(n,k) = ak_n_k_2n(n,k);
            
        elseif k>=(2*n)
            akn(n,k) = ak_k_great_eq_2n(n,k);
        end
    end
end

bkn = zeros(20,b-1);
for n = 1:20
    for k = 1:b-1
        
        if k<=n
            bkn(n,k) = bk_k_less_eq_n(n,k);
            
        elseif n<k && k<2*n
            bkn(n,k) = bk_n_k_2n(n,k);
            
        elseif k>=(2*n)
            bkn(n,k) = bk_k_great_eq_2n(n,k);
        end
    end
end

% mas_noadj = (4*0.7*dt).*akn;
% 
% msd_adj = mas_noadj./bkn;

mas_noadj = msd.*akn; % MSD of Adjusted trajectory that is not rescaled

msd_adj = mas_noadj./bkn; % MSD of adjusted trajectory that has been rescaled

tk = akn./bkn; % Rescaled time

figure
loglog(1:b-1,mas_noadj(4,:),...
    1:b-1,mas_noadj(8,:),...
    1:b-1,mas_noadj(12,:),...
    1:b-1,mas_noadj(16,:),...
    1:b-1,mas_noadj(20,:));
axis([0 100 0 5])
% 4, 8, 12, 16, and 20 are wavelet spans of the Haar Wavelet
% We want the MSD that hasnt been rescaled (mas_noadj)to increase for 
% short times, but flatten out for times >2n.




% Ignore the code below I was testing some variables
% figure
% loglog(1:b-2,MSDn(4,:),...
%     1:b-2,MSDn(8,:),...
%     1:b-2,MSDn(12,:),...
%     1:b-2,MSDn(16,:),...
%     1:b-2,mas_noadj(20,:));
% axis([0 100 0 5])

% 
% figure
% plot(tk(2,:),msd_adj(2,:))

