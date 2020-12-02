close all

clear

delete bmsimdata.xlsx

clc

%% Declarations

a = 2; % Number of Particles

b = 300; % Steps

dt = 0.02;

pixel = 512;

delta = rand(a,1).^1; % Must be in dimension of (a,1)

vidObj1 = VideoWriter('brownianmotionsim'); % include extension in string
vidObj1.FrameRate = 20; % set frame rate (e.g., 20)
open(vidObj1) % opens video for recording

%% Brownian Motion

x = BMFxVariDelta(a,b,delta,dt,pixel);
y = BMFxVariDelta(a,b,delta,dt,pixel);

% for j = 1:a
%     for k = 1:b
%     plot(x(j,1:k),y(j,1:k))
%     axis([0,512,0,512])
%     drawnow
%     end
% end
figure
for h = 1:a
    plot(x(h,1:b),y(h,1:b))
    axis([0,pixel,0,pixel])
    drawnow
    hold on
    currentFrame = getframe(gcf); % saves plot in variable
    writeVideo(vidObj1,currentFrame); % writes frame to video
    %     pause(dt)
end
hold off
title('2D Brownian Motion')
xlabel('x')
ylabel('y')

% while true
currentFrame = getframe(gcf); % saves plot in variable
writeVideo(vidObj1,currentFrame); % writes frame to video
% end

close(vidObj1)

%% Shift
xgood = x;
ygood = y;
shift = 0.15;
for kdx = 1:a
    %shift = rand; Creates random shift. You can assign one shift to all
    %particles or you can choose to have some particles have one shift
    %while others are random simply by using the main shift decleration and
    %by creating multiple for loops each with different kidx to indicate
    %which particles are assigned what
    for k = 2:b
        
        xgood(kdx,k) = xgood(kdx,k)+(k-1)*shift;
        ygood(kdx,k) = ygood(kdx,k)+(k-1)*shift;
    end
end
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

%%%% Generate Different Shift
% Inthe end we want to compare the correct version

%% Write Data

x2 = x';
y2 = y';

numidx = [0:a*b - 1]';
xdir = x2(:,1);

for k = 2:a
    
    xdir = [xdir;x2(:,k)];
    
    
end

ydir = y2(:,1);

for k = 2:a
    
    ydir = [ydir;y2(:,k)];
    
    
end

stepnum = [0:b-1]';

for k = 2:a
    
    stepnum = [stepnum;[0:b-1]'];
    
end

particle = 0 * ones(b,1);

for k = 1:a-1
    
    particle = [particle;k*ones(b,1)];
    
end


body_cells = num2cell([numidx,particle,stepnum,xdir,ydir]);
col_header = {'','particle','frame','x','y'};

sheets_final = [col_header;body_cells];
% xlswrite('bmsimdata.xlsx',sheets_final)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% xydir = [x(1,:);y(1,:)];
%
% for k = 2:a
%
%     xydir = [xydir;x(k,:);y(k,:)];
%
% end
%
% stepnum2 = [0:b-1];
%
%
% particle2 = 0 * ones(2,1);
%
% for k = 1:a-1
%
%     particle2 = [particle2;k*ones(2,1)];
% end
%
%
% col_header2 = {'particle'};
%
% sheets_final2 = [{''},num2cell(stepnum2);...
%     num2cell([particle2,xydir])];
% xlswrite('bmsim2.xlsx',sheets_final2)
%




%% Read Data
% partdata = readtable('06222020_100nm_G0F_Slide1_5ugml_02.tif.xlsx');
% endpart = partdata.particle(end);
% x = partdata.x;
% y = partdata.y;
% 
% %%%% Store as a cell!!
% 
% 
% figure
% for h = 0:endpart
%     k = find(partdata.particle == h);
%     %%% You can now do whatever you want ^^ indexes all of the x and y
%     %%% values
%     plot(x(k),y(k))
%     axis([0,pixel,0,pixel])
%     drawnow
%     hold on
%     
%     pause(0.5)
%     
% end
% hold off
% 
