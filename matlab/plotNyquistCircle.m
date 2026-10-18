function plotNyquistCircle(H,f,frange,titleStr)

%% Select resonance

idx = f>=frange(1) & f<=frange(2);

Hmode = H(idx);

x = real(Hmode);
y = imag(Hmode);

%% Fit Circle

[xc,yc,R] = circleFit(x,y);

theta = linspace(0,2*pi,400);

xfit = xc + R*cos(theta);
yfit = yc + R*sin(theta);

%% RMS fitting error

dist = sqrt((x-xc).^2 + (y-yc).^2);

rmse = sqrt(mean((dist-R).^2));

%% Plot

figure
hold on

plot(x,y,'b','LineWidth',2)

plot(xfit,yfit,'r--','LineWidth',2)

plot(xc,yc,'ko','MarkerFaceColor','k')

axis equal
grid on

xlabel('Real(H)')
ylabel('Imag(H)')
title(titleStr)

legend('Experimental FRF','Circle Fit','Center','Location','northeast')

% text(min(x),max(y),sprintf('R = %.3e\nRMSE = %.3e',R,rmse),...
%     'FontSize',11,...
%     'BackgroundColor','w')

%% Circle fit quality

dist = sqrt((x-xc).^2 + (y-yc).^2);

RMSE = sqrt(mean((dist-R).^2));

SSres = sum((dist-R).^2);

SStot = sum((dist).^2);

R2 = 1 - SSres/SStot;
sprintf('Radius = %.4e\nRMSE = %.3e\nR^2 = %.4f',...
R,RMSE,R2)
end