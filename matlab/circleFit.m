function [xc,yc,R] = circleFit(x,y)

x = x(:);
y = y(:);

A = [2*x 2*y ones(length(x),1)];
b = x.^2 + y.^2;

params = A\b;

xc = params(1);
yc = params(2);
c  = params(3);

R = sqrt(c + xc^2 + yc^2);

end