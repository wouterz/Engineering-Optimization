close all;
clc;
clear

% zfun1 = @(x1,x2) (x1-0.75).^2 + (x2-2).^2;
% zfun2 = @(x1,x2) (x1-2.5).^2 + (x2-1.5).^2;


zfun1 = @(x1,x2) (x1-2).^2 + (x2-5).^2;
zfun2 = @(x1,x2) (x1-4.5).^2 + (x2-8.5).^2;

% zhandle = fcontour(zfun1)

[x,y] = meshgrid(-10:0.1:10,-10:0.1:10);

z1 = zfun1(x,y);
z2 = zfun2(x,y);

% contour(x,y,z1, [0:3:100])
% hold on
% contour(x,y,z2,[0:3:100])
% grid on

mesh(x,y,z1)
hold on
mesh(x,y,z2)

