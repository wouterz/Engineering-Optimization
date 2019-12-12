clear

spacesteps = 20;
th_start = 28/180*pi;
th_end = 115/180*pi;
[th,h] = meshgrid(linspace(th_start, th_end, spacesteps), linspace(1, 4, spacesteps));
A=9.3;
f = 1./(A./h-h.*cot(th)+2*h./sin(th));
figure;
contour(h,th*180/pi,f,'k');
xlabel('h');
ylabel('\theta (degrees)');

figure;
mesh(h, th*180/pi, f);
xlabel('h');
ylabel('\theta (degrees)');
zlabel('1/p');


opt = optimset();
[xopt, fopt] = fminsearch(@func, [1,1], opt, A);

sprintf('h opt: %g', xopt(1))
sprintf('theta opt %g', fopt*pi*180)

function y = func(x_0, A)
y = 1./(A./x_0(1)-x_0(1).*cot(x_0(2))+2*x_0(1)/sin(x_0(2)));
end