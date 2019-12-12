P = 1;
E = 2;
h = 2;

x0 = [0,0];

options = optimset('Algorithm', 'active-set');
lb = [0, 0];
ub = [];
[x,f] = fmincon(@NonLinF, x0, [], [], [], [], ...
    lb, ub, @NonLinConst, options, P, E, h);


sprintf('%g - %g - %g  - %g', x)

function y = NonLinF(x, p, e, h)
y = ((p*h)/e)*(1/(x(1)+sqrt(2)*x(2)))
end

function [C,Ceq] = NonLinConst(x, p, e, h)
conU1 = (p*(x(2)+sqrt(2)*x(1))) ...
    /(sqrt(2)*x(1)^2+2*x(1)*x(2));
conU2 = (p) / (x(1)+ sqrt(2)*x(2));
conL = -p*(x(2)/(sqrt(2)*x(1)^2+2*x(1)*x(2)));
C(1) = conU1-17.5;
C(2) = conU2-17.5;
C(3) = conL-12;

Ceq = [];
end
