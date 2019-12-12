clear
clc
P = 1;
E = 2;
h = 2;

options = optimoptions(options,'PlotFcn',{@gaplotbestf,@gaplotmaxconstr}, ...
    'Display','iter');

% ObjectiveFunction = @simple_fitness;
nvars = 2;    % Number of variables
LB = [0 0];   % Lower bound
UB = [1 13];  % Upper bound
% ConstraintFunction = @simple_constraint;
[x,fval] = ga(NonLinF,nvars,[],[],[],[],LB,UB, ...
    NonLinConst, options, P, E, h)
% options = optimoptions(@ga,'MutationFcn',@mutationadaptfeasible);
% X0 = [0.5 0.5]; % Start point (row vector)
% options.InitialPopulationMatrix = X0;
% 
% % Next we run the GA solver.
% [x,fval] = ga(NonLinF,nvars,[],[],[],[],LB,UB, ...
%     NonLinConst,options)
% 
% % Next we run the GA solver.
% [x,fval] = ga(NonLinF,nvars,[],[],[],[],LB,UB, ...
%     NonLinConst,options)

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

