clear
clc

J = 0.1475;    %Inertial
B = 0.51;      %Damping
Tf = 2.2613;   %Friction


ObjectiveFunction = @simple_fitness;
nvars = 5;    % Number of variables


%%% Constraints
% inEqCon * x <= inEqConVal
inEqCon = [];
inEqConVal = [];
% eqCon * x = eqConVal
eqCon = [];
eqConVal = [];

lb = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001];   % Lower bound
ub = [1000 10 100 1 1]; % Upper bound


options = optimoptions(@gamultiobj);
options = optimoptions(options, 'PlotFcn',{@gaplotpareto,@gaplotscorediversity}); %Visualization
% options = optimoptions(options,'FunctionTolerance',1e-3,'MaxStallGenerations',150); %Stopping criteria
options = optimoptions(options, 'PopulationSize', 200);

[x,fval] = fminimax(@funcGaMult,nvars, ...
    inEqCon,inEqConVal,eqCon,eqConVal, ...
    lb,ub, ...
    options)

fprintf('The number of points on the Pareto front was: %d\n', size(x,1));

function y = funcGaMult(x) %Fitness Function
model = 'SBW_SystemModel';
% open_system(model);

in = Simulink.SimulationInput(model);

in = in.setVariable('x',x);
out = sim(in);
y(1) = out.Error(end);
y(2) = out.Consumption(end);
end
