clear
clc

LB = [0.0001, 0.0001, 0.0001, 0.0001, 0.0001];   % Lower bound
UB = [1000 10 100 1 1]; % Upper bound


x = OrderedList()
%% Order
er bound
UB = [1000 10 100 1 1]; % Upper bound


x = OrderedList();
%% Order

x_n_avg = sum(x(:


%% Reflection



%% Expansion



%% Contraction



%% Shrink




%% Function
function y = objectiveFunction(x)
J = 0.1475;    %Inertial
B = 0.51;      %Damping
Tf = 2.2613;   %Friction

model = 'SBW_SystemModel';
in = Simulink.SimulationInput(model);
in = in.setVariable('x',x);
out = sim(in);
y = out.Error(end);
end
