clear
clc

% Plot two objective functions on the same axis
% x = -10:0.5:10;
% f1 = (x+2).^2 - 10;
% f2 = (x-2).^2 + 20;
% plot(x,f1);
% hold on;
% plot(x,f2,'r');
% grid on;
% title('Plot of objectives ''(x+2)^2 - 10'' and ''(x-2)^2 + 20''');

numberOfVariables = 1;

%%% Constraints
% inEqCon * x <= inEqConVal
inEqCon = [];
inEqConVal = [];
% eqCon * x = eqConVal
eqCon = [];
eqConVal = [];

%%% Bounds
lb = -1.5;
ub = 0;

options = optimoptions(@gamultiobj,'PlotFcn',{@gaplotpareto,@gaplotscorediversity}); %Visualization
options = optimoptions(options,'FunctionTolerance',1e-3,'MaxStallGenerations',150); %Stopping criteria

[x,fval] = gamultiobj(@funcGaMult,numberOfVariables, ...
    inEqCon,inEqConVal,eqCon,eqConVal, ...
    lb,ub, ...
    options);

fprintf('The number of points on the Pareto front was: %d\n', size(x,1));

function y = funcGaMult(x) %Fitness Function
    y(1) = (x+2).^2 - 10;
    y(2) = (x-2).^2 + 20;
end
