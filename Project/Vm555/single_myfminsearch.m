clear 
clc

func = @singleObjFunc;

% Simulink model parameters
J = 0.1475;    %Inertial
B = 0.51;      %Damping
Tf = 2.2613;   %Friction

x_0 = [200 5 50 0.5 0.5]';
n = numel(x_0);

alph = 1;      %alpha rho - reflection parameter (a > 0, default 1)
beta = 0.5;    %beta psi - contraction parameter (0 < b < 1, default 0.5)
gamma = 2;      %gamma chi = expansion paramater (y > 1, default 2)
delta = 0.5;    %delta sigma - shrinkage parameter (0 < d < 1, default 0.5)
eps = 1e-6;

%% Initial simplex setup - step along each dimension
extra_points = repmat(x_0, 1,n);
extra_points = extra_points + extra_points .* (0.05 * eye(n)); % simplex 
s = [x_0 extra_points];
F = zeros(1,n+1);
%matrix approach seems to not work so use forloop
% columns = num2cell(s);
% cellfun(func,columns)
for i=1:n+1
    F(i) = func(s(:,i));
end



%% Begin loop
round = 0;
last_action = 'initial simplex';
while true
    round = round+1;
    % Sort on ascending function value
    [F, index] = sort(F);
    s = s(:,index);
    
    %% Evaluate termination
%     sprintf('round %d:  x: [%s], F: %g, eps: %g, last action: %s', round, 'todo', F(n+1),max(abs(F(1)-F(2:n+1))), last_action)
    g=sprintf('%.4f\t\t', s(:,n+1));
    fprintf('%d \t%g \t%s\n', round, max(abs(F(1)-F(2:n+1))), g)
    if max(abs(F(1)-F(2:n+1))) <= eps % epsilon between worst points is satisfied
        disp('Terminated')
        break
    end
    
    %% Compute reflection point
    centroid = mean(s(:,1:n), 2); % mean over axis 2
    x_refl = (1+alph)*centroid - alph*s(:,n+1);
    F_refl = func(x_refl);
    
    
    if F_refl < F(1)
        % replace s(n+1) with best of either x_refl or x_exp
        
        %% Expansion point
        x_exp = (1 + gamma)*centroid - gamma*s(:,n+1);
        F_exp = func(x_exp);
        
        if F_exp < F(1)
            s(:,n+1) = x_exp;
            F(n+1) = F_exp;
        else
            s(:,n+1) = x_refl;
            F(n+1) = F_refl;
        end
        
    else
        if F_refl < F(n)
            s(:,n+1) = x_refl;
            F(n+1) = F_refl;
        else 
            %% Contraction
            if F_refl < F(n+1)
                s(:,n+1) = x_refl;
            end
            x_cont = beta * s(:,n+1) + (1 - beta) * centroid;
            F_cont = func(x_cont);
            if F_cont < F(n+1)
                s(:,n+1) = x_cont;
                F(n+1) = F_cont;
            else
                %% Shrink
                x_min = s(:,1);
                F_min = F(1);
                for i_shrink=1:n
%                     s(:,i_shrink) = (s(:,i_shrink)+x_min)/2;
                    s(:,i_shrink) = delta*x_min + (1-delta)*s(:,i_shrink);

                end
            end  
        end
        
    end
end
g=sprintf('%.4f\t', s(:,n+1));
sprintf('round %d: x: %s\tF:%g', round, g, F(n+1))
