clear 
clc

f = @singleObjFunc;
dila=0;


%% Parameters

% Simulink model parameters
J = 0.1475;    %Inertial
B = 0.51;      %Damping
Tf = 2.2613;   %Friction

nrParam = 5;
initial_point = [500 1 20 0 0];

alph = 1;      %alpha - reflection parameter (a > 0, default 1)
beta = 0.5;    %beta - contraction parameter (0 < b < 1, default 0.5)
gamma = 2;      %gamma = expansion paramater (y > 1, default 2)
delta = 0.5;    %delta - shrinkage parameter (0 < d < 1, default 0.5)
l = 1;

eps = 1e-2;

%% Algorithm
for K=1:intmax
    S=zeros(1,nrParam+1);
    sdila=0;% For the 6th step
    
    
    %% Initial Simplex
    D1=l*(sqrt(nrParam+1)+nrParam-1)/sqrt(2)*nrParam;
    D2=l*(sqrt(nrParam+1)-1)/sqrt(2)*nrParam;
    x=zeros(nrParam,nrParam+1);
    x(:,1)=initial_point;
    
    % Generate n+1 vertices around an the initial point
    for k=2:nrParam+1
        for i=1:nrParam
            if i==k
                x(i,k)=initial_point(i)+D1;
            else
                x(i,k)=initial_point(i)+D2;
            end
        end
    end
    x_old=x;

    for KK=1:intmax 
        % Compute the values of the vertices
        F = zeros(1,nrParam+1);
        for i=1:nrParam+1
            F(i)=f(x(:,i));
        end
        
        %         Fbackup=F;
        F_min=min(F);   % Best point        
        F_max=max(F);   % Worst point
        
        if nrParam>1
            [ignore, index] = max(F); 
            F(index) = -Inf;
            Fmax2 = max(F);
        end
        
        %% 4- Convergence Control

        b=2;
        for j=1:nrParam+1
            for k=b:nrParam+1
                Diff=x(:,j)-x(:,k);
                u=0;
                for s=1:nrParam
                    u=u+Diff(s)^2;
                end
                Length=sqrt(u);
                if Length<eps && abs(f(x(:,j))-f(x(:,k))) < eps
                    for i=1:nrParam+1
                        if f(x(:,i))==F_min
                            xmin=x(:,i);
                            break
                        end
                    end
                    dila=1;
                    xmin
                    F_min = f(xmin)
                    imin=i;
                    break 
                end
            end
            b=b+1;
            if b>nrParam+1
                break
            end
            if dila==1
                break
            end
            
        end
        if dila==1
            break
        end

        %% 5- Renewal of the worst corner by adjusting the step scale.
        for i=1:nrParam+1
            if f(x(:,i))==F_max
                x_max=x(:,i);
                i_max=i;
                break
            end
        end
        xs=0;
        for i=1:nrParam+1
            if i==i_max
                xs=xs;
            else
                xs=xs+x(:,i);
            end
        end
        xs=xs/nrParam;
        xnew=x_max+2*(xs-x_max);
        
        
        F_new=f(xnew);
        
        if nrParam==1
            if F_min<=F_new && F_new<F_max
                theta=beta;
            else if F_new>=F_max
                    theta=-beta;
                else if F_new<=F_min
                        theta=gamma;
                    end
                end
            end
        else if F_min<F_new && F_new<Fmax2
            theta=alph;
        else if Fmax2<=F_new && F_new<F_max
                theta=beta;
            else if F_new>=F_max
                    theta=-beta;
                else if F_new<=F_min
                        theta=gamma;
                    end
                end
            end
            end
        end
        x(:,i_max)=x(:,i_max)+(1+theta)*(xs-x(:,i_max));

        %% Shrink

        x_diff=x-x_old;
        for i=1:nrParam+1
            if x_diff(:,i)==zeros(nrParam,1)
                S(i)=S(i)+1;
            end
        end
        m=round(1.65*nrParam+.05*nrParam^2);
        for i=1:nrParam+1
            if S(i)==m
                sdila=1;
                break
            end
        end
        if sdila==1 || dila==1
            break
        end

    end

    if dila==1
        break
    end
    
    % lambda = delta * lambda
    l=delta*l;

end