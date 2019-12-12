%Dila T�rkmen - 2014
%Nelder and Mead simplex algorithm
clear 
clc
% fname=input('Please enter the name of the function you have defined:','s');
fname = 'singleObjFunc';
f=str2func(fname);
dila=0;


%% Parameters

J = 0.1475;    %Inertial
B = 0.51;      %Damping
Tf = 2.2613;   %Friction
% n=input('enter the number of dimensions of the problem,n:');
% x0=input('enter a starting point vector, x0=[....]:');
n = 5;
x0 = [500 1 20 0 0];
% alpha=input('alpha(recommended value=1):');
alpha2 = 1;
% beta=input('beta(recommended value=0.5):');
beta = 0.5;
% gamma=input('gamma(recommended value=2):');
gamma = 2;
% delta=input('delta(recommended value=0.5):');
delta = 0.5;
% lambda=input('lambda(recommended value=1):');
lambda = 1;
% eps1=input('eps1 i.e.:10^-6:');
eps1 = 1e-2;
% eps2=input('eps2 i.e.:10^-6:');
eps2 = 1e-2;


K=1;%Counter is set to 1.

for K=1:inf %%%%%%%%%%%%%%%-Will be returned if a corner is constant for m times.-%%%%%%
    S=zeros(1,n+1);sdila=0;% For the 6th step
    
    %% 2-Generation of the initial simplex.

    D1=lambda*(sqrt(n+1)+n-1)/sqrt(2)*n;
    D2=lambda*(sqrt(n+1)-1)/sqrt(2)*n;
    x=zeros(n,n+1);
    x(:,1)=x0;
    for k=2:n+1
        for i=1:n
            if i==k
                x(i,k)=x0(i)+D1;
            else
                x(i,k)=x0(i)+D2;
            end
        end
    end
    xbackup=x;
    % Simplex corners are defined.

    for KK=1:inf %%%%%%%%%%%%%%%-Will be continue here if the corners are changing-%%%%%%
        %% 3- Findings
        F=zeros(1,n+1);
        for i=1:n+1
%             x
%             i
%             x(:,i)
%             f(x(:,i))
            F(i)=f(x(:,i));
        end
        Fmin=min(F);
        Fbackup=F;
        Fmax=max(F);
        if n>1
            [ignore, index] = max(F); 
            F(index) = -Inf;
            Fmax2 = max(F);
        end
        
        %% 4- Convergence Control

        b=2;
        for j=1:n+1
            for k=b:n+1
                Diff=x(:,j)-x(:,k);
                u=0;
                for s=1:n
                    u=u+Diff(s)^2;
                end
                Length=sqrt(u);
                if Length<eps1 && abs(f(x(:,j))-f(x(:,k)))<eps2
                    for i=1:n+1
                        if f(x(:,i))==Fmin
                            xmin=x(:,i);
                            break
                        end
                    end
                    dila=1;
                    xmin
                    Fmin=f(xmin)
                    imin=i;
                    break 
                end
            end
            b=b+1;
            if b>n+1
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
        for i=1:n+1
            if f(x(:,i))==Fmax
                xmax=x(:,i);imax=i;
                break
            end
        end
        xs=0;
        for i=1:n+1
            if i==imax
                xs=xs;
            else
                xs=xs+x(:,i);
            end
        end
        xs=xs/n;
        xnew=xmax+2*(xs-xmax);
        Fnew=f(xnew);
        if n==1
            if Fmin<=Fnew && Fnew<Fmax
                theta=beta;
            else if Fnew>=Fmax
                    theta=-beta;
                else if Fnew<=Fmin
                        theta=gamma;
                    end
                end
            end
        else if Fmin<Fnew && Fnew<Fmax2
            theta=alpha2;
        else if Fmax2<=Fnew && Fnew<Fmax
                theta=beta;
            else if Fnew>=Fmax
                    theta=-beta;
                else if Fnew<=Fmin
                        theta=gamma;
                    end
                end
            end
            end
        end
        x(:,imax)=x(:,imax)+(1+theta)*(xs-x(:,imax));

        %% 6- Shrinkage of the simplex

        change=x-xbackup;
        for i=1:n+1
            if change(:,i)==zeros(n,1)
                S(i)=S(i)+1;
            end
        end
        m=round(1.65*n+.05*n^2);
        for i=1:n+1
            if S(i)==m
                sdila=1;
                break
            end
        end
        if sdila==1 || dila==1
            break
        end

    end% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% if the corners change %%%%%%%%%%%%%%%%%

    if dila==1
        break
    end
    lambda=delta*lambda;

end% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% if a corner is constant for m times %%%%%%%%%%%%