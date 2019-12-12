clc;
clear;
close all;

tic;

J = 0.1475;    %Inertial
B = 0.51;      %Damping
Tf = 2.2613;   %Friction

%% Problem Definition

CostFunction= @multi_obj_func;      % Cost Function

nVar=5;             % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix

VarMin=0;          % Lower Bound of Variables
VarMax=1;          % Upper Bound of Variables


%% MOPSO Parameters

MaxIt=25;           % Maximum Number of Iterations

nPop=50;            % Population Size

nRep=10;            % Repository Size

w=0.5;              % Inertia Weight
wdamp=0.99;         % Intertia Weight Damping Rate
c1=1;               % Personal Learning Coefficient
c2=2;               % Global Learning Coefficient

alpha=0.1;          % Inflation Rate

beta=2;             % Leader Selection Pressure
gamma=2;            % Deletion Selection Pressure

mu=0.1;             % Mutation Rate

%% Initialization

empty_particle.Position=[];
empty_particle.Velocity=[];
empty_particle.Cost=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
empty_particle.IsDominated=false;

pop=repmat(empty_particle,nPop,1);
i = 1;
while i <= nPop
    %Generate Random Solution
    pop(i).Position=unifrnd(VarMin,VarMax,VarSize);
    %Initalize Velocity
    pop(i).Velocity=zeros(VarSize);
    %Evaulation
    pop(i).Cost=CostFunction(pop(i).Position);
 
    
    %Update Personal Best to its current location
    pop(i).Best.Position=pop(i).Position;
    %Update personal best cost
    pop(i).Best.Cost=pop(i).Cost;
    
    % Forice initial population to have valid cost
    if pop(i).Cost(1) < 9
        i= i+1;
    end
    
end

% Determine Domination
pop=evaluateDomination(pop);

rep=pop(~[pop.IsDominated]);

%% MOPSO Main Loop

for it=1:MaxIt
    
    for i=1:nPop
        
        leader=rep(randi([1 numel(rep)]));
        
        pop(i).Velocity = w*pop(i).Velocity ...
            +c1*rand(VarSize).*(pop(i).Best.Position-pop(i).Position) ...
            +c2*rand(VarSize).*(leader.Position-pop(i).Position);
        
        pop(i).Position = pop(i).Position + pop(i).Velocity;
        
        pop(i).Position = max(pop(i).Position, VarMin);
        pop(i).Position = min(pop(i).Position, VarMax);
        
        pop(i).Cost = CostFunction(pop(i).Position);
        
        % Dismiss particles with error over 9
%         if pop(i).Cost(1) > 9
%             continue
%         end
        
        if dominates(pop(i),pop(i).Best)
            pop(i).Best.Position=pop(i).Position;
            pop(i).Best.Cost=pop(i).Cost;
            
        elseif dominates(pop(i).Best,pop(i))
            % Do Nothing
            
        else
            if rand<0.5
                pop(i).Best.Position=pop(i).Position;
                pop(i).Best.Cost=pop(i).Cost;
            end
        end
        
    end
    
    % Dismiss particles with error over 9
%     a = arrayfun(@(n) n.Cost(1) < 9, pop);
%     pop = pop(a);
    
    % Add Non-Dominated Particles to REPOSITORY
    rep=[rep
         pop(~[pop.IsDominated])]; %#ok
    
    % Determine Domination of New Resository Members
    rep=evaluateDomination(rep);
    
    % Keep only Non-Dminated Memebrs in the Repository
    rep=rep(~[rep.IsDominated]);
    
    % remove random particle if rep is full
    while numel(rep)>nRep
        rep(randi([1 numel(rep)]))=[];        
    end
    
    % Plot Costs
    figure(1);
    xlabel('Error');
    ylabel('Consumption');
    %Plot current population
%     pop_costs=[pop.Cost];
%     plot(pop_costs(1,:),pop_costs(2,:),'ko');
%     hold on;
    %Plot archive of points
    rep_costs=[rep.Cost];
    plot(rep_costs(1,:),rep_costs(2,:),'r*');
%     hold off;
    grid on

    
    % Show Iteration Information
    disp(['Iteration ' num2str(it) ': Number of Rep Members = ' num2str(numel(rep))]);
    
    % Damping Inertia Weight
    w=w*wdamp;
    
end

toc

function y = multi_obj_func(x)
model = 'SBW_SystemModel';
% model = 'test2';
J = 0.1475;    %Inertial
B = 0.51;      %Damping
Tf = 2.2613;   %Friction

in = Simulink.SimulationInput(model);
in = in.setVariable('x',x);

out = sim(in);
% y = out.simout(end);
f1 = out.Error(end);
f2 = out.Consumption(end);
y = [f1
    f2];
end

function b = dominates(x,y)
    if isstruct(x)
        x=x.Cost;
    end
    if isstruct(y)
        y=y.Cost;
    end
    % Atleast one of the inequalities must be strict, others can be equal
    % or smaller
    b = any(x<y) && all(x<=y);

end

function population = evaluateDomination(population)
    n=numel(population);
    % Reset all dominated values to false before beginning
    for i=1:n
        population(i).IsDominated=false;
    end
    
    % Pairwise evaluation of each combination
    for i=1:n-1
        for j=i+1:n
            if dominates(population(i),population(j))
               population(j).IsDominated=true;
            end
            
            if dominates(population(j),population(i))
               population(i).IsDominated=true;
            end
            
        end
    end
end