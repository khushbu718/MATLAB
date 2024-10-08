tic
clc;
clear;
close all;
%% Problem Definition
CostFunction=@funpso1;     % Cost Function
nVar=2;                   % Number of Decision Variables
VarSize=[1 nVar];         % Size of Decision Variables Matrix
VarMin=[0.1 0.1];          % Lower Bound of Variables
VarMax= [2 3];             % Upper Bound of Variables
%% PSO Parameters
MaxIt=100;       % Maximum Number of Iterations
nPop=100;        % Population Size (Swarm Size)
% PSO Parameters
w=1;            % Inertia Weight
wdamp=0.99;     % Inertia Weight Damping Ratio
c1=1.5;         % Personal Learning Coefficient
c2=2.0;         % Global Learning Coefficient
% If you would like to use Constriction Coefficients for PSO,
% uncomment the following block and comment the above set of parameters.
% % Constriction Coefficients
% phi1=2.05;
% phi2=2.05;
% phi=phi1+phi2;
% chi=2/(phi-2+sqrt(phi^2-4*phi));
% w=chi;          % Inertia Weight
% wdamp=1;        % Inertia Weight Damping Ratio
% c1=chi*phi1;    % Personal Learning Coefficient
% c2=chi*phi2;    % Global Learning Coefficient
% Velocity Limits
%total_time_start = tic;
VelMax=0.1*(VarMax-VarMin);
VelMin=-VelMax;
%% Initialization
empty_particle.Position=[];
empty_particle.Cost=[];
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];
particle=repmat(empty_particle,nPop,1);
GlobalBest.Cost=inf;
for i=1:nPop
    
    % Initialize Position
    %particle(i).Position=unifrnd(VarMin,VarMax,VarSize);
    particle(i).Position(1)=unifrnd(VarMin(1),VarMax(1));
    particle(i).Position(2)=unifrnd(VarMin(2),VarMax(2));
    
    
    % Initialize Velocity
    particle(i).Velocity=zeros(VarSize);
    
    % Evaluation
    particle(i).Cost=CostFunction(particle(i).Position);
    
    % Update Personal Best
    %particle(i).Best.Position=particle(i).Position;
    particle(i).Best.Position(1)=(particle(i).Position(1));
    particle(i).Best.Position(2)=particle(i).Position(2);
    
    particle(i).Best.Cost=particle(i).Cost;
    
    % Update Global Best
    if particle(i).Best.Cost<GlobalBest.Cost
        
        GlobalBest=particle(i).Best;
        
    end
    
end
BestCost=zeros(MaxIt,1);

% Start the timer



%% PSO Main Loop
for it=1:MaxIt
    
    for i=1:nPop
        
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(GlobalBest.Position-particle(i).Position);
        
        % Apply Velocity Limits
        particle(i).Velocity = max(particle(i).Velocity,VelMin);
        particle(i).Velocity = min(particle(i).Velocity,VelMax);
        
        % Update Position
        %particle(i).Position = particle(i).Position + particle(i).Velocity;
        particle(i).Position(1) = (particle(i).Position(1) + particle(i).Velocity(1));
        particle(i).Position(2) = particle(i).Position(2) + particle(i).Velocity(2);
        
        
        % Velocity Mirror Effect
        IsOutside=(particle(i).Position<VarMin | particle(i).Position>VarMax);
        particle(i).Velocity(IsOutside)=-particle(i).Velocity(IsOutside);
        
        % Apply Position Limits
        particle(i).Position = max(particle(i).Position,VarMin);
        particle(i).Position = min(particle(i).Position,VarMax);
        
        % Evaluation
        particle(i).Cost = CostFunction(particle(i).Position);
        
        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            
            % Update Global Best
            if particle(i).Best.Cost<GlobalBest.Cost
                
                GlobalBest=particle(i).Best;
                
            end
            
        end
        
    end
    
    BestCost(it)=GlobalBest.Cost;
    
    % Stop the timer and get the elapsed time


    
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
   
   
    w=w*wdamp;
    
end
BestSol = GlobalBest;
%total_time_elapsed = toc(total_time_start);
 H=[BestCost]';
%disp(['Total Time: ', num2str(total_time_elapsed), ' seconds']);
%% Results
% [x,fval,exitflag] = particleswarm(CostFunction,nVar)
toc
figure;
set(gca,'FontSize',14,'FontName','Times New Roman','FontWeight','Bold');
%ylim([20 30]) 
%xlim([1 6])
ylabel('Best Cost','FontSize',14,'FontName','Times New Roman', 'Color','k')
xlabel('Iteration', 'FontSize',14,'FontName','Times New Roman','Color','k')

plot(BestCost,'LineWidth',2);
semilogy(BestCost,'LineWidth',2);
xlabel('Iteration');
ylabel('Best Cost');
grid on;
