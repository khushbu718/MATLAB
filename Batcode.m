% Parameters
tic
n = 100;                % Number of bats
maxIterations = 100;   % Maximum number of iterations
A = 0.5;               % Loudness
r = 0.5;               % Pulse rate
Qmin = 0;              % Minimum frequency
Qmax = 2;              % Maximum frequency
lowerBound = [0.1 0.1];     % Lower search space boundary
upperBound = [2 3.8];       % Upper search space boundary

% Initialize bats
position = zeros(n, 2);    % Bat positions (2 variables)
velocity = zeros(n, 2);    % Bat velocities (2 variables)
frequency = zeros(n, 1);   % Bat frequencies
loudness = zeros(n, 1);    % Bat loudness
pulseRate = zeros(n, 1);   % Bat pulse rates
bestPosition = zeros(1, 2);% Best bat position (2 variables)
bestFitness = inf;         % Best bat fitness

% Initialize bats randomly-- all these parameters are at time = 0 s
for i = 1:n
    position(i,:) = lowerBound + (upperBound - lowerBound).* rand;
    velocity(i,:) = zeros(1,2);
    frequency(i) = Qmin + (Qmax - Qmin) * rand;
    loudness(i) = A;
    pulseRate(i) = r;
end
% Create the figure and axes
% figure;
% ax = gca;
% ax.XLim = [min(position(:,1))-1, max(position(:,1))+1];
% ax.YLim = [min(position(:,2))-1, max(position(:,2))+1];
% hold on;
% quiverPlot = quiver(position(:,1), position(:,2), velocity(:,1), velocity(:,2), 0.5, 'b', 'LineWidth', 1.5);
% scatterPlot = scatter(position(:,1), position(:,2), 50, 'r', 'filled');
% hold off;
% xlabel('X');
% ylabel('Y');
% title('Position and Velocity Vectors');
%legend('Velocity', 'Position');
%% Bat code
  bestFitnessValues = zeros(maxIterations, 1);   
for t = 1:maxIterations
    % Update bats
    for i = 1:n
        % Generate a random solution
        epsilon = randn(1, 2);
        
        % Update bat frequency
        frequency(i) = Qmin + (Qmax - Qmin) * rand;
        
        % Update bat velocity
        velocity(i,:) = velocity(i,:) + (position(i,:) - bestPosition) * frequency(i);
        
        % Update bat position
        position(i,:) = position(i,:) + velocity(i,:);
        
        % Apply random walk
        position(i,:) = position(i,:) + epsilon;
        
        % Apply bounds
        position(i,:) = max(position(i,:), lowerBound);
        position(i,:) = min(position(i,:), upperBound);
        
        % Evaluate fitness
        fitness = funbat1(position(i,:));
        
        % Update bat loudness and pulse rate
        if rand > pulseRate(i)
            loudness(i) = A * rand;
            pulseRate(i) = r * (1 - exp(-0.1 * t));
        end
        
        % Update the best bat
        if fitness < bestFitness
            bestFitness = fitness;
            bestPosition = position(i,:);
        end
   
    end
    %    set(scatterPlot, 'XData', position(:,1), 'YData', position(:,2));
    %   set(quiverPlot, 'XData', position(:,1), 'YData', position(:,2), 'UData', velocity(:,1), 'VData', velocity(:,2));
    %  drawnow;
    
    % Display best fitness value at each iteration
    disp(['Iteration' num2str(t) ': Best Fitness = ' num2str(bestFitness)]);
    % disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
    
   bestFitnessValues(t) = bestFitness;
end
H1=[bestFitnessValues]';
% figure;
% plot(bestFitness,'r', 'LineWidth',2);
% %semilogy(bestFitness,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost');
% grid on;
toc
% Plot the convergent graph
figure;
plot(1:maxIterations, bestFitnessValues, 'r', 'LineWidth', 2);
xlabel('Iteration');
ylabel('Best Fitness');
title('Convergence Graph');
grid on;

% Display the final best solution
%disp('------------------------');
%disp('Final Result:');
disp(['Best Fitness:' num2str(bestFitness)]);
disp(['Best Position: ' num2str(bestPosition)]);