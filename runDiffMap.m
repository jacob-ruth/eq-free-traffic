%% The function runDiffMap takes in data from traffic simulations,
% builds a diffusion map, compares the eigen directions given by the
% diffusion map eigenvectors, and plots relevent eigenvectors scaled by
% corresponding eigenvalues
function [vec,val,eps] = runDiffMap(data,k)

% build the data from the simulation results
allTime = data; % data should be passed in as headways

stdev = std(allTime, 0, 1);             % find the standard deviations

% calcuate the pairwise distances between data points
D = zeros(length(allTime));
for r = 1:length(allTime)
    for c = 1:length(allTime)
        D(r,c) = norm(allTime(:,r)-allTime(:,c));
    end
end

eps = median(D(:))/3; % choose epsilon for the kernel based on the pairwise distances

t = 3;      % number of timesteps to evolve the diffusion map
[vec,val] = diffusionMap(eps,D,k);          % calculate the diffusion map

% calculate how unique each eigen direction is
%{
r = zeros(k, 1);
r(1) = 1;
for j = 2:k
    r(j) = linearFit(vec,j);
end
% display the eigen direction computation results
fprintf('The r_k values are: \n');
disp(r);
%}

% evolve the diffusion map t times
%val = val^t;

%Change this to see the relationship between the eigenfunction at eigDisp
%and the standard deviation.
%{
eigDisp = 1;
stry = sprintf('Eigenvector %i',eigDisp);
figure;
hold on;
scatter(stdev, vec(:,eigDisp),'b.')
scatter(newstdev, dvec(eigDisp, :),500,'r.')
xlabel('\sigma');
ylabel(stry);
%}


% Plot relationship between first three eigenvectors
%{
figure;
hold on;
scatter(vec(:,1), vec(:,2),100,'b.');
scatter(dvec(1,:),dvec(2,:),100,'r.');
hold off;
xlabel('Eigenvector 1');
ylabel('Eigenvector 2');
%}


    % plotEps creates a log-log plot of the number of data points that are
    % less than epsilon vs. epsilon
    function plotEps(distances)
        % create the values of epsilon to test from 0 to maxEps by stepSize
        stepSize = .0001;
        maxEps = 1;
        epsilon = 0:stepSize:maxEps;
        % count the number of points that are less than each epsilon
        L = zeros(size(epsilon));
        for iEps = 1:length(epsilon)
            curEps = epsilon(iEps);
            L(iEps) = sum(sum(distances<curEps,1),2);
        end
        % plot the log-log plot of L(epsilon) vs epsilon
        figure;
        loglog(epsilon,L);
        xlabel('\epsilon','FontSize',24);
        ylabel('L(\epsilon)','FontSize',20);
    end

    % getHeadways returns headways between 60 cars given their positions
    function hways = getHeadways(v)
        futureCars = circshift(v,[-1,0]);
        hways = mod(futureCars - v, 60);
    end
end