%% The function runDiffMap takes in data from traffic simulations,
% builds a diffusion map, compares the eigen directions given by the
% diffusion map eigenvectors, and plots relevent eigenvectors scaled by
% corresponding eigenvalues
% INPUT:
% data - simulation data for building the diffusion map
% k - the number of eigenvectors to return
% RETURNS:
% vec - the eigenvectors from the diffusion map
% val - the eigenvalues from the diffusion map
% eps - the value of epsilon used by the diffusion map
function [vec,val,eps] = runDiffMap(data,k)

% build the data from the simulation results
allTime = data;

% calcuate the pairwise distances between data points
D = squareform(pdist(allTime'));
eps = 5*median(D(:)); % choose epsilon for the kernel based on the pairwise distances

[vec,val] = diffusionMap(eps,D,k);          % calculate the diffusion map

% calculate how unique each eigen direction is
%{
r = zeros(k, 1);
r(1) = 1;
for j = 2:k
    r(j) = linearFit(vec,j);
end
%}

    % plotEps creates a log-log plot of the number of data points that are
    % less than epsilon vs. epsilon
    % distances - the pairwise distances between data points
    function plotEps(distances)
        % create the values of epsilon to test from 0 to maxEps by stepSize
        epsilon = logspace(-3,2,1000);
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
end