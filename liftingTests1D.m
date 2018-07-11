%% tests the identity R(L) = I by taking coordinates in the diffusion map 
% embedding, lifting them, restricting them, and looking at the difference

% load diffusion map data
load('30dataAligned.mat', 'allData', 'alignData', 'evals', 'evecs', 'eps');
vel = allData(:,end);
allData = alignData;                                 % aligned headway data

numRestrict = 50;
samplePoints = floor(linspace(1,size(evecs,1),numRestrict)); % points to lift and restrict
orig = zeros(numRestrict, 1);
restricted = zeros(numRestrict, 1);
restrictDiff = zeros(numRestrict, 1);

for i=1:numRestrict
    disp(i);
    iRestrict = samplePoints(i);
    orig(i) = evecs(iRestrict);
    lifted = newLift(orig(i), evecs, evals, eps, vel(iRestrict), allData);  % lift the profile
    restricted(i) = diffMapRestrictAlt(getHeadways(lifted(1:30), 60),evals,evecs,allData,eps); % restrict the lifted profile
    restrictDiff(i) = norm(restricted(i) - orig(i));    % compute the difference from the original embedding
end

figure; % plot the differences
save('liftedTestData.mat', 'restricted', 'restrictDiff');
hold on;
scatter(orig,restricted, 200, restrictDiff, '.');
colorbar;
title('Lifted and Restricted Points Colored by Distance from Original Embedding', 'fontsize', 18);
xlabel('\psi_1', 'fontsize', 20);
ylabel('\psi_2', 'fontsize', 20);