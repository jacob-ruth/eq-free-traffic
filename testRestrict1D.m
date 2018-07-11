%% tests how accurate the restriction operator is when a data point is
% removed from the diffusion map

%% load diffusion map data
load('30dataAligned.mat', 'allData', 'alignData', 'evals', 'evecs', 'eps');
allData = alignData;

numRestrict = 50;
samplePoints = floor(linspace(1,size(evecs,1),numRestrict)); % points to lift and restrict
orig = zeros(numRestrict, 1);                   % to store the data
restricted = zeros(numRestrict, 1);
restrictDiff = zeros(numRestrict, 1);
origRestricted = zeros(numRestrict,1);
%nystrom error comparison, calculated the following way:
%   'correct' coordinates are given by the evec coord
%   approximated coordinates are given by diff map restrict, when the value
%   is removed from the data set
%   better approximated coordinates are given by diff map restrict, when
%   the value is still included in the data set
for i = 1:numRestrict
    disp(i);
    iRestrict = samplePoints(i);
    allData2 = allData;
    allData2(:, i) = [];                    % remove this point from the data set
    origRestricted(i) = evecs(iRestrict, :);
    evecs2 = evecs;
    evecs2(i) = []; 
    restricted(i) = diffMapRestrictAlt(allData(:,iRestrict),evals,evecs2,allData2,eps); % restrict this point
    restrictDiff(i) = norm((origRestricted(i) - restricted(i)));
end

figure; % plot the differences
scatter(evecs(samplePoints),restricted, 200, restrictDiff, '.');
colorbar;
title('Restricted Points Colored by Distance from Original Embedding', 'fontsize', 12);
xlabel('Original Coordinate', 'fontsize', 12);
ylabel('Restricted Coordinate', 'fontsize', 12);
