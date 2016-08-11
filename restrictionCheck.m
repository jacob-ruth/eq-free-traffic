load('bigDataMap.mat', 'evecs', 'evals', 'eps', 'allData');
numRestrict = 500;
restrictDiff = zeros(numRestrict, 2);
restrictDiffS = zeros(numRestrict, 2);
%nystrom error comparison, calculated the following way:
%   'correct' coordinates are given by the evec coord
%   approximated coordinates are given by diff map restrict, when the value
%   is removed from the data set
%   better approximated coordinates are given by diff map restrict, when
%   the value is still included in the data set
for iRestrict = 1:numRestrict
    evecs2 = [evecs(1:iRestrict - 1, :) ; evecs(iRestrict + 1 : end, :)]; 
    allData2 = [allData(:, 1:iRestrict - 1) , allData(:,iRestrict + 1 : end)]; 
    restricted = diffMapRestrict(allData(:,iRestrict),evals,evecs2,allData2,eps);
    restrictedS = diffMapRestrict(allData(:,iRestrict),evals,evecs,allData,eps);
    origRestricted = evecs(iRestrict, :);
    restrictDiff(iRestrict,:) = origRestricted - restricted';
    restrictDiffS(iRestrict,:) = origRestricted - restrictedS';
end
disp('i did a thing');
