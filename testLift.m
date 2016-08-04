function testLift(evecs, evals, eps, allData)
tic;
steves = linspace(min(evecs), max(evecs), 100);
stevesLR  = steves;
lifts = zeros(100, 60);
for i = 1 : 100    
    i
    lifts(i,:) = diffMapLift(steves(i), evecs, evals, eps, allData);
    stevesLR(i) = diffMapRestrict(lifts(i,:)', evals, evecs, allData, eps);
end
norm(steves - stevesLR);
toc
end