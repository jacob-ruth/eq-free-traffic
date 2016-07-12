function res = linearFit2(evecs,k)
m = size(evecs,1);
PHI = [ones(m,1) evecs(:,1:k-1)]; %eigvectors are now columns
phi = evecs(:,k);
pdists = pdist(PHI(:,2:end));
ereg = median(pdists)/3;
sqpdists = squareform(pdists);

closeish = zeros(size(phi));

for i = 1:m
    curdists = sqpdists(i,:);
    dists = [curdists(1:i-1) , curdists(i+1:end)];
    W = diag(kernel(dists));
    curPHI = [PHI(1:i-1,:) ; PHI(i+1:end,:)];
    curphi = [phi(1:i-1,:) ; phi(i+1:end,:)];
    
    closeish(i) = PHI(i,:)*(((curPHI'*W*curPHI)\curPHI')*W*curphi);
end

res = norm(phi-closeish)/norm(phi);

    function dist = kernel(d)
        dist = exp(-d.^2/ereg^2);
    end
end