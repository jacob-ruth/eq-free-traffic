%% restriction operator based on diffusion map
% newData   - point to get restriction values for (column vector)
% evals     - square matrix of eigenvalues
% evecs     - matrix with eigenvectors in columns
% origData  - matrix of data in columns used to compute eigenvalues/vectors
% eps       - length scale for kernel
%
% returns diffusion map embedding for new data using the Nystrom extension
% in a column vector
function pnew = diffMapRestrict(newData,evals,evecs,origData,eps)
% pdist2 much faster than for loop. Finds pairwise distances between newData and each entry in origData.
dist = pdist2(newData',origData')';
w = basicKernel(dist);      % apply Gaussian kernel to weight smaller distances more heavily
k = (1/sum(w))*w;           % normalize weights to add to 1
pnew = (evecs' * k)./diag(evals);   % new embedding given by weighted sum of eigencoordinates already in the data

    function af = basicKernel(s)
        af = exp(-s.^2/eps^2);
    end
end

