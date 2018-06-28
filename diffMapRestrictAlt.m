%% restriction operator based on diffusion map
% newData   - point to get restriction values for (column vector)
% evals     - square matrix of eigenvalues
% evecs     - matrix with eigenvectors in columns
% origData  - matrix of data in columns used to compute eigenvalues/vectors
% eps       - length scale for kernel
%
% returns diffusion map embedding for new data using the Nystrom extension
% in a column vector
function pnew = diffMapRestrictAlt(newData,evals,evecs,origData,eps,fullData)
newDist = pdist2(newData',origData')';     % calculate the pairwise distances between newData and origData
newKernel = basicKernel(newDist);

if fullData
    % Obtains the sums of the kernel for the distances between any element
    % the diffusion map matrix and all others without having to compute it
    % every time
    load('sumsSelfKernel.mat', 'sumsSelfKernel');
else
    selfDist = pdist2(origData',origData')';
    selfKernel = basicKernel(selfDist);
    sumsSelfKernel = sum(selfKernel);
end

% Nystrom extension calculation
try
    distRatio = (newKernel' ./ sqrt(sumsSelfKernel)) / sqrt(sum(newKernel));
catch
    disp('oops');
end
pnew = (evecs' * distRatio') ./ diag(evals);

    function af = basicKernel(s)
        af = exp(-s.^2/eps^2);
    end

end