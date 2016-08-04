function [vec,val] = diffusionMap(epsilon,distMatrix,k)
    %find Markov matrix and its eigenvalues/eigenvectors
    A = basicKernel(distMatrix);
    M = markovify(A);
    [eigenvec,eigenval] = eigs(M,k+1);
    
    %sort eigenvalues and eigenvectors
    seval = diag(sort(diag(eigenval),'descend'));
    [c,ind] = sort(diag(eigenval),'descend');
    sevec = eigenvec(:,ind);
    seval = sparse(seval);
    eigsigns = sevec(1,:)./abs(sevec(1,:));
    sevec = sevec * diag(eigsigns);
    
    %return first k eigenvalues/vectors
    val = seval(2:end,2:end);
    vec = sevec(:,2:end);    
    % create a Markov matrix
    function M = markovify(af)
        M = zeros(size(af));
        rsum = sum(af,2);
        for r = 1:length(rsum)
            M(r,:) = af(r,:)./rsum(r);
        end
    end

    function SM = smarkovify(af)
        d = diag(sum(af,2).^-.5);
        SM = d*af*d;
    end

    % kernel function
    function af = basicKernel(s)
        af = exp(-s.^2/epsilon^2);
    end

    function fancy = pKernel(s)
        af = basicKernel(s);
        sqrtp = sqrt(sum(af,2)); % column vector of sums across rows
        pmatrix = sqrtp * sqrtp';
        fancy = af./pmatrix;
    end
end