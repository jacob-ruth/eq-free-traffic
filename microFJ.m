function F = microFJ(sys1,sys2,M,N,L, v0, tau)
%% parameters and variables
% N - discretization points
% M - number of cars (M must divide N)
% L - length of road
% tau - momentum constant;
k = N/M;
mu = sys1(end);
c = sys1(end - 1);     % wavespeed

u = sys1(1:N);
u0 = sys2(1:N);

%% operators & differentiation matrices

D  = fourdif(N,1)*2*pi/M;
D2 = fourdif(N,2)*(2*pi/M)^2;

e0    = ones(N,1);
shift = sparse([[N-k+1:N] [1:(N-k)]],1:N,e0,N,N);                               % shift matrix
e = sparse(1:N,[1:N],e0,N,N);                                                   % identity

z = sparse(N,N);

LN = -c^2*tau*D2+c*D;                              % linear comp. of DE
LD = D;                                                                %  phase condition matrix
ln_consv = sparse(ones(M,1),k*[1:M],ones(M,1),1,N);     % conservation law vector

%% Function
F = [ LN*u + optimalVelocity(shift*u, v0) - optimalVelocity(u, v0) + mu.*ones(N,1);
    ln_consv*u - L; ...
    (LD*u)'*(u-u0)...
    ];
%% Jacobian computation
%{
if nargout > 1
    J0 = diag(optimalVHway(shift*u,v0));
    J1 = diag(optimalVHway(u,v0));
    Jc = (-2*c*tau*D2+D)*u;
    Jmu = ones(N,1);

    J  = sparse(...
        [LN+J0-J1, Jc, Jmu;...
        ln_consv, 0, 0;...
        (LD*u)'+ u'*LD - u0'*LD, 0, 0]...
        );
end
%}

%% optimal veloctiy function
    function v = optimalVelocity(headway,v0)
        h = 1.2;
        v = v0 * (tanh(headway - h) + tanh(h));
    end

%% derivative of optimal velocity function
    function v = optimalVHway(headway,v0)
        h = 1.2;
        v = v0 * (sech(headway - h).^2);
    end
end

