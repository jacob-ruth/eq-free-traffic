
function [F,J] = my_sys0(v,v0,par)

% system:
%
% c^2*tau*d" -c*d' = tanh(d(x+1) - h(x+1)) - tanh(d(x) - h(x)) %headways
% -c*alpha*h'      = -(h - c*beta*d') % can ignore
% v = [d; h; c]

%% parameters and variables 

N = par.N;          % discretization points
M = par.M;          % number of cars (M must divide N)   
L = par.L;          % length of road
vel = par.vel;      % v0 velocity
tau = par.tau;
k = N/M;

mu = v(end);
c = v(end-1);     % wavespeed

u = v(1:N);
u0 = v0(1:N);

% w = v(2*N+1:end-3);
% w0 = v0(2*N+1:end-3);

%% operators & differentiation matrices

D  = fourdif(N,1)*2*pi/M;
D2 = fourdif(N,2)*(2*pi/M)^2;

e0    = ones(N,1);
shift = sparse([[N-k+1:N] [1:(N-k)]],1:N,e0,N,N);                               % shift matrix
e = sparse(1:N,[1:N],e0,N,N);                                                   % identity

z = sparse(N,N);

LN = c^2*tau*D2-c*D;                              % linear comp. of DE
LD = D;                                                                %  phase condition matrix
ln_consv = sparse(ones(M,1),k*[1:M],ones(M,1),1,N);     % conservation law vector

e2     = [e,-e];
shift2 = [shift, -shift];

%% Function
  F = [ LN*u + optimalVelocity(shift*u, vel) - optimalVelocity(u, vel) ; ...
%     - [zeros(N,1); sparse(N,1)] ...
%     + [tanh(e2*u); sparse(N,1)] ...
%     - [tanh(shift2*u); sparse(N,1)] + mu*[ones(N,1); sparse(N,1)]; ...
    ln_consv*u - L; ... % - L?
    (LD*u)'*(u-u0)...
    ];

%% Jacobian computation
if nargout > 1
    J0 = diag(optimalVHway(shift*u,vel));
    J1 = diag(optimalVHway(u,vel));
    Jc = (-2*c*tau*D2+D)*u;
    Jmu = ones(N,1);
    
    J  = sparse(...
        [LN+J0-J1, Jc, Jmu;...
        ln_consv, 0, 0;...
        (LD*u)'+ u'*LD - u0'*LD, 0, 0]...
        );
end

    % optimal veloctiy function
    function v = optimalVelocity(headway,v0)
        h = 1.2;
        v = v0 * (tanh(headway - h) + tanh(h));
    end

    function v = optimalVHway(headway,v0)
        h = 1.2;
        v = v0 * (sech(headway - h).^2);
    end
end

