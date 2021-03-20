function [R, t,tau,e, W0,flag] = BackThenForward(r1, eta, PH, hIG, gamma2,N,M,b )
% The proposed algorithm to solve P4. For more details, please refer to Algorithm P4 (Theorem 1 and Appendix A).

% R is the maximum achievable rate with the optimal solution
% t denotes the time slots for data backscattering
% tau denotes the time slots for data forwarding
% e denotes the energy consumed for data forwarding at the gateways
% W0 denotes the optimal energy beamforming during the EH phase
% flag is the indicator to show whether the bisection method works
%% Initialization
%%%  Initialize Lagranian multipliers
mu = 1e3 * rand(N,1);
% the range of zeta as shown in Theorem 1
zeta = min(r1) / 2;  % 0 < zeta < r1;

% compute A as shown in Theorem 1
A = zeros(M,M);
for j=1:N
    A = A + eta * PH * mu(j) * hIG(:,j) * hIG(:,j)';
end

for j =1:M
    A(j,j) = real(A(j,j));
end

% eigenvalue decomposition of A as shown in Appendix A
[ UB, SB, VB ] = svd( A );
% find the maximum eigenvalue value
[ lambda, idx_max ] = max( diag( SB ) );

rho = lambda; % rho^* = \lambad_{A,1}

% u_{A,1} is the unit-norm eigenvector
w0 = UB(:,idx_max); % the energy beamforming desgin during b0;
% refer to Eq. (1) of this paper
W0 = b * w0 * w0';

for j=1:M
    W0(j,j) = real(W0(j,j));
end

% the energy allocation
e = zeros(N,1);
for j = 1:N
    % refer to Eq. (2) of this paper
    e(j) = eta * PH * real(trace(hIG(:,j) * hIG(:,j)' * W0));
end


% the parameters for the bisetion method
StartPo = 0;  % the minimum of \xi is 1;
EndPo = 1e100;  % the maximum of \xi is the postive infinite
NN = 1e3; % the maximum number of iterations
eps_step = 1e-5; % the step size
eps_abs = 1e-5; % the toleration

z = zeros(N,1);
ww = 1 - zeta./r1;
for k=1:N
    % the bisection method to find the optimal z, where z_i^* is the unique solution of (1- \zeta/R_i) f(z_i) = \zeta
    % for more details, please refer to Theorem 1.
    [z(k),flag] = bisection1(StartPo, EndPo, NN, eps_step, eps_abs, ww(k), zeta);
    % if flag is equal to one, the bisection method cannot work. Then, we ignore the results obtain for this time of channel generation.
    if flag == 1
        t = zeros(N,1);
        tau = zeros(N,1);
        R = 0;
        return;
    end
end

% refer to Eq. (3)
tau = min( e.* gamma2 ./ z ,1-b);
% refer to Appendix A, B is defined below Eq. (8)
B = A - rho * diag(ones(M,1));


%%% compute the dual function above Theorem 1
gMin = -sum(mu .* e) + rho * b - zeta * b + zeta - zeta * sum(tau) + trace(B*W0);
for j = 1:N
    if tau(j) == 0
        gMin = gMin + 0;
    else
        gMin = gMin + ww(j) .* tau(j) * log2(1+ gamma2(j) .* e(j) ./ tau(j));
    end
end

% update Lagrangian multipliers
mu_opt = mu;
zeta_opt = zeta;


%% Main Loop of Algorithm 1
% this loop is used to find the optimal Lagrangian multipliers
iter =1;
CountNum =0;
Iter = 200; % the maximum number of iterations
while iter <= Iter
    
    % sub-gradient of Lagrangian Multipliers
    Temp = 0;
    for j=1:N
        if tau(j) == 0
            Temp = Temp + 0;
        else
            Temp = Temp + 1/r1(j) * tau(j) * log2(1+ gamma2(j) * e(j)/ tau(j));
        end
    end
    zetaSG = 1- b - sum(tau) - Temp ;
    muSG = zeros(N,1);
    for j=1:N
        muSG(j) = eta * PH * real(trace(hIG(:,j) * hIG(:,j)' * W0)) - e(j);
    end
    
    
    % step size
    % note that the step size can be further designed to acheive better peformance of the convergence speed.
    stepZeta = zeta /20;
    stepMu = mu /20;
    
    % update Lagrangian multipliers
    muNew = max(1e-5, mu - stepMu .* muSG); % lambda_i > 0;
    zetaNew = min( max(0, zeta - stepZeta * zetaSG), min(r1));
    
    % compute A as shown in Theorem 1
    A = zeros(M,M);
    for j =1:N
        A = A + eta * PH * muNew(j) * hIG(:,j) * hIG(:,j)';
    end
    for j =1:M
        A(j,j) = real(A(j,j));
    end
    
    % eigenvalue decomposition of A as shown in Appendix A
    [ UB, SB, VB ] = svd( A );
    % find the maximum eigenvalue value
    [ lambda, idx_max ] = max( diag( SB ) );
    
    rhoNew = lambda;  % rho^* = \lambad_{A,1}
    
    % u_{A,1} is the unit-norm eigenvector
    w0 = UB(:,idx_max); % the energy beamforming desgin during b0;
    % refer to Eq. (1) of this paper
    W0 = b * w0 * w0';
    for j=1:M
        W0(j,j) = real(W0(j,j));
    end
    
    
    
    z = zeros(N,1);
    wwNew = 1 - zetaNew./r1;
    for k=1:N
        % the bisection method to find the optimal z, where z_i^* is the unique solution of (1- \zeta/R_i) f(z_i) = \zeta
        % for more details, please refer to Theorem 1.
        [z(k),flag] = bisection1(StartPo, EndPo, NN, eps_step, eps_abs, wwNew(k), zetaNew);
        % if flag is equal to one, the bisection method cannot work. Then, we ignore the results obtain for this time of channel generation.
        if flag == 1
            t = zeros(N,1);
            tau = zeros(N,1);
            R = 0;
            return;
        end
    end
    
    % the energy allocation
    e = zeros(N,1);
    for j = 1:N
        % refer to Eq. (2)
        e(j) = eta * PH * real(trace(hIG(:,j) * hIG(:,j)' * W0));
    end
    
    % refer to Eq. (3)
    tau = min( e.* gamma2 ./ z ,1-b);
    % refer to Appendix A, B is defined below Eq. (8)
    B = A - rhoNew * diag(ones(M,1));
    
    %%% compute the dual function above Theorem 1
    gTemp = -sum(muNew .* e) + rhoNew * b - zetaNew * b + zetaNew - zetaNew * sum(tau) + trace(B*W0);
    for j = 1:N
        if tau(j) == 0
            gTemp = gTemp + 0;
        else
            gTemp = gTemp + wwNew(j) .* tau(j) * log2(1+ gamma2(j) .* e(j) ./ tau(j));
        end
    end
    
    % solving the dual problem
    gminTemp = gMin;
    if gminTemp <= gTemp
        gMin = gminTemp;
        
    else
        gMin = gTemp;
        mu_opt = muNew;
        zeta_opt = zetaNew;
        CountNum = CountNum +1;
    end
    iter = iter +1;
    
    % the condition of convergence
    if ( iter >Iter )
        iter = 1e10; % to break the itration
    else
        % update the Lagrangian multipliers
        mu = muNew;
        zeta = zetaNew;
    end
    
end

%% Compute the optimal results with the optimal Lagrangian multipliers

% compute A as shown in Theorem 1
A = zeros(M,M);
for j=1:N
    A = A + eta * PH * mu_opt(j) * hIG(:,j) * hIG(:,j)';
end

for j =1:M
    A(j,j) = real(A(j,j));
end

% eigenvalue decomposition of A as shown in Appendix A
[ UB, SB, VB ] = svd( A );
% find the maximum eigenvalue value
[ lambda, idx_max ] = max( diag( SB ) );


% u_{A,1} is the unit-norm eigenvector
w0 = UB(:,idx_max); % the energy beamforming desgin during t0;

% refer to Eq. (1) of this paper
W0 = b * w0 * w0';
for j=1:M
    W0(j,j) = real(W0(j,j));
end

z = zeros(N,1);
ww = 1 - zeta_opt./r1;
for k=1:N
    % The bisection method is used to find the optimal z, where z_i^* is the unique
    % solution of $(1- \zeta^{*} / bar{R}_i) f(z_i) = \zeta^{*}$.
    % For more details, please refer to Theorem 1.
    z(k) = bisection1(StartPo, EndPo, NN, eps_step, eps_abs, ww(k), zeta_opt);
end

% the energy allocation
e = zeros(N,1);
for j = 1:N
    % refer to Eq. (2)
    e(j) = eta * PH * real(trace(hIG(:,j) * hIG(:,j)' * W0));
end

% refer to Eq. (3)
tau = min( e.* gamma2 ./ z ,1-b);

t =zeros(N,1);
for j=1:N
    if tau(j) == 0
        t(j) = 0;
    else
        % compute the time for data backscattering based on the equation
        % above (P3)
        t(j) = 1./r1(j) * tau(j) * log2(1+ gamma2(j) * e(j) / tau(j) );
    end
end

% compute the achievable rate
R = 0;
for j =1:N
    if tau(j) == 0
        R = R + 0;
    else
        R = R + tau(j) * log2(1+ gamma2(j) * e(j) / tau(j));
    end
end













