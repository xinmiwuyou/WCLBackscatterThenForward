% Bin Lyu et al., "Backscatter then Forward: A Relaying Scheme for Batteryless IoT Networks", accepted by IEEE WCL.
% By Bin Lyu
% Dec. 21, 2019.
% Sum-rate versus distance between HAP and G_i.

clear;
% parameters setting
NoiseH = 1e-10; % -70 dBm the noise power at the HAP.
NoiseG = 1e-10; % -70 dBm the noise power at the relays.
N = 5; % the number of WSNs
eta = 0.7; % the energy efficiency
Count = 500; % the number of random channel generations
alpha = 1; % the refelction coefficient
kappa = 3; % path-loss exponent
tol = 1e-3; % tolerence for golden section search method.
PHdBm = 20; % dBm
PH = 1e-3 .* 10.^(PHdBm./10); % the transmit power of the HAP

Dug = 2; % the distance between each pair of gateway and user.
DigSet = 4:2:12; % the set of distances between the HAP and the Gateways





%% setting the number of antennas
M = 5; % the number of antennas

% the set to store all obtained results for each transmit tpower and each
% channel generations for the proposed scheme
RProSumM5 = zeros(Count,length(DigSet));
% the set to store all obtained results for each transmit tpower and each
% channel generations for the benchmark with random energy beamforming
REBSumM5 = zeros(Count,length(DigSet));
% the set to store all obtained results for each transmit tpower and each
% channel generations for the benchmark with equal time allocation
RETSumM5 = zeros(Count,length(DigSet));

% the set to store the average results with M=5 for the proposed scheme
RProAveM5 = zeros(1, length(DigSet));
% the set to store the average results with M=5 for the benchmark with random energy beamforming
REBAveM5 = zeros(1, length(DigSet));
% the set to store the average results with M=5 for the benchmark with equal time allocation
RETAveM5 = zeros(1,length(DigSet));


for ii = 1: length(DigSet)
    Dig = DigSet(ii);
    Diu = Dug + Dig; % the distance between the HAP and the Users
    
    % it is used to compute the number of that the bisection method fails or the CVX fails.
    countNumPro = 0; 
    countNumEB = 0;
    countNumET = 0;
    for i =1: Count
        
        
        %  random channel generation
        hIU = sqrt(Diu^(-kappa)/2) * ( randn(M,N) + 1i*randn(M,N)); % the downlink channels between the HAP and the users.
        hIG = sqrt(Dig^(-kappa)/2) * ( randn(M,N) + 1i*randn(M,N)); % the downlink channels between the HAP and the gateways.
        gIU = sqrt(Dug^(-kappa)/2) * ( randn(1,N) + 1i*randn(1,N)); % the uplink channels between the users and the gateways.
        gIG = sqrt(Dig^(-kappa)/2) * ( randn(M,N) + 1i*randn(M,N)); % the uplink channels between the gateways and the HAP.
        
        
        
        % the SNR as shown in the constraint C1 or C6
        gammaIG = zeros(N,1);
        for j =1:N
            % following Lemma 1.
            gammaIG(j) = PH .* abs(gIU(:,j))^2 .* norm(hIU(:,j))^2 .* alpha ./ NoiseG;
        end
        % the partial SNR shown in the constraint C2 or C7
        gammaIH = zeros(N,1);
        for j=1:N
            gammaIH(j) = norm(gIG(:,j))^2 ./ NoiseH;
        end
        
        %% proposed scheme with M =5
        % the achievable rate from each user to its gateway
        r1 = log2(1+gammaIG);
        
        %%% golden search method to find the optimal EH time
        % the setting search range
        mm = 0;
        nn = 1;
        cc = (-1+sqrt(5))/2;
        % the number of searching times.
        countGS = 0;
        % the loop of the golden search method.
        while 1
            
            countGS = countGS + 1;
            
            % x1 and x2 are used to limit the time range of the duration of EH phase,i.e., b_0
            x1 = nn - cc * (nn-mm);
            x2 = mm + cc * (nn-mm);
            
            % the function "BackThenForward" is the proposed algorithm (Algorithm 1) to solve P4, and obtain the optimal solution based on Eqs.(1)-(3).
            
            % R1 and R2 are the maximum achievable rate with the optimal solution
            % t1 and t2 denote the time slots for data backscattering
            % tau1 and tau2 denote the time slots for data forwarding
            % e1 and e2 denote the energy consumed for data forwarding at the gateways
            % W01 and W02 denotes the optimal energy beamforming during the EH phase
            % flagl  and flagu are the indicators to show whether the bisection method works
            [R1, t1,tau1,e1, W01,flagl] = BackThenForward(r1, eta, PH, hIG, gammaIH,N,M,x1);
            [R2, t2,tau2,e2, W02,flagu] = BackThenForward(r1, eta, PH, hIG, gammaIH,N,M,x2);
            
            % filag1 ==1 or flagu ==1, it indicates that the biection
            % method cannot find the optimal solution. So, we ignore the
            % results for this iteration (channel generation).
            if flagl ==1 || flagu ==1
                RProSumM5(i,ii) = 0;
                countNumPro = countNumPro +1; % to compute the number that the bisection method fails
                break; % ignore this iteration.
            end
            
            % update the golden search method.
            if R1 <= R2
                mm = x1;
            else
                nn = x2;
            end
            % the condition that the golden search method stops updating.
            if abs(nn-mm) < tol
                RProSumM5(i,ii) = (R1 + R2)/2; % the optimal achievable rate for each iteration under the given transmit power for M=5;
                break;
            end
        end
        
        
        %% the benchmark with random energy beamforming for M =5
        
        % random energy beamforming design
        temp0 = 0;
        for j=1:N
            temp0 = temp0 + hIG(:,j);
        end
        W0Ben = temp0 / N; % the random energy beamforming for the EH phase
        
        % With the random energy beamforming design, we only need to optimize the optimal time allocation
        
        % the CVX (http://cvxr.com/cvx/) is used to compute the optimal time allocation directly.
        cvx_begin quiet
        %%%% define the optimization variables %%%%
        variable t(N,1) nonnegative % the duration of each user for data backscattering
        variable b nonnegative % the duration of the EH phase
        variable tau(N,1) nonnegative % the duration of each user for data forwarding
        variable e(N,1) nonnegative % the energy consumed for data forwarding
        variable R(N,1) nonnegative % the achievable rate of all users
        
        %%% solve the optimization problem, which is similar to P3 of this paper %%%%%%%%
        maximize sum(R) % the objective function 
        % the constraints
        subject to 
        b + sum(t) + sum(tau) <=1;
        R - t .* log2(1 + gammaIG) <= 0;
        R - (-rel_entr(tau, tau + gammaIH .* e )  / log(2)) <= 0;
        for j=1:N
            e(j) - eta .* PH * abs(hIG(:,j)' * W0Ben )^2 .* b <= 0;
        end
        t - ones(N,1) <= 0;
        tau - ones(N,1) <= 0;
        b - 1 <= 0;
        cvx_end
        
        % if CVX cannot find the optimal solution, we ignore this iteration
        if isnan(cvx_optval) == 1
            countNumEB = countNumEB+1; % the total times that CVX fails
            REBSumM5(i,ii) = 0; % if CVX fails, set the achievable rate to zero. Note that this setting will not affect the final results
        else 
            % the CVX can find the optimal solution, store this solution.
            REBSumM5(i,ii) = sum(R);
        end
        
        
        %% the benchmark with the equal time allocation for M=5
        
        % the equal time allocation
        time = 1/(2*N+1);
        b = time;
        t = ones(N,1) * time;
        tau = ones(N,1) * time;
        
        % With the equal time allocation, we only need to optimize the energy beamforming vector (matrix) during the EH phase
        
        cvx_begin quiet
        
        %%%% define the optimization variables %%%%
        variable W(M,M) complex semidefinite % the energy beamforming matrix during the EH phase
        variable e(N,1) nonnegative % the energy comsumed for data forwarding
        variable R(N,1) nonnegative % the achievable rate
        
       %%% solve the optimization problem, which is similar to P3 of this paper %%%%%%%%
        maximize sum(R) % the objective function
        % the constraints
        subject to
        R - t .* log2(1 + gammaIG) <= 0;
        R - (-rel_entr(tau, tau + gammaIH .* e )  / log(2)) <= 0;
        for j=1:N
            e(j) - eta .* PH * real(trace(hIG(:,j) * hIG(:,j)' * W)) <= 0;
        end
        real(trace(W)) - b <= 0 ;
        
        cvx_end
        
       % if CVX cannot find the optimal solution, we ignore this iteration
        if isnan(cvx_optval) == 1
            countNumET = countNumET+1; % the total times that CVX fails
            RETSumM5(i,ii) = 0; % if CVX fails, set the achievable rate to zero. Note that this setting will not affect the final results
        else
        % the CVX can find the optimal solution, store this solution.
            RETSumM5(i,ii) = sum(R);
        end
        
        
    end
    % the average achievable rate with M=5 for the proposed scheme
    RProAveM5(1,ii) = sum(RProSumM5(:,ii) ) / (Count- countNumPro);
    % the average achievable rate with M=5 for the benchmark with random energy beamforming
    REBAveM5 (1,ii) = sum(REBSumM5(:,ii)) / (Count-countNumEB);
    % the average achievable rate with M=5 for the benchmark with equal time allocation
    RETAveM5(ii) = sum(RETSumM5(:,ii)) / (Count-countNumET);
end


%% setting the number of antennas
M = 10; % the number of antennas


% the set to store all obtained results for each transmit tpower and each
% channel generations for the proposed scheme
RProSumM10 = zeros(Count,length(DigSet));
% the set to store all obtained results for each transmit tpower and each
% channel generations for the benchmark with random energy beamforming
REBSumM10 = zeros(Count,length(DigSet));
% the set to store all obtained results for each transmit tpower and each
% channel generations for the benchmark with equal time allocation
RETSumM10 = zeros(Count,length(DigSet));

% the set to store the average results with M=5 for the proposed scheme
RProAveM10 = zeros(1, length(DigSet));
% the set to store the average results with M=5 for the benchmark with random energy beamforming
REBAveM10 = zeros(1, length(DigSet));
% the set to store the average results with M=5 for the benchmark with equal time allocation
RETAveM10 = zeros(1,length(DigSet));

for ii = 1: length(DigSet)
    Dig = DigSet(ii);
    Diu = Dug + Dig; % the distance between the HAP and the Users
    
    % it is used to compute the number of that the bisection method fails or the CVX fails.
    countNumPro = 0; 
    countNumEB = 0;
    countNumET = 0;
    for i =1: Count

        % the channel generation
        hIU = sqrt(Diu^(-kappa)/2) * ( randn(M,N) + 1i*randn(M,N)); % the downlink channels between the HAP and the users.
        hIG = sqrt(Dig^(-kappa)/2) * ( randn(M,N) + 1i*randn(M,N)); % the downlink channels between the HAP and the gateways.
        gIU = sqrt(Dug^(-kappa)/2) * ( randn(1,N) + 1i*randn(1,N)); % the uplink channels between the users and the gateways.
        gIG = sqrt(Dig^(-kappa)/2) * ( randn(M,N) + 1i*randn(M,N)); % the uplink channels between the gateways and the HAP.
        
        
        
        % the SNR as shown in the constraint C1 or C6
        gammaIG = zeros(N,1);
        for j =1:N
            % following Lemma 1.
            gammaIG(j) = PH .* abs(gIU(:,j))^2 .* norm(hIU(:,j))^2 .* alpha ./ NoiseG;
        end
        % the partial SNR shown in the constraint C2 or C7
        gammaIH = zeros(N,1);
        for j=1:N
            gammaIH(j) = norm(gIG(:,j))^2 ./ NoiseH;
        end
        
        %% proposed scheme with M =10
        % the achievable rate from each user to its gateway
        r1 = log2(1+gammaIG);
        
        %%% golden search method to find the optimal EH time
        % the setting search range
        mm = 0;
        nn = 1;
        cc = (-1+sqrt(5))/2;
        % the number of searching times.
        countGS = 0;
        % the loop of the golden search method.
        while 1
            
            countGS = countGS + 1;
            
            % x1 and x2 are used to limit the time range of the duration of EH phase,i.e., b_0
            x1 = nn - cc * (nn-mm);
            x2 = mm + cc * (nn-mm);
            
            % the function "BackThenForward" is the proposed algorithm (Algorithm 1) to solve P4, and obtain the optimal solution based on Eqs.(1)-(3).
            
            % R1 and R2 are the maximum achievable rate with the optimal solution
            % t1 and t2 denote the time slots for data backscattering
            % tau1 and tau2 denote the time slots for data forwarding
            % e1 and e2 denote the energy consumed for data forwarding at the gateways
            % W01 and W02 denotes the optimal energy beamforming during the EH phase
            % flagl  and flagu are the indicators to show whether the bisection method works
            [R1, t1,tau1,e1, W01,flagl] = BackThenForward(r1, eta, PH, hIG, gammaIH,N,M,x1);
            [R2, t2,tau2,e2, W02,flagu] = BackThenForward(r1, eta, PH, hIG, gammaIH,N,M,x2);
            
            % filag1 ==1 or flagu ==1, it indicates that the biection
            % method cannot find the optimal solution. So, we ignore the
            % results for this iteration (channel generation).
            if flagl ==1 || flagu ==1
                RProSumM10(i,ii) = 0;
                countNumPro = countNumPro +1; % to compute the number that the bisection method fails
                break; % ignore this iteration.
            end
            
            % update the golden search method.
            if R1 <= R2
                mm = x1;
            else
                nn = x2;
            end
            % the condition that the golden search method stops updating.
            if abs(nn-mm) < tol
                RProSumM10(i,ii) = (R1 + R2)/2; % the optimal achievable rate for each iteration under the given transmit power for M=5;
                break;
            end
        end
        
        
        %% the benchmark with random energy beamforming for M =10
        
      % random energy beamforming design
        temp0 = 0;
        for j=1:N
            temp0 = temp0 + hIG(:,j);
        end
        W0Ben = temp0 / N; % the random energy beamforming for the EH phase
        
        % With the random energy beamforming design, we only need to optimize the optimal time allocation
        
        % the CVX (http://cvxr.com/cvx/) is used to compute the optimal time allocation directly.
        cvx_begin quiet
        %%%% define the optimization variables %%%%
        variable t(N,1) nonnegative % the duration of each user for data backscattering
        variable b nonnegative % the duration of the EH phase
        variable tau(N,1) nonnegative % the duration of each user for data forwarding
        variable e(N,1) nonnegative % the energy consumed for data forwarding
        variable R(N,1) nonnegative % the achievable rate of all users
        
        %%% solve the optimization problem, which is similar to P3 of this paper %%%%%%%%
        maximize sum(R) % the objective function 
        % the constraints
        subject to 
        b + sum(t) + sum(tau) <=1;
        R - t .* log2(1 + gammaIG) <= 0;
        R - (-rel_entr(tau, tau + gammaIH .* e )  / log(2)) <= 0;
        for j=1:N
            e(j) - eta .* PH * abs(hIG(:,j)' * W0Ben )^2 .* b <= 0;
        end
        t - ones(N,1) <= 0;
        tau - ones(N,1) <= 0;
        b - 1 <= 0;
        cvx_end
        
        % if CVX cannot find the optimal solution, we ignore this iteration
        if isnan(cvx_optval) == 1
            countNumEB = countNumEB+1; % the total times that CVX fails
            REBSumM10(i,ii) = 0; % if CVX fails, set the achievable rate to zero. Note that this setting will not affect the final results
        else 
            % the CVX can find the optimal solution, store this solution.
            REBSumM10(i,ii) = sum(R);
        end
        
        
  %% the benchmark with the equal time allocation for M=10
        
        % the equal time allocation
        time = 1/(2*N+1);
        b = time;
        t = ones(N,1) * time;
        tau = ones(N,1) * time;
        
        % With the equal time allocation, we only need to optimize the energy beamforming vector (matrix) during the EH phase
        
        cvx_begin quiet
        
        %%%% define the optimization variables %%%%
        variable W(M,M) complex semidefinite % the energy beamforming matrix during the EH phase
        variable e(N,1) nonnegative % the energy comsumed for data forwarding
        variable R(N,1) nonnegative % the achievable rate
        
       %%% solve the optimization problem, which is similar to P3 of this paper %%%%%%%%
        maximize sum(R) % the objective function
        % the constraints
        subject to
        R - t .* log2(1 + gammaIG) <= 0;
        R - (-rel_entr(tau, tau + gammaIH .* e )  / log(2)) <= 0;
        for j=1:N
            e(j) - eta .* PH * real(trace(hIG(:,j) * hIG(:,j)' * W)) <= 0;
        end
        real(trace(W)) - b <= 0 ;
        
        cvx_end
        
       % if CVX cannot find the optimal solution, we ignore this iteration
        if isnan(cvx_optval) == 1
            countNumET = countNumET+1; % the total times that CVX fails
            RETSumM10(i,ii) = 0; % if CVX fails, set the achievable rate to zero. Note that this setting will not affect the final results
        else
        % the CVX can find the optimal solution, store this solution.
            RETSumM10(i,ii) = sum(R);
        end
        
        
    end
    % the average achievable rate with M=10 for the proposed scheme
    RProAveM10(1,ii) = sum(RProSumM10(:,ii) ) / (Count- countNumPro);
    % the average achievable rate with M=10 for the benchmark with random energy beamforming
    REBAveM10 (1,ii) = sum(REBSumM10(:,ii)) / (Count-countNumEB);
    % the average achievable rate with M=10 for the benchmark with equal time allocation
    RETAveM10(ii) = sum(RETSumM10(:,ii)) / (Count-countNumET);
end


plot(DigSet, RProAveM5,'r-s','LineWidth',2);
hold on;
plot(DigSet, REBAveM5, 'k:s','LineWidth',2);
plot(DigSet, RETAveM5, 'b-.s','LineWidth',2);
plot(DigSet, RProAveM10,'r-o','LineWidth',2);
plot(DigSet, REBAveM10 , 'k:o','LineWidth',2);
plot(DigSet, RETAveM10, 'b-.o','LineWidth',2);

xlabel('Distance between HAP and G_i (m)');
ylabel('Achievable sum-rate (bits/s/Hz)');
legend('Proposed scheme, M=5', 'Random EM designs, M=5', 'Equal Time Allocation, M=5', 'Proposed scheme, M=10', 'Random EM designs, M=10', 'Equal Time Allocation, M=10');








