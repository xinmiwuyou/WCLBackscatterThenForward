function [z,flag] = bisection1(a, b, NN, eps_step, eps_abs, nu, lambda)
% Find the unique solution $z_i^{*}$ of $(1-\zeta^* / bar{R}_i ) f(z_i) = \zeta^*$. For more details, 
% please refer to Theorem 1.

% a is the start point
% b is the end point
% NN is the maximal times of search
% eps_step is the search step
% eps_abs the accuracy of search
% flag is the indicator that shows if the parameter setting is feasible.

if abs(nu * fun1(a)- log(2)* lambda) <= 1e-5
    z = a;
    flag = 0;
    return;
elseif abs( nu * fun1(b)- log(2)* lambda) <= 1e-5
    z = b;
    flag = 0;
    return;
elseif ( nu * fun1(a)- log(2)* lambda) * (nu * fun1(b)- log(2)* lambda) > 0
    %     error('fun(a) and fun(b) dont have opposite signs');
    flag = 1; % the indicator shows that the parameter setting is not feasible.
    z = (a+b)/2; % just random assignment
    return;
end

for k= 1: NN
    c = (a+b)/2;
    % Check if we found a root or not
    % we should continue with:
    %          [a, c] if f(a) and f(c) have opposite signs, or
    %          [c, b] if f(c) and f(b) have opposite signs.
    if abs( nu * fun1(c)- log(2)* lambda) <= 1e-5
        z = c;
        flag = 0;
        return;
    elseif ( nu * fun1(c)- log(2)* lambda) * ( nu * fun1(a)- log(2)* lambda) < 0
        b = c;
    else
        a = c;
    end
    
    % If |b - a| < eps_step, check whether or not
    %       |f(a)| < |f(b)| and |f(a)| < eps_abs and return 'a', or
    %       |f(b)| < eps_abs and return 'b'.
    if (b-a < eps_step)
        if (abs( nu * fun1(a)- log(2)* lambda) < abs( nu * fun1(b)- log(2)* lambda) && abs( nu * fun1(a)- log(2)* lambda) < eps_abs )
            z = a;
            flag = 0;
            return;
        elseif (abs( nu * fun1(b)- log(2)* lambda) < eps_abs )
            z = b;
            flag = 0;
            return;
        end
    end
    
end
error('the method did not converge');
end

