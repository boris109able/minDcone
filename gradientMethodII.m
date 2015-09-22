function[mu, J] = gradientMethodII(g,B,lambda,s,musc)
%% compute the minimum value of J in the set Q
%% minimize \|g-lambda*B'*musc-lambda*B_{s,Omega}'*\tilde{mu_s}\|_2^2
%% such that \|\tilde{mu_s}\|_{\infty}\leq 1
%% return the minima mu

s_length = numel(s);
n = size(B,2);
a_old = zeros(s_length,1);
alpha_old = 0.5;
b_old = a_old;
q = 0.5*(1-cos(pi/n));
L = 8*lambda^2;
J_old = 100000000;
J_new = computeCost(g,B,lambda,s,musc,a_old);
J(1) = J_new;
eps = 0.0001;
while(abs(J_new-J_old)>eps)
    a_new = threshold(b_old-computeDiff(g,B,lambda,s,musc,b_old)/L,1);
    eqn = [1 alpha_old^2-q -alpha_old^2];
    r = roots(eqn);
    r = r(r>0);
    alpha_new = r(1);
    beta = alpha_old*(1-alpha_old)/(alpha_old^2+alpha_new);
    b_new = a_new+beta*(a_new-a_old);
    J_old = J_new;
    J_new = computeCost(g,B,lambda,s,musc,a_new);
    J(length(J)+1)=J_new;
    a_old = a_new;
    b_old = b_new;
    alpha_old = alpha_new;
end
mu = a_old;
end