function[cost]=computeCost(g,B,lambda,s,musc,tmus)
%% compute J(mu)
cost = sum((g-lambda*B'*musc-lambda*B(s,:)'*tmus).^2);
end