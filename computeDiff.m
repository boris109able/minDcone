function[diff]=computeDiff(g,B,lambda,s,musc,tmus)
%% compute J^{(1)}(mu)
diff = -2*lambda*B(s,:)*(g-lambda*B'*musc-lambda*B(s,:)'*tmus);
end