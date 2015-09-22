n=100;

x_star = zeros(n,1);
nRun = 10000;

B = zeros(n-1,n);
for i=1:n-1
    B(i,i)=1;
    B(i,i+1)=-1;
end
lambda_star = zeros(n,1);
optimal_value = zeros(n,1);
tolerance = 1;

optimal_value(1)=1;
optimal_value(n)=n;

for K=1:2
    fprintf('Computing sparisity level %d\n', K);
    x_star(1:n)=1;
    x_star(n-K+1:2:n)=-1;
    musc = zeros(n-1,1);
    musc(n-K:2:end)=1;
    musc(n-K+1:2:end)=-1;
    s = 1:n-K-1;
    s_length = length(s);  
    cost_pre = 1000000;
    chance = 0;
    for point=1:400
        lambda = 0.01*point;
        fprintf('Computing point number %d: lambda is %f\n', point, lambda);
        avgJ = 0;
        g = normrnd(0,1,n,nRun);
        avg = zeros(nRun,1);
        a = zeros(s_length,nRun);
        parfor i=1:nRun
            [a(:,i),J] = gradientMethodII(g(:,i),B,lambda,s,musc);
            avgJ(i) = computeCost(g(:,i),B,lambda,s,musc,a(:,i));
        end
        cost_now = sum(avgJ)/nRun;
        fprintf('Expectation is %f...\n', cost_now);
        if point==1
            optimal_value(K+1) = cost_now;
            lambda_star(K+1) = lambda;
        else
            if cost_now < optimal_value(K+1)
                optimal_value(K+1) = cost_now;
                lambda_star(K+1) = lambda;
                chance=0;
            else
                chance = chance+1;
            end   
        end
        if chance>tolerance 
            break;
        end         
        cost_pre = cost_now;
    end  
end