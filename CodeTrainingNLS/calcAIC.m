function [logL, AIC, SSE, varRes] = calcAIC(res, k)
    N = length(res); 
%     N = 10000;
    SSE = sum(res(1:N)'*res(1:N)); MSE = SSE/N;
    varRes = var(res(1:N));
    logL = -0.5*SSE/varRes - (N/2)*log(varRes);
%     logSS = log(MSE);
%     AIC = N*log(MSE) + 2*k;
    AIC = -2*logL + log(N)*k;
end