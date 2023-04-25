function[prob] = probs(PGA,mu_ln_PGA,sigma_ln_PGA)
% Una vez obtenidas mu_ln y sigma_ln se puede calcular la probabilidad de
% excedencia como 1 - P_NoExcedencia
prob = 1 - normcdf(log(PGA),mu_ln_PGA,sigma_ln_PGA);
end