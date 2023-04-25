function [mu,sigma] = campbell_1997_AH(Mw,r,F,S_SR,S_HR)
% Obtiene mu_ln y sigma_ln de AH utilizando Campbell 1997


% GMPE para fuente A
% media lognormal
mu = -3.512 + 0.904*Mw - 1.328*log(sqrt(r + (0.149*exp(0.647*Mw))^2)) + F*(1.125 - 0.112*log(r) -0.0957*Mw) + S_SR*(0.440 -0.171*log(r)) + S_HR*(0.405 -0.222*log(r));

% AH
AH = exp(mu);

% sigma lognormal
if AH < 0.068
    sigma = 0.55;
elseif and(AH <= 0.21, AH>=0.068)
    sigma = 0.173 - 0.140*log(AH);
elseif AH > 0.21
    sigma = 0.39;
end


end