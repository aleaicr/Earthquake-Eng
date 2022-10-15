function [CDF,PDF] = gutenberg_richter_acotada(m,beta,m0,mmax)
% Ley de recurrencia de Gutenberg-Richter Acotada (por arriba y abajo)

CDF = (1-exp(-beta*(m-m0)))/(1-exp(-beta*(mmax-m0)));
PDF = (beta*exp(-beta*(m-m0)))/(1-exp(-beta*(mmax-m0)));

end