function [R,Phi] = R_Miranda_Aluvial(T,mu)
    Phi = 1 + 1/(12*T-mu*T) - 2/(5*T)*exp(-2*(log(T)-1/5)^2);
    R = max((mu-1)/Phi + 1,1);
end