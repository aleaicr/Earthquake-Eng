function [Uppg,T] = T_Fourier(dt,uppg)
% T_fourier aplica la transofrmada de fourier discreta a un registro uddg con
% paso temporal dt, y da como resultado F(uddg) y los periodos asociados
    n = length(uppg);                                                       % Cantidad de datos
    Uppg = fft(uppg);                                                       % Transformada de fourier discreta
    Fs = 1/dt;
    s = -Fs/2:Fs/n:(Fs/2-Fs/n);
    for i = 1:length(s)
        T(i,1) = 1/s(1,i);
    end
end