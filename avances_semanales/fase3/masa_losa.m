function [Mlloc] = masa_losa(rho,espesor,Af, L)
%MASA_LOSA Matriz que se adiciona a la matriz de masa de un elemento de
%pórtico para considerar la masa del entrepiso
%   La masa se coloca en los grados de libertad traslacionales
Mlloc = rho*espesor*Af*L/2*...
    [   1   0   0   0   0   0
        0   0   0   0   0   0
        0   0   0   0   0   0
        0   0   0   1   0   0
        0   0   0   0   0   0
        0   0   0   0   0   0];
end

