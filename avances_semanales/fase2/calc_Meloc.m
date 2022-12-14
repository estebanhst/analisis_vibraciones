function Meloc = calc_Meloc(tipo, L, A, I, rho, masa)
% Gran parte de esta implementación se deriva del trabajo del profesor 
% Diego Andrés Álvarez Marín https://github.com/diegoandresalvarez/elementosfinitos/blob/master/codigo/vigas/Timoshenko_python/K_exacta_portico_TE.ipynb
% Calcula la matriz de masa consistente M de un EF de pórtico o de cercha en 
% coordenadas locales. 
% LA MATRIZ DE MASA QUEDA EN kg si rho está en kg/m³ y lo demás en m, m²
% PARAMETROS:
% tipo    = 'EE' (pórtico)
%         = 'RR' (cercha)
%         = 'RE' (rotula/empotrado)
%         = 'ER' (empotrado/rotula)
% A      área  
% I      inercia
% L      longitud de la barra    
L2  = L^2;
% La matriz de masa consistente M se separa en dos clases de matrices
% M = M_rhoA + M_rhoI
% M_rhoA es la asociada a la Inercia Traslacional
% M_rhoI es la asociada a la Inercia Rotacional
% matriz de rigidez local expresada en el sistema de coordenadas locales
if strcmp(masa, 'condensada')
    Meloc = rho*A*L/2*...
        [   1   0   0   0   0   0
            0   1   0   0   0   0
            0   0   0   0   0   0
            0   0   0   1   0   0
            0   0   0   0   1   0
            0   0   0   0   0   0];
elseif strcmp(masa, 'consistente')
    switch tipo
        case 'EE'
            M_rhoA = rho*A*L/420*...
                [   140     0       0       70      0       0
                    0       156     22*L    0       54      -13*L
                    0       22*L    4*L2    0       13*L    -3*L2
                    70      0       0       140     0       0
                    0       54      13*L    0       156     -22*L
                    0       -13*L   -3*L2   0       -22*L   4*L2];
            M_rhoI = rho*I/30/L*...
                [   0       0       0       0       0       0
                    0       36      3*L     0       -36     3*L
                    0       3*L     4*L2    0       -3*L    -L2
                    0       0       0       0       0       0
                    0       -36     -3*L    0       36      -3*L
                    0       3*L     -L2     0       -3*L    4*L2];
        case 'RR'
            M_rhoA = rho*A*L/420*...
                [   140     0       0       70      0       0
                    0       140     0       0       70      0
                    0       0       0       0       0       0
                    70      0       0       140     0       0
                    0       70      0       0       140     0
                    0       0       0       0       0       0];
            M_rhoI = rho*I/L*...
                [   0       0       0       0       0       0
                    0       1       0       0       -1      0
                    0       0       0       0       0       0
                    0       0       0       0       0       0
                    0       -1      0       0       1       0
                    0       0       0       0       0       0];
        case 'ER'
            M_rhoA = rho*A*L/840*...
                [   280     0       0       140     0       0
                    0       408     72*L    0       117     0
                    0       72*L    16*L2   0       33*L    0
                    140     0       0       280     0       0
                    0       117     33*L    0       198     0
                    0       0       0       0       0       0];
            M_rhoI = rho*I/30/L*...
                [   0       0       0       0       0       0
                    0       36      6*L     0       -36     0
                    0       6*L     6*L2    0       -6*L    0
                    0       0       0       0       0       0
                    0       -36     -6*L    0       36      0
                    0       0       0       0       0       0];
        case 'RE'
            M_rhoA = rho*A*L/840*...
                [   280     0       0       140     0       0
                    0       198     0       0       117     -33*L
                    0       0       0       0       0       0
                    140     0       0       280     0       0
                    0       117     0       0       408     -72*L2
                    0       -33*L   0       0       -72*L   16*L2];
            M_rhoI = rho*I/30/L*...
                [   0       0       0       0       0       0
                    0       36      0       0       -36     6*L
                    0       0       0       0       0       0
                    0       0       0       0       0       0
                    0       -36     0       0       36      -6*L
                    0       6*L     0       0       -6*L    6*L2];
        otherwise
            error('Tipo de elemento no soportado.')
    end
    Meloc = M_rhoA+M_rhoI;
    %Meloc = M_rhoA;
end
end