function Keloc = calc_Keloc(tipo, L, A, E, I)
% Calcula la matriz de rigidez local de un EF de pórtico o de cercha en 
% coordenadas locales
% PARAMETROS:
% tipo    = 'EE' (pórtico)
%         = 'RR' (cercha)
%         = 'RE' (rotula/empotrado)
%         = 'ER' (empotrado/rotula)
% A      área  
% E      módulo de elasticidad
% I      inercia
% L      longitud de la barra
    
AE  = A*E;       
EI  = E*I;       
L2  = L^2;
L3  = L^3;

% matriz de rigidez local expresada en el sistema de coordenadas locales
switch tipo
    case 'EE'
        Keloc = [   AE/L    0           0           -AE/L   0           0     
                    0       12*EI/L3    6*EI/L2     0       -12*EI/L3   6*EI/L2
                    0       6*EI/L2     4*EI/L      0       -6*EI/L2    2*EI/L 
                    -AE/L   0           0           AE/L    0           0      
                    0       -12*EI/L3   -6*EI/L2    0       12*EI/L3    -6*EI/L2
                    0       6*EI/L2     2*EI/L      0       -6*EI/L2    4*EI/L ];
    case 'RR'
        k=AE/L;
        Keloc = [   k   0   0  -k   0   0     
                    0   0   0   0   0   0
                    0   0   0   0   0   0
                   -k   0   0   k   0   0      
                    0   0   0   0   0   0
                    0   0   0   0   0   0 ];
    case 'ER'
        Keloc = [   AE/L    0           0           -AE/L   0           0     
                    0       3*EI/L3     3*EI/L2     0       -3*EI/L3    0
                    0       3*EI/L2     3*EI/L      0       -3*EI/L2    0 
                    -AE/L   0           0           AE/L    0           0      
                    0       -3*EI/L3    -3*EI/L2    0       3*EI/L3     0
                    0       0           0           0       0           0 ];
    case 'RE'
        Keloc = [   AE/L    0           0   -AE/L   0           0     
                    0       3*EI/L3     0   0       -3*EI/L3    3*EI/L2 
                    0       0           0   0       0           0 
                    -AE/L   0           0   AE/L    0           0      
                    0       -3*EI/L3    0   0       3*EI/L3     -3*EI/L2
                    0       3*EI/L2     0   0       -3*EI/L2    3*EI/L  ];
    otherwise
        error('Tipo de elemento no soportado.')
end
end