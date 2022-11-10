function Te = calc_Te(x1, y1, x2, y2)
%calc_Te Matriz de transformación para un elemento con 3 gld por nodo
%       PARAMETROS:
%    (x1,y1) y (x2,y2) son las coordenadas de los extremos de la barra
L = hypot(x2-x1,y2-y1);
c = (x2-x1)/L; % coseno de la inclinación
s = (y2-y1)/L; % seno de la inclinación
Te = [  c   s   0   0   0   0
       -s   c   0   0   0   0
        0   0   1   0   0   0
        0   0   0   c   s   0
        0   0   0  -s   c   0
        0   0   0   0   0   1];
end