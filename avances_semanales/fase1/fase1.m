%{
================ PROYECTO FINAL DE ANÁLISIS DE VIBRACIONES ================
--------------------------------| FASE 1 |---------------------------------

Implementación del método de rigideces para análisis de estructuras
reticulares.

    Nelson Esteban Hernandez Soto
    Santiago Restrepo Villa
    Daniel Rojas Chica
    Maria Fernanda Villegas Loaiza
===========================================================================
%}
clear
clc
close all
%% === UNIDADES ===
% kN y m
%% == CONSTANTES ==
NL1 = 1; NL2 = 2;
X = 1; Y = 2; TH = 3;
g = 9.80665; % m/s²

nombre_archivo = 'entrada.xlsx';
xy_nod      = readtable(nombre_archivo, 'Sheet','xy_nod');
elementos   = readtable(nombre_archivo, 'Sheet','elementos');
prop_mat    = readtable(nombre_archivo, 'Sheet','prop_mat');
prop_sec    = readtable(nombre_archivo, 'Sheet','prop_sec');
carga_punt  = readtable(nombre_archivo, 'Sheet','carga_punt');
xy          = xy_nod{:, ["x" "y"]};
nno         = size(xy,1); % número de nodos

% Grados de libertad
ngdl    = 3*nno;
% fila = nodo
% col1 = gdl en direccion x
% col2 = gdl en direccion y
% col3 = gdl en direccion angular antihoraria
gdl     = reshape(1:ngdl, 3, nno)';

LaG     = elementos{:,["NL1" "NL2"]};
nelem   = size(LaG,1);
mat     = elementos{:,"material"};
sec     = elementos{:,"seccion"};
tipo    = elementos{:,"tipo"};

nmat    = size(prop_mat,1);
nsec    = size(prop_sec,1);
E   = prop_mat{:,"E"};
rho = prop_mat{:,"rho"};
A   = prop_sec{:,"A"};
I   = prop_sec{:,"I"};

f = zeros(ngdl,1);
ncp = size(carga_punt,1); % número de cargas puntuales
for cp = 1:ncp
    f(gdl(carga_punt{cp,"nodo"},carga_punt{cp,"direccion"})) = carga_punt{cp,"fuerza"};
end

%% Se dibuja la estructura junto con su numeracion
figure(1); 
hold on;
for e = 1:nelem
   line(xy(LaG(e,:),X), xy(LaG(e,:),Y));
   
   % Calculo la posicion del centro de gravedad de la barra
   cgx = (xy(LaG(e,NL1),X) + xy(LaG(e,NL2),X))/2;
   cgy = (xy(LaG(e,NL1),Y) + xy(LaG(e,NL2),Y))/2;   
   h = text(cgx, cgy, num2str(e)); set(h, 'Color', [0 0 1]);
end

axis equal
grid minor
plot(xy(:,X), xy(:,Y), 'ro');
text(xy(:,X), xy(:,Y), num2str((1:nno)'));
title('Numeración de la estructura');
