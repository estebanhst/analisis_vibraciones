%{
================ PROYECTO FINAL DE ANÁLISIS DE VIBRACIONES ================
--------------------------------| FASE 3 |---------------------------------

Implementación de un método de análisis de balance de energía.

    Nelson Esteban Hernandez Soto
    Santiago Restrepo Villa
    Daniel Rojas Chica
    Maria Fernanda Villegas Loaiza
===========================================================================
%}
clear
clc
%% === UNIDADES ===
% kN y m
%% == CONSTANTES ==
NL1 = 1; NL2 = 2;
X = 1; Y = 2; TH = 3;
g = 9.80665; % m/s²
% matriz_masa = 'condensada';
matriz_masa = 'consistente';

%nombre_archivo = 'entrada.xlsx';
nombre_archivo = 'entrada.xlsx';
xy_nod      = readtable(nombre_archivo, 'Sheet','xy_nod');
elementos   = readtable(nombre_archivo, 'Sheet','elementos');
prop_mat    = readtable(nombre_archivo, 'Sheet','prop_mat');
prop_sec    = readtable(nombre_archivo, 'Sheet','prop_sec');
carga_punt  = readtable(nombre_archivo, 'Sheet','carga_punt');
carga_distr = readtable(nombre_archivo, 'Sheet','carga_distr');
restricciones = readtable(nombre_archivo, 'Sheet','restricciones');
acelerograma=readtable(nombre_archivo, 'Sheet','acelerograma');
p=acelerograma{:,"aceleracion"};         %Vector de registros del acelerograma del sismo
xy          = xy_nod{:, ["x" "y"]};
nno         = int16(size(xy,1)); % número de nodos

% Grados de libertad
ngdl    = 3*nno;
% fila = nodo
% col1 = gdl en direccion x
% col2 = gdl en direccion y
% col3 = gdl en direccion angular antihoraria
gdl     = reshape(1:ngdl, 3, nno)';

LaG     = elementos{:,["NL1" "NL2"]};
nelem   = int16(size(LaG,1));
mat     = elementos{:,"material"};
sec     = elementos{:,"seccion"};
tipo    = elementos{:,"tipo"};
% Para la masa de entrepiso
Af      = elementos{:,"Af"};
espesor = elementos{:,"espesor"};
mat_losa= elementos{:,"material_losa"};

nmat    = int16(size(prop_mat,1));
nsec    = int16(size(prop_sec,1));
E   = prop_mat{:,"E"};
rho = prop_mat{:,"rho"};
A   = prop_sec{:,"A"};
I   = prop_sec{:,"I"};

%% Se dibuja la estructura junto con su numeracion
dibujar_numeracion(xy, LaG);

% vector de fuerzas nodales equivalentes global
f = zeros(ngdl,1);
ncp = int16(size(carga_punt,1)); % número de cargas puntuales
for cp = 1:ncp
    f(gdl(carga_punt{cp,"nodo"},carga_punt{cp,"direccion"})) = carga_punt{cp,"fuerza"};
end

% fuerzas distribuidas aplicadas sobre los elementos en coordenadas locales
b1 = carga_distr{:,"b1"};
b2 = carga_distr{:,"b2"};
q1 = carga_distr{:,"q1"};
q2 = carga_distr{:,"q2"};

%% Separo la memoria
K   = zeros(ngdl);      % matriz de rigidez global
M   = zeros(ngdl);      % matriz de masa consistente global
Ke  = cell(nelem,1);    % matriz de rigidez local en coordenadas globales
Me  = cell(nelem,1);    % matriz de masa consistente en coordenadas globales
T   = cell(nelem,1);    % matriz de transformacion de coordenadas
idx = cell(nelem,1);    % almacena los 6 gdls de las barras
fe  = cell(nelem,1);    % fuerzas nodales equivalentes globales de cada elemento 
x1 = xy(LaG(:,NL1),X); y1 = xy(LaG(:,NL1),Y);
x2 = xy(LaG(:,NL2),X); y2 = xy(LaG(:,NL2),Y);
L = hypot(x2-x1,y2-y1);
%% ensamblo la matriz de rigidez global (K) y vector de fuerzas global (f)
for e = 1:nelem  % para cada elemento
   % saco los 6 gdls del elemento e
   idx{e} = [gdl(LaG(e,NL1),:) gdl(LaG(e,NL2),:)];
   
   % matriz de transformacion de coordenadas para la barra e
   T{e} = calc_Te(x1(e), y1(e), x2(e), y2(e));
         
   % matriz de rigidez local expresada en el sistema de coordenadas locales
   % para el elemento e en kN/m
   Kloc = calc_Keloc(tipo{e}, L(e), A(sec(e)), E(mat(e)), I(sec(e)));

   % matriz de masa consistente expresada en el sistema de coordenadas locales
   % para el elemento e en kg=N*s²/m -> se convierte a kN*s²/m
   Mloc = calc_Meloc(tipo{e}, L(e), A(sec(e)), I(sec(e)), rho(mat(e)), matriz_masa)/1000;
   % matriz de masa adicional para las vigas o elementos que reciban el
   % peso del entrepiso
   if isnan(Af(e))
       M_losa = zeros(size(Mloc));
   else
       % kN*s²/m
       M_losa = masa_losa(rho(mat_losa(e)),espesor(e),Af(e),L(e))/1000;
   end
   % Inclusión de las fuerzas por peso propio
   wx = rho(mat(e))*A(sec(e))*g*(y2(e)-y1(e))/L(e)/1000; % kN/m
   b1(e) = b1(e)-wx;
   b2(e) = b2(e)-wx;
   wy = rho(mat(e))*A(sec(e))*g*(x2(e)-x1(e))/L(e)/1000; % kN/m
   q1(e) = q1(e)-wy;
   q2(e) = q2(e)-wy;

   % vector de fuerzas nodales equivalentes en coordenadas locales
   feloc = calc_feloc(tipo{e}, L(e), b1(e), b2(e), q1(e), q2(e));

   % Cambio a coordenadas globales
   Ke{e} = T{e}'*Kloc*T{e};
   Me{e} = T{e}'*(Mloc+M_losa)*T{e};
   fe{e} = T{e}'*feloc;
   
   K(idx{e},idx{e}) = K(idx{e},idx{e}) + Ke{e}; % sumo Ke{e} a K global
   M(idx{e},idx{e}) = M(idx{e},idx{e}) + Me{e}; % sumo Me{e} a M global
   f(idx{e})        = f(idx{e})        + fe{e}; % sumo a f global
end

%% Restricciones y grados de libertad con desplazamientos conocidos
nres = int16(size(restricciones, 1)); % número de restricciones
c = zeros(nres,1, 'int16');
for i=1:nres
    c(i) = gdl(restricciones{i,"nodo"}, restricciones{i,"direccion"});
end

% restricciones conocidas
ac = restricciones{:,"valor"};

% grados de libertad desconocidos
d = setdiff((1:ngdl)',c);

%% Extraigo las submatrices y especifico las cantidades conocidas
% f = vector de fuerzas nodales equivalentes
% q = vector de fuerzas nodales de equilibrio del elemento
% a = desplazamientos

%| qd |   | Kcc Kcd || ac |   | fd |    recuerde que siempre qc=0
%|    | = |         ||    | - |    |
%| qc |   | Kdc Kdd || ad |   | fc |    en este caso en particular fd=0

Kcc = K(c,c); Kcd = K(c,d); fd = f(c);
Kdc = K(d,c); Kdd = K(d,d); fc = f(d);

%% resuelvo el sistema de ecuaciones
ad = Kdd\(fc - Kdc*ac);    % calculo desplazamientos desconocidos
qd = Kcc*ac + Kcd*ad - fd; % calculo fuerzas de equilibrio desconocidas

% armo los vectores de desplazamientos (a) y fuerzas (q)
a = zeros(ngdl,1);  q = zeros(ngdl,1);   % separo la memoria
a(c) = ac;   a(d) = ad; % desplazamientos 
q(c) = qd;  %q(d) = qc; % fuerzas nodales de equilibrio

%% Fuerzas internas en cada elemento
qe_loc  = cell(nelem,1);
qe_glob = cell(nelem,1);
for e = 1:nelem % para cada barra
   %fprintf('\n\n Fuerzas internas para elemento %d en coord. globales = \n', e);
   qe_glob{e} = Ke{e}*a(idx{e,:}) - fe{e};
   %disp(qe_glob{e})
   
   %fprintf('\n\n Fuerzas internas para elemento %d en coord. locales = \n', e);
   qe_loc{e} = T{e}*qe_glob{e};
   %disp(qe_loc{e});   
end

%% imprimo los resultados
format short
disp('Desplazamientos nodales                                            ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
vect_mov = reshape(a,3,nno)'; % vector de movimientos
for i = 1:nno
   fprintf('Nodo %3d: u = %12.4g mm, v = %12.4g mm, theta = %12.4g rad \n', ...
      i, 1000*vect_mov(i,X), 1000*vect_mov(i,Y), vect_mov(i,TH));
end

disp(' ');
disp('Fuerzas nodales de equilibrio (solo imprimo los diferentes de cero)')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
qq = reshape(q,3,nno)';
for i = 1:nno   
   if ~all(abs(qq(i,:) - [0 0 0]) < 1e-5)
      fprintf('Nodo %3d qx = %12.4g kN, qy = %12.4g kN, mom = %12.4g kN*m\n', ...
         i, qq(i,X), qq(i,Y), qq(i,TH));
   end
end

%% Dibujar la estructura y su deformada
esc_def    = 10;               % escalamiento de la deformada
esc_faxial = 0.00005;           % escalamiento del diagrama de axiales
esc_V      = 0.001;           % escalamiento del diagrama de cortantes
esc_M      = 0.001;           % escalamiento del diagrama de momentos

%xdef = xnod + esc_def*vect_mov(:,[X Y]);

figure(2); hold on; title(sprintf('Deformada exagerada %d veces', esc_def));    xlabel('x, m'); ylabel('y, m'); axis equal
%figure(3); hold on; title('Fuerza axial [kN]');      xlabel('x, m'); ylabel('y, m'); axis equal
%figure(4); hold on; title('Fuerza cortante [kN]');   xlabel('x, m'); ylabel('y, m'); axis equal
%figure(5); hold on; title('Momento flector [kN-m]'); xlabel('x, m'); ylabel('y, m'); axis equal

for e = 1:nelem
   dibujar_graficos(tipo{e},...
      A(sec(e)), E(mat(e)), I(sec(e)), ...
      x1(e),y1(e), x2(e),y2(e), b1(e), b2(e), q1(e), q2(e),...
      qe_loc{e}, T{e}*a(idx{e,:}), ...
      esc_def, esc_faxial, esc_V, esc_M);
end

%% ANÁLISIS MODAL
n = 12; % número de modos a tener en cuenta
Kdd = K(d,d); Mdd = M(d,d);
[Phi, lams] = eig(Kdd, Mdd);
omega = sqrt(diag(lams)); % frecuencias angulares rad/s
[omega,iModo] = sort(omega);% ordena las frecuencias
T_mod = 2*pi./omega;            % s - periodo de vibración
f_mod = 1./T_mod;           % frecuencias Hz
Phi = Phi(:,iModo);         % ordena los modos según el orden de las frecuencias
% participación modal
alfa = Phi'*Mdd*ones(size(Phi(:,1)));
M_mod_efectiva = alfa.^2;
participacion_masa = M_mod_efectiva/sum(M_mod_efectiva);

%% Matriz de Amortiguamiento (Método de Rayleigh, Libro de Daryl L. Logan):
zeta=0.05;  %Amortiguamiento
w1=omega(1,1);  %Frecuencia del primer modo de vibrar
w2=omega(2,1);  %Frecuencia del segundo modo de vibrar

%Después de resolver el sistema de ecuaciones para hallar alfa0 y alfa1
a0=((2*w1*w2)/(w2^2-w1^2))*(w2*zeta-w1*zeta);
a1=(2/(w2^2-w1^2))*(w2*zeta-w1*zeta);

%Ensamblaje Matriz de Amortiguamiento:
% a0*M corresponde al amortiguamiento externo y a1K es el
% amortiguamiento material o interno.
C=a0*M+a1*K; % kN*s/m
Cdd=C(d,d);

%% Análisis dinámico. El acelerograma se ingresa en el Excel en m/s² o se define como función
% n: número de modos a considerar
% COMENTAR PARA USAR EL SISMO INGRESADO EN EL EXCEL
%{ USANDO FUNCIÓN SINUSOIDAL
% FUNCIÓN SINUSOIDAL
t = (0:0.01:15)'; % s
v = 2*pi; % frecuencia angular de la función sinusoidal
aceleracion = 0.5*sin(v*t);
% los últimos 5 segundos se hacen 0 las aceleraciones
%aceleracion(t>10) = 0;
%}
acelerograma = table(t, aceleracion, 'VariableNames',{'tiempo', 'aceleracion'});
u = zeros(ngdl, size(acelerograma,1)); p = zeros(size(u));
r = zeros(ngdl,1); r(gdl(:,X))=1;
[u(d,:), p(d,:), B_E] = balance_energia(Mdd, Cdd, Kdd, fc, Phi, n, acelerograma, r(d));
% u: desplazamientos en metros
% p: fuerza en kN
% cada fila corresponde a un grado de libertad. Las aceleraciones del sismo
% fueron colocadas únicamente en los grados de libertad horizontales
% (dirección 1)
figure
title('Desplazamientos en el último piso')
xlabel('Tiempo (s)')
ylabel('Desplazamiento (mm)')
nodos_ultimo_piso = find(xy(:,Y)==max(xy(:,Y)))';
hold on
legend
for i = nodos_ultimo_piso   
    plot(acelerograma{:,'tiempo'}', u(gdl(i,X),:)*1000, 'DisplayName', sprintf('Nodo %d',i))
end

figure
title('Desplazamientos en los nodos, dirección X')
xlabel('Tiempo (s)')
ylabel('Desplazamiento (mm)')
nodos_ultimo_piso = find(xy(:,Y)==max(xy(:,Y)))';
hold on
legend
for i = 1:nno   
    plot(acelerograma{:,'tiempo'}', u(gdl(i,X),:)*1000, 'DisplayName', sprintf('Nodo %d',i))
end