function [u, p, E] = balance_energia(m, c, k, Fs, Phi, n, acelerograma, r, nombre_archivo)
% En esta función se realiza el análisis dinámico y el balance de energía
% de la estructura ingresada
% -------------------------------------------------------------------------
% m: matriz de masa
% c: matriz de amortiguamiento
% k: matriz de rigidez
% Fs: vector de fuerzas estáticas aplicadas a la estructura
% Phi: matriz modal. Las filas son los modos, ordenadas de menor frecuencia
% a mayor frecuencia.
% n: número de modos a considerar
% c: vector de aceleraciones para el análisis dinámico
% r: vector que indica la dirección en que se considera el sismo. Las
% posiciones de los grados de libertad correspondientes llevan 1, las otras
% 0.
% nombre_archivo: para el título de la gráfica.

% CONSTANTES MÉTODO DE NEWMARK
% Ver Capítulo 16, Tabla 16.2.2 Dinámica de Anil K. Chopra
gammaN = 1/2;
%betaN = 1/6; % aceleración lineal
betaN = 1/4; % promedio de aceleración constante
ngdl = size(k,1); % número de grados de libertad
Phi = Phi(:,1:n); % matriz modal para el número de modos a considerar
M = Phi'*m*Phi; % matriz de masa por modos
C = Phi'*c*Phi; % matriz de amortiguamiento por modos
K = Phi'*k*Phi; % matriz de rigidez por modos
t = acelerograma{:,1}'; % s
a = acelerograma{:,2}'; % m/s²
dt = t(2)-t(1); % delta de tiempo
np = length(t); % número de puntos
% se inicializa el vector de fuerzas aplicadas por efecto de las aceleraciones
p = zeros(ngdl, np); 
for i=1:np % se rellena dicho vector según los grados de libertad en que se aplica la fuerza sísmica
    p(:,i)=-m*r*a(i)+Fs;
end

u = zeros(ngdl,np);  % desplazamiento inicial en cada gdl
du_t = zeros(ngdl,np);  % velocidad inicial en cada gdl
d2u_t = zeros(ngdl,np);  % aceleración inicial en cada gdl
% desplazamientos, velocidades y aceleraciones en coordenadas generalizadas
% por los modos
q = zeros(n,np); dq_t = zeros(n,np); d2q_t = zeros(n,np);
for j=1:n
    phi = Phi(:,j);
    q(j,1)=(phi'*m*u(:,1))/(phi'*m*phi);
    dq_t(j,1)=(phi'*m*du_t(:,1))/(phi'*m*phi);
end
P0 = Phi'*p(:,1);
% Se inicializan los vectores que contendrán los distintos tipos de energía
D = zeros(1,np); % kN.m energía disipada por amortiguamiento viscoso
Tr = zeros(1,np); % kN.m energía cinética relativa del sistema
R = zeros(1,np); % kN.m energía absorbida como resultado de las deformaciones de
                 % elementos estructurales.
Ir = zeros(1,np); % kN.m energía sísmica relativa de entrada
S = zeros(1,np); % kN.m trabajo hecho por las cargas presísmicas
% Se soluciona M*d2q_t0=P0-C*dq_t0-K*q0 para d2q_t
d2q_t(:,1) = M\(P0-C*dq_t(:,1)-K*q(:,1));
% variables de apoyo para el método de Newmark
a1 = 1/(betaN*dt^2)*M + gammaN/(betaN*dt)*C;
a2 = 1/(betaN*dt)*M + (gammaN/betaN-1)*C;
a3 = (1/(2*betaN)-1)*M+dt*(gammaN/(2*betaN)-1)*C;
K_hat = K+a1;
% Se programa el método de Newmark y la ecuación explícita de balance de
% energía relativa
% Tr(t)+D(t)+R(t)=Ir(t)+S(t)
for i=1:(np-1)
    % ANÁLISIS DINÁMICO (desplazamientos, velocidades y aceleraciones)
    P_hat = Phi'*p(:,i+1)+a1*q(:,i)+a2*dq_t(:,i)+a3*d2q_t(:,i);
    q(:,i+1) = K_hat\P_hat;
    dq_t(:,i+1) = gammaN/(betaN*dt)*(q(:,i+1)-q(:,i))+(1-gammaN/betaN)*dq_t(:,i)+dt*(1-gammaN/(2*betaN))*d2q_t(:,i);
    d2q_t(:,i+1) = 1/(betaN*dt^2)*(q(:,i+1)-q(:,i))-(1/(betaN*dt))*dq_t(:,i)-(1/(2*betaN)-1)*d2q_t(:,i);
    u(:,i+1) = Phi*q(:,1+i);
    du_t(:,i+1) = Phi*dq_t(:,1+i);
    d2u_t(:,i+1) = Phi*d2q_t(:,1+i);
    % BALANCE ENERGÍA
    Tr(i+1) = 1/2*du_t(:,i+1)'*m*du_t(:,i+1);
    D(i+1) = D(i)+1/2*(du_t(:,i)+du_t(:,i+1))'*c*(u(:,i+1)-u(:,i));
    % R = U + H. Se programa la fórmula de H.
    R(i+1) = R(i)+1/2*(u(:,i+1)-u(:,i))'*(k*u(:,i)+k*u(:,i+1));
    Ir(i+1) = Ir(i)-1/2*(u(:,i+1)-u(:,i))'*m*r*(a(i)+a(i+1));
    S(i+1) = S(i)+(u(:,i+1)-u(:,i))'*Fs;
end
% Vector de energías
E = [Tr; D; R; Ir; S];
% Error de balance de energía relativo
EBEr = abs(Ir+S-Tr-D-R)./abs(Ir)*100; % porcentaje
% GRÁFICOS
figure
subplot(2,3,[1,2,3])
plot(t, E')
title(sprintf('Balance de Energía %s', nombre_archivo))
xlabel('Tiempo (s)')
ylabel('kN m')
legend({'Tr(t): Energía cinética relativa del sistema', ...
    'D(t): Energía disipada por el amortiguamiento viscoso', ...
    'R(t): Energía absorbida por la deformación de los elementos estructurales', ...
    'Ir(t): Energía sísmica de entrada', ...
    'S(t): Trabajo realizado por cargas sísmicas preaplicadas'}, ...
    'Location', 'best')
subplot(2,3,4)
plot(t, [Tr' D' R' (Tr+D+R)'])
title('Interna')
xlabel('Tiempo (s)')
ylabel('kN m')
legend({'Tr(t)', 'D(t)', 'R(t)', 'Tr(t)+D(t)+R(t)'}, 'Location', 'best')
subplot(2,3,5)
plot(t, [Ir' S' (Ir+S)'])
title('Externa')
xlabel('Tiempo (s)')
ylabel('kNm')
legend({'Ir(t)', 'S(t)','Ir(t)+S(t)'}, 'Location', 'best')
subplot(2,3,6)
% Dado que inicialmente el balance aún se está calibrando, se grafica el
% error a partir de la tercera iteración.
plot(t(3:end), EBEr(3:end))
title('Error de balance de energía')
xlabel('Tiempo (s)')
ylabel('%')
legend({'EBEr %'}, 'Location', 'best')
end
