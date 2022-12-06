function [u, p, E] = balance_energia(m, c, k, Fs, Phi, n, acelerograma, r)
% CONSTANTES MÉTODO DE NEWMARK
% Ver Capítulo 16, Tabla 16.2.2 Dinámica de Anil K. Chopra
% iota es uno si en dicha fila se considera la aceleración, en otro caso, 0
gammaN = 1/2;
%betaN = 1/6; % aceleración lineal
betaN = 1/4; % promedio de aceleración constante
ngdl = size(k,1);
Phi = Phi(:,1:n);
M = Phi'*m*Phi;
C = Phi'*c*Phi;
K = Phi'*k*Phi;
t = acelerograma{:,1}';
a = acelerograma{:,2}';
dt = t(2)-t(1);
np = length(t);
p = zeros(ngdl, np);
% for i=1:3:ngdl
%     a_t(i,:)=a;
% end
for i=1:np
    p(:,i)=-m*r*a(i)+Fs;
end

u = zeros(ngdl,np);  % desplazamiento inicial en cada gdl
du_t = zeros(ngdl,np);  % desplazamiento inicial en cada gdl
d2u_t = zeros(ngdl,np);  % desplazamiento inicial en cada gdl
du0_t = zeros(ngdl,1); % velocidad inicial en cada gdl
q = zeros(n,np); dq_t = zeros(n,np); d2q_t = zeros(n,np);
for j=1:n
    phi = Phi(:,j);
    q(j,1)=(phi'*m*u(:,1))/(phi'*m*phi);
    dq_t(j,1)=(phi'*m*du0_t)/(phi'*m*phi);
end
P0 = Phi'*p(:,1);
D = zeros(1,np); % energía disipada por amortiguamiento viscoso
Tr = zeros(1,np); % energía cinética relativa del sistema
R = zeros(1,np); % energía absorbida como resultado de las deformaciones de
                 % elementos estructurales.
Ir = zeros(1,np); % energía sísmica relativa de entrada
S = zeros(1,np); % trabajo hecho por las cargas presísmicas
% Se soluciona M*d2q_t0=P0-C*dq_t0-K*q0 para d2q_t
d2q_t(:,1) = M\(P0-C*dq_t(:,1)-K*q(:,1));

a1 = 1/(betaN*dt^2)*M + gammaN/(betaN*dt)*C;
a2 = 1/(betaN*dt)*M + (gammaN/betaN-1)*C;
a3 = (1/(2*betaN)-1)*M+dt*(gammaN/(2*betaN)-1)*C;
K_hat = K+a1;
% Se programa el método de Newmark y la ecuación explícita de balance de
% energía 
% Tr(t)+D(t)+R(t)=Ir(t)+S(t)
for i=1:(np-1)
    P_hat = Phi'*p(:,i+1)+a1*q(:,i)+a2*dq_t(:,i)+a3*d2q_t(:,i);
    q(:,i+1) = K_hat\P_hat;
    dq_t(:,i+1) = gammaN/(betaN*dt)*(q(:,i+1)-q(:,i))+(1-gammaN/betaN)*dq_t(:,i)+dt*(1-gammaN/(2*betaN))*d2q_t(:,i);
    d2q_t(:,i+1) = 1/(betaN*dt^2)*(q(:,i+1)-q(:,i))-(1/(betaN*dt))*dq_t(:,i)-(1/(2*betaN)-1)*d2q_t(:,i);
    u(:,i+1) = Phi*q(:,1+i);
    du_t(:,i+1) = Phi*dq_t(:,1+i);
    d2u_t(:,i+1) = Phi*d2q_t(:,1+i);
    Tr(i+1) = 1/2*du_t(:,i+1)'*m*du_t(:,i+1);
    D(i+1) = D(i)+1/2*(du_t(:,i)+du_t(:,i+1))'*c*(u(:,i+1)-u(:,i));
    R(i+1) = R(i)+(u(:,i+1)-u(:,i))'*k*u(:,i+1);
    Ir(i+1) = Ir(i)-1/2*(u(:,i+1)-u(:,i))'*m*r*(a(i)+a(i+1));
    S(i+1) = S(i)+(u(:,i+1)-u(:,i))'*Fs;
end
E = [Tr; D; R; Ir; S];
figure
subplot(2,3,[1,2,3])
plot(t, E')
title('Balance de Energía')
xlabel('Tiempo (s)')
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
legend({'Tr(t)', 'D(t)', 'R(t)', 'Tr(t)+D(t)+R(t)'}, 'Location', 'best')
subplot(2,3,5)
plot(t, [Ir' S' (Ir+S)'])
title('Externa')
xlabel('Tiempo (s)')
legend({'Ir(t)', 'S(t)','Ir(t)+S(t)'}, 'Location', 'best')
subplot(2,3,6)
plot(t, (Tr+D+R)-(Ir+S))
title('Diferencia')
xlabel('Tiempo (s)')
legend({'Tr(t)+D(t)+R(t)-Ir(t)-S(t)'}, 'Location', 'best')
end
