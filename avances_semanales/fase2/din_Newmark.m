function [u, p] = din_Newmark(m, c, k, Phi, n, acelerograma, comp)
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
    p(:,i)=-m*comp*a(i);
end


u = zeros(ngdl,np);  % desplazamiento inicial en cada gdl
du0_t = zeros(ngdl,1); % velocidad inicial en cada gdl
q = zeros(n,np); dq_t = zeros(n,np); d2q_t = zeros(n,np);
for j=1:n
    phi = Phi(:,j);
    q(j,1)=(phi'*m*u(:,1))/(phi'*m*phi);
    dq_t(j,1)=(phi'*m*du0_t)/(phi'*m*phi);
end
P0 = Phi'*p(:,1);
% Se soluciona M*d2q_t0=P0-C*dq_t0-K*q0 para d2q_t
d2q_t(:,1) = M\(P0-C*dq_t(:,1)-K*q(:,1));

a1 = 1/(betaN*dt^2)*M + gammaN/(betaN*dt)*C;
a2 = 1/(betaN*dt)*M + (gammaN/betaN-1)*C;
a3 = (1/(2*betaN)-1)*M+dt*(gammaN/(2*betaN)-1)*C;
K_hat = K+a1;
for i=1:(np-1)
    P_hat = Phi'*p(:,i+1)+a1*q(:,i)+a2*dq_t(:,i)+a3*d2q_t(:,i);
    q(:,i+1) = K_hat\P_hat;
    dq_t(:,i+1) = gammaN/(betaN*dt)*(q(:,i+1)-q(:,i))+(1-gammaN/betaN)*dq_t(:,i)+dt*(1-gammaN/(2*betaN))*d2q_t(:,i);
    d2q_t(:,i+1) = 1/(betaN*dt^2)*(q(:,i+1)-q(:,i))-(1/(betaN*dt))*dq_t(:,i)-(1/(2*betaN)-1)*d2q_t(:,i);
    u(:,i+1) = Phi*q(:,1+i);
end
end