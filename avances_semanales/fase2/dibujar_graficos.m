function dibujar_graficos(tipo, A, E, I, x1,y1, x2,y2, b1,b2, q1,q2, qe, ae, esc_def, esc_faxial, esc_V, esc_M)
%{
Esta función dibuja el elemento de pórtico/cercha deformado junto con sus 
respectivos diagramas de fuerza axial, fuerza cortante y momento flector.
El diagrama de momento flector se grafica en el lado opuesto de la fibra
a tracción
PARAMETROS DE ENTRADA (junto con algunos ejemplos):
tipo    = 'EE' (pórtico)
        = 'RR' (cercha)
        = 'RE' (rotula/empotrado)
        = 'ER' (empotrado/rotula)
A = area
E = E
I = Ix local
(x1,y1) y (x2,y2) son las coordenadas de los puntos iniciales
b1, b2    carga axial    en x1 y x2
q1, q2    carga vertical en x1 y x2  
qe = [ 0.01,        % U1, V1, M1 reacciones del nodo 1 en coord. locales
      -0.01,
       0.04,
      -0.01,        % U2, V2, M2 reacciones del nodo 2 en coord. locales
       0.02,
      -0.07 ]
ae = [ 0.01,        % u1, v1, t1 desplazamientos nodo 1 en coord. locales
      -0.01,
       0.04,
      -0.01,        % u2, v2, t2 desplazamientos nodo 2 en coord. locales
       0.02,
      -0.07 ]
esc_def    = 10     % escalamiento de la deformada
esc_faxial = 10     % escalamiento del diagrama de axiales
esc_V      = 10     % escalamiento del diagrama de cortantes
esc_M      = 10     % escalamiento del diagrama de momentos
con_lineas = {True/False} % graficar lineas verticales en diagramas de M/V/A
%}
X = 1; Y = 2;
X1 = 1; Y1 = 2; M1 = 3; X2 = 4; Y2 = 5; M2   = 6; 
npuntos = 1001;
L = hypot(x2-x1, y2-y1);
x = linspace(0, L, npuntos);
s = x;

ae = num2cell(ae);
% se asignan los desplazamientos en coordenadas locales
[u1, v1, t1, u2, v2, t2] = ae{:};

% se calculan mediante fórmulas las fuerzas, el momento y los desplazamientos
%q = q1 - (x*(q1 - q2))/L
switch tipo
    case 'EE'
        axial = ((b1 - b2)*x.^2)/(2*L) - b1*x + (2*L^2*b1 + L^2*b2 - 6*A*E*u1 + 6*A*E*u2)/(6*L);
        V     = q1*x - (7*L^4*q1 + 3*L^4*q2 - 240*E*I*v1 + 240*E*I*v2 - 120*E*I*L*t1 - 120*E*I*L*t2)/(20*L^3) - (x.^2*(q1 - q2))/(2*L);
        M     = (3*L^5*q1 + 2*L^5*q2 - 10*L^2*q1*x.^3 + 30*L^3*q1*x.^2 + 10*L^2*q2*x.^3 - 21*L^4*q1*x - 9*L^4*q2*x - 240*E*I*L^2*t1 - 120*E*I*L^2*t2 - 360*E*I*L*v1 + 360*E*I*L*v2 + 720*E*I*v1*x - 720*E*I*v2*x + 360*E*I*L*t1*x + 360*E*I*L*t2*x)/(60*L^3);
        u     = (b1*x.^3 - b2*x.^3 - 3*L*b1*x.^2 + 2*L^2*b1*x + L^2*b2*x + 6*A*E*L*u1 - 6*A*E*u1*x + 6*A*E*u2*x)/(6*A*E*L);
        v     = (5*L^3*q1*x.^4 - L^2*q1*x.^5 - 7*L^4*q1*x.^3 + 3*L^5*q1*x.^2 + L^2*q2*x.^5 - 3*L^4*q2*x.^3 + 2*L^5*q2*x.^2 + 120*E*I*L^3*v1 + 240*E*I*v1*x.^3 - 240*E*I*v2*x.^3 + 120*E*I*L*t1*x.^3 + 120*E*I*L^3*t1*x + 120*E*I*L*t2*x.^3 - 360*E*I*L*v1*x.^2 + 360*E*I*L*v2*x.^2 - 240*E*I*L^2*t1*x.^2 - 120*E*I*L^2*t2*x.^2)/(120*E*I*L^3);
        %t    = (20*L^3*q1*x.^3 - 5*L^2*q1*x.^4 - 21*L^4*q1*x.^2 + 5*L^2*q2*x.^4 - 9*L^4*q2*x.^2 + 6*L^5*q1*x + 4*L^5*q2*x + 120*E*I*L^3*t1 + 720*E*I*v1*x.^2 - 720*E*I*v2*x.^2 - 720*E*I*L*v1*x + 720*E*I*L*v2*x + 360*E*I*L*t1*x.^2 - 480*E*I*L^2*t1*x + 360*E*I*L*t2*x.^2 - 240*E*I*L^2*t2*x)/(120*E*I*L^3);
    case 'RR'
        axial = ((b1 - b2)*x.^2)/(2*L) - b1*x + (2*L^2*b1 + L^2*b2 - 6*A*E*u1 + 6*A*E*u2)/(6*L);
        V     = q1*x - (L*q2)/6 - (L*q1)/3 - (x.^2*(q1 - q2))/(2*L);
        M     = -(x*(L - x)*(2*L*q1 + L*q2 - q1*x + q2*x))/(6*L);
        u     = (b1*x.^3 - b2*x.^3 - 3*L*b1*x.^2 + 2*L^2*b1*x + L^2*b2*x + 6*A*E*L*u1 - 6*A*E*u1*x + 6*A*E*u2*x)/(6*A*E*L);
        v     = (3*q2*x.^5 - 3*q1*x.^5 - 20*L^2*q1*x.^3 - 10*L^2*q2*x.^3 + 15*L*q1*x.^4 + 8*L^4*q1*x + 7*L^4*q2*x + 360*E*I*L*v1 - 360*E*I*v1*x + 360*E*I*v2*x)/(360*E*I*L);
        %t    = (8*L^4*q1 + 7*L^4*q2 - 15*q1*x.^4 + 15*q2*x.^4 - 60*L^2*q1*x.^2 - 30*L^2*q2*x.^2 - 360*E*I*v1 + 360*E*I*v2 + 60*L*q1*x.^3)/(360*E*I*L);
    case 'RE'
        axial = ((b1 - b2)*x.^2)/(2*L) - b1*x + (2*L^2*b1 + L^2*b2 - 6*A*E*u1 + 6*A*E*u2)/(6*L);
        V     = q1*x - (11*L^4*q1 + 4*L^4*q2 - 120*E*I*v1 + 120*E*I*v2 - 120*E*I*L*t2)/(40*L^3) - (x.^2*(q1 - q2))/(2*L);
        M     = -(x*(33*L^4*q1 + 12*L^4*q2 + 20*L^2*q1*x.^2 - 20*L^2*q2*x.^2 - 360*E*I*v1 + 360*E*I*v2 - 60*L^3*q1*x - 360*E*I*L*t2))/(120*L^3);
        u     = (b1*x.^3 - b2*x.^3 - 3*L*b1*x.^2 + 2*L^2*b1*x + L^2*b2*x + 6*A*E*L*u1 - 6*A*E*u1*x + 6*A*E*u2*x)/(6*A*E*L);
        v     = (10*L^3*q1*x.^4 - 2*L^2*q1*x.^5 - 11*L^4*q1*x.^3 + 2*L^2*q2*x.^5 - 4*L^4*q2*x.^3 + 3*L^6*q1*x + 2*L^6*q2*x + 240*E*I*L^3*v1 + 120*E*I*v1*x.^3 - 120*E*I*v2*x.^3 + 120*E*I*L*t2*x.^3 - 120*E*I*L^3*t2*x - 360*E*I*L^2*v1*x + 360*E*I*L^2*v2*x)/(240*E*I*L^3);
        %t    = (3*L^6*q1 + 2*L^6*q2 - 10*L^2*q1*x.^4 + 40*L^3*q1*x.^3 - 33*L^4*q1*x.^2 + 10*L^2*q2*x.^4 - 12*L^4*q2*x.^2 - 120*E*I*L^3*t2 - 360*E*I*L^2*v1 + 360*E*I*L^2*v2 + 360*E*I*v1*x.^2 - 360*E*I*v2*x.^2 + 360*E*I*L*t2*x.^2)/(240*E*I*L^3);
    case 'ER'
        axial = ((b1 - b2)*x.^2)/(2*L) - b1*x + (2*L^2*b1 + L^2*b2 - 6*A*E*u1 + 6*A*E*u2)/(6*L);
        V     = q1*x - (16*L^4*q1 + 9*L^4*q2 - 120*E*I*v1 + 120*E*I*v2 - 120*E*I*L*t1)/(40*L^3) - (x.^2*(q1 - q2))/(2*L);
        M     = -((L - x)*(20*L^2*q2*x.^2 - 7*L^4*q2 - 20*L^2*q1*x.^2 - 8*L^4*q1 + 360*E*I*v1 - 360*E*I*v2 + 40*L^3*q1*x + 20*L^3*q2*x + 360*E*I*L*t1))/(120*L^3);
        u     = (b1*x.^3 - b2*x.^3 - 3*L*b1*x.^2 + 2*L^2*b1*x + L^2*b2*x + 6*A*E*L*u1 - 6*A*E*u1*x + 6*A*E*u2*x)/(6*A*E*L);
        v     = (10*L^3*q1*x.^4 - 2*L^2*q1*x.^5 - 16*L^4*q1*x.^3 + 8*L^5*q1*x.^2 + 2*L^2*q2*x.^5 - 9*L^4*q2*x.^3 + 7*L^5*q2*x.^2 + 240*E*I*L^3*v1 + 120*E*I*v1*x.^3 - 120*E*I*v2*x.^3 + 120*E*I*L*t1*x.^3 + 240*E*I*L^3*t1*x - 360*E*I*L*v1*x.^2 + 360*E*I*L*v2*x.^2 - 360*E*I*L^2*t1*x.^2)/(240*E*I*L^3);
        %t    = (40*L^3*q1*x.^3 - 10*L^2*q1*x.^4 - 48*L^4*q1*x.^2 + 10*L^2*q2*x.^4 - 27*L^4*q2*x.^2 + 16*L^5*q1*x + 14*L^5*q2*x + 240*E*I*L^3*t1 + 360*E*I*v1*x.^2 - 360*E*I*v2*x.^2 - 720*E*I*L*v1*x + 720*E*I*L*v2*x + 360*E*I*L*t1*x.^2 - 720*E*I*L^2*t1*x)/(240*E*I*L^3);
    otherwise
        error('Tipo de elemento no soportado')
end
%theta = atan(t)  % Angulo de giro [rad]

% rotación de la solución antes de dibujar
ang = atan2(y2-y1,x2-x1);
T       = [ cos(ang) -sin(ang)
            sin(ang) cos(ang)];
      
% Dibujar deformada
figure(2)
pos = T*[   s+esc_def*u
            esc_def*v];

xx = pos(X,:) + x1;
yy = pos(Y,:) + y1;

plot([x1 x2], [y1 y2], 'b-', xx, yy, 'r-','LineWidth',2);

% % Dibujar los diagramas de fuerza axial 
% figure(3)
% pos = T*[ s; esc_faxial*axial ]; % escalamiento del diagrama
% 
% ss = pos(X,:) + x1;
% aa = pos(Y,:) + y1;
% 
% plot([x1 x2], [y1 y2], 'b-', [x1 ss x2], [y1 aa y2], 'r-','LineWidth',2);
% text(ss(1),   aa(1),   num2str(-qe(X1)));
% text(ss(end), aa(end), num2str(+qe(X2)));

% % Dibujar los diagramas de fuerza cortante
% figure(4)
% pos = T*[ s; esc_V*V ]; % escalamiento del diagrama
% 
% ss = pos(X,:) + x1;
% vv = pos(Y,:) + y1;
% 
% plot([x1 x2], [y1 y2], 'b-', [x1 ss x2], [y1 vv y2], 'r-','LineWidth',2);
% text(ss(1),   vv(1),   num2str(+qe(Y1)));
% text(ss(end), vv(end), num2str(-qe(Y2)));

% %% Dibujar los diagramas de momento flector
% figure(5)
% pos = T*[ s; esc_M*M ]; % escalamiento del diagrama
% 
% ss = pos(X,:) + x1;
% mm = pos(Y,:) + y1;
% 
% plot([x1 x2], [y1 y2], 'b-', [x1 ss x2], [y1 mm y2], 'r-','LineWidth',2);
% text(ss(1),   mm(1),   num2str(-qe(M1)));
% text(ss(end), mm(end), num2str(+qe(M2)));
% [minM,idminM] = min(M); text(ss(idminM), mm(idminM), num2str(minM));
% [maxM,idmaxM] = max(M); text(ss(idmaxM), mm(idmaxM), num2str(maxM));
end