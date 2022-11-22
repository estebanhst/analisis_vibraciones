function dibujar_numeracion(xy, LaG)
X = 1; NL1 = 1; Y=2; NL2=2; nno = size(xy,1); nelem = size(LaG,1);
figure; 
hold on;
for e = 1:nelem
   line(xy(LaG(e,:),X), xy(LaG(e,:),Y), 'LineWidth',2);
   
   % Calculo la posicion del centro de gravedad de la barra
   cgx = (xy(LaG(e,NL1),X) + xy(LaG(e,NL2),X))/2;
   cgy = (xy(LaG(e,NL1),Y) + xy(LaG(e,NL2),Y))/2;   
   h = text(cgx, cgy, num2str(e)); set(h, 'Color', [0 0 1]);
end

axis equal
grid minor
plot(xy(:,X), xy(:,Y), 'ro');
text(xy(:,X), xy(:,Y), num2str((1:nno)'));
title('NUMERACIÓN DE LA ESTRUCTURA');
end