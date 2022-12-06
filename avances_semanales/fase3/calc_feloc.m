function feloc = calc_feloc(tipo, L, b1, b2, q1, q2)
%calc_feloc Calcula el vector de fuerzas nodales equivalentes para cargas
%distribuidas en un elemento
switch tipo
    case 'EE'
        X1 =  (L*(2*b1 + b2))/6;
        Y1 =  (L*(7*q1 + 3*q2))/20;
        M1 =  (L^2 * (3*q1 + 2*q2))/60;
        X2 =  (L*(b1 + 2*b2))/6;
        Y2 =  (L*(3*q1 + 7*q2))/20;
        M2 = -(L^2 * (2*q1 + 3*q2))/60;
    case 'RR'
        X1 = (L*(2*b1 + b2))/6;
        Y1 = (L*(2*q1 + q2))/6;
        M1 = 0;
        X2 = (L*(b1 + 2*b2))/6; 
        Y2 = (L*(q1 + 2*q2))/6;
        M2 = 0;
    case 'RE'
        X1 = (L*(2*b1 + b2))/6;
        Y1 = (L*(11*q1 + 4*q2))/40;
        M1 = 0;
        X2 = (L*(b1 + 2*b2))/6;
        Y2 = (L*(9*q1 + 16*q2))/40;
        M2 = -(L^2*(7*q1 + 8*q2))/120;
    case 'ER'
        X1 = (L*(2*b1 + b2))/6;
        Y1 = (L*(16*q1 + 9*q2))/40;
        M1 = (L^2*(8*q1 + 7*q2))/120;
        X2 = (L*(b1 + 2*b2))/6;
        Y2 = (L*(4*q1 + 11*q2))/40;
        M2 = 0;
    otherwise       
        error('Tipo de elemento no soportado')
end
feloc = [X1 Y1 M1 X2 Y2 M2]';
end