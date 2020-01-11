%% --------------------------------------------------------------------
% Este programa genera el archivo bmp para un circuito con degradado y pianos
% a partir de las coordenadas significativas del circuito.

% Programa desarrollado por Ruben Espino San Jose
% Equipo Puma Pride de Robotica de Competicion.
%%

clc
clear variables
close all
warning ("off", "Octave:divide-by-zero");
time = clock;
fprintf('\n---------- CIRCUIT MAKER PROGRAM ----------\n');
fprintf('started at %d:%d:%.1f\n', time(4), time(5), time(6));

%% PARAMETROS DE CONFIGURACION
pix_mm = 0.25; % resolucion en pixeles/mm
mm_pix = 1/pix_mm;
blanco = [255 255 255];
negro = [0 0 0];
rojo = [255 0 0];
verde = [15 155 15];
fondo = verde; % color de fondo en RGB
color_piano_1 = rojo;
color_piano_2 = blanco;
color_exterior_pista = 0.15 * 255; % tonalidad exterior de gris en la pista (0 = negro, 255 = blanco)
color_interior_pista = 0.4 * 255; % tonalidad interior de gris en la pista (0 = negro, 255 = blanco)
ancho_pista = 400; % ancho de la carretera en mm
ancho_piano = 60; % ancho de cada piano en mm
marcas_salida = 2; % numero de posiciones de salida = numero de robots
separacion_salida = 500; % separacion entre marcas de salida

representar_trazado_central = 1;
representar_trazado_limite = 0;
generar_circuito = 0;
mostrar_circuito = 0;

%% COORDENADAS DEL CIRCUITO
[dim origen tramos] = coord_nascar();
%[dim origen tramos] = coord_nascar_inv();
%[dim origen tramos] = coord_nascar_vert();
%[dim origen tramos] = coord_mgw2015();
%[dim origen tramos] = coord_mgw2015_inv();
%[dim origen tramos] = coord_cosmobot2017();
%[dim origen tramos] = coord_pruebas2017();
%[dim origen tramos] = coord_alcaniz();

%% Calcular parametros de la trayectoria principal
[m n] = size(tramos);
tramos = [tramos zeros(m,8)];

for i=1:m
    if(tramos(i,1)==0) % recta
        switch(i)
            case 1 % primer tramo (recta)
                tramos(i,3) = origen(3);
                tramos(i,4) = origen(4);
                tramos(i,5) = origen(1);
                tramos(i,6) = origen(2);
                
                xdir = tramos(i,3);
                ydir = tramos(i,4);
                x0 = tramos(i,5);
                y0 = tramos(i,6);
                long = tramos(i,2);
                
                tramos(i,7) = (long/(sqrt(1+(ydir/xdir)^2)) + x0) * ((xdir>=0) - (xdir<0));
                tramos(i,8) = (long/(sqrt(1+(xdir/ydir)^2)) + y0) * ((ydir>=0) - (ydir<0));

            case m-1 % penultimo tramo (recta)
                tramos(i,5) = tramos(i-1,7);
                tramos(i,6) = tramos(i-1,8);

                alpha = atan(origen(4)/origen(3)) + pi*(origen(3)<0);
                alpha = alpha + sum(tramos(1:i,1))*pi/180;
                alpha = alpha + 2*pi*(alpha<=-pi) - 2*pi*(alpha>pi);
                alpha = alpha*(abs(abs(alpha)-pi/2)>1e-3) + (pi/2)*((alpha>0)-(alpha<0))*(abs(abs(alpha)-pi/2)<=1e-3);

                if(alpha~=pi/2 && alpha~=-pi/2)
                    xdir = (abs(alpha)<pi/2) - (abs(alpha)>pi/2);
                    ydir = xdir*tan(alpha);
                    ydir = ydir * (abs(ydir)>1e-10);
                else
                    ydir = (alpha==pi/2) - (alpha==-pi/2);
                    xdir = 0;
                end
                tramos(i,3) = xdir;
                tramos(i,4) = ydir;

                x1 = origen(1);
                y1 = origen(2);
                xdir1 = origen(3);
                ydir1 = origen(4);
                beta1 = atan(ydir1/xdir1) + pi*(xdir1<0);
                x2 = tramos(i,5);
                y2 = tramos(i,6);
                xdir2 = tramos(i,3);
                ydir2 = tramos(i,4);
                beta2 = atan(ydir2/xdir2) + pi*(xdir2<0);
                beta2 = beta2*(abs(abs(beta2)-pi/2)>1e-3) + (pi/2)*((beta2>0)-(beta2<0))*(abs(abs(beta2)-pi/2)<=1e-3);

                if(mod(beta1-beta2,pi)==0)
                    beta3 = beta1 + pi/2;
                    if(beta3~=pi/2 && beta3~=-pi/2)
                        xdir3 = (abs(beta3)<pi/2) - (abs(beta3)>pi/2);
                        ydir3 = xdir3*tan(beta3) * ((xdir3>=0) - (xdir3<0));

                        if(mod(beta3,pi)==0)
                            xf = x2;
                            yf = origen(2);
                        else
                            xf = (y1 - y2 - ydir3/xdir3*x1 + ydir2/xdir2*x2)/(ydir2/xdir2 - ydir3/xdir3);
                            yf = ydir2/xdir2*xf + y2 - ydir2/xdir2*x2;
                        end
                    else
                        xf = origen(1);
                        yf = y2;
                    end
                else
                    A = [ydir2 (-xdir2); ydir1 (-xdir1)];
                    b = [(x2*ydir2 - y2*xdir2); (x1*ydir1 - y1*xdir1)];
                    corte = zeros(1,2);

                    for t=1:2
                        Ab = [A(:,1:t-1),b,A(:,t+1:2)];
                        corte(t) = det(Ab)/det(A);
                    end

                    dist = sqrt((corte(1)-x1)^2 + (corte(2)-y1)^2);
                    if(beta2~=pi/2 && beta2~=-pi/2)
                        a = 1;
                        b = -2*corte(1);
                        c = corte(1)^2 - (dist^2/(1+(ydir2/xdir2)^2));
                        p = [a b c];
                        raices = roots(p);

                        omega = (beta1 - beta2)*(tramos(i+1,1)>0) + (beta2 - beta1)*(tramos(i+1,1)<0);
                        omega = omega + 2*pi*(omega<0);
                        if(omega > pi)
                            xf = min(raices)*(xdir2 < 0) + max(raices)*(xdir2 > 0);
                        else
                            xf = max(raices)*(xdir2 < 0) + min(raices)*(xdir2 > 0);
                        end

                        yf = ydir2/xdir2*xf + y2 - ydir2/xdir2*x2;
                    else
                        xf = x2;
                        yf = corte(2) + dist*((beta2==-pi/2) - (beta2==pi/2));
                    end
                end

                tramos(i,7) = xf;
                tramos(i,8) = yf;
                
                x0 = tramos(i,5);
                y0 = tramos(i,6);
                xf = tramos(i,7);
                yf = tramos(i,8);
                
                long = sqrt((xf-x0)^2 + (yf-y0)^2);
                tramos(i,2) = long;

                gamma1 = beta2 + pi/2*(tramos(i+1,1)<0) - pi/2*(tramos(i+1,1)>0);
                gamma1 = gamma1 + 2*pi*(gamma1<=-pi) - 2*pi*(gamma1>pi);
                gamma1 = gamma1*(abs(abs(gamma1)-pi/2)>1e-3) + (pi/2)*((gamma1>0)-(gamma1<0))*(abs(abs(gamma1)-pi/2)<=1e-3);
                tramos(i+1,9) = gamma1;
                gamma2 = beta1 + pi/2*(tramos(i+1,1)<0) - pi/2*(tramos(i+1,1)>0);
                gamma2 = gamma2 + 2*pi*(gamma2<=-pi) - 2*pi*(gamma2>pi);
                gamma2 = gamma2*(abs(abs(gamma2)-pi/2)>1e-3) + (pi/2)*((gamma2>0)-(gamma2<0))*(abs(abs(gamma2)-pi/2)<=1e-3);
                tramos(i+1,10) = gamma2;

                x0 = tramos(i,7);
                y0 = tramos(i,8);
                xf = origen(1);
                yf = origen(2);

                if(mod(gamma1-gamma2,pi)==0)
                    xc = (x0+xf)/2;
                    yc = (y0+yf)/2;
                else
                    if(gamma1~=pi/2 && gamma1~=-pi/2)
                        xdir1 = (abs(gamma1)<pi/2) - (abs(gamma1)>pi/2);
                        ydir1 = xdir1*tan(gamma1);
                    else
                        ydir1 = (gamma1==pi/2) - (gamma1==-pi/2);
                        xdir1 = 0;
                    end
                    xdir1 = xdir1 * (abs(xdir1)>1e-10);
                    ydir1 = ydir1 * (abs(ydir1)>1e-10);

                    if(gamma2~=pi/2 && gamma2~=-pi/2)
                        xdir2 = (abs(gamma2)<pi/2) - (abs(gamma2)>pi/2);
                        ydir2 = xdir2*tan(gamma2);
                    else
                        ydir2 = (gamma2==pi/2) - (gamma2==-pi/2);
                        xdir2 = 0;
                    end
                    xdir2 = xdir2 * (abs(xdir2)>1e-10);
                    ydir2 = ydir2 * (abs(ydir2)>1e-10);

                    if(xdir1~=0 && xdir2~=0)
                        xc = (y0-yf-x0*ydir1/xdir1+xf*ydir2/xdir2)/(ydir2/xdir2-ydir1/xdir1);
                        yc = (xc-x0)*ydir1/xdir1 + y0;
                    else
                        if(xdir1==0)
                            xc = x0;
                            yc = (xc-xf)*ydir2/xdir2 + yf;
                        else
                            xc = xf;
                            yc = (xc-x0)*ydir1/xdir1 + y0;
                        end
                    end
                end

                tramos(i+1,3) = xc;
                tramos(i+1,4) = yc;

                radio = sqrt((xc-x1)^2 + (yc-y1)^2);
                tramos(i+1,2) = radio;

            case m % ultimo tramo (recta)
                tramos(i,3) = origen(3);
                tramos(i,4) = origen(4);
                tramos(i,7) = origen(1);
                tramos(i,8) = origen(2);
                
                x0 = tramos(i-1,5);
                y0 = tramos(i-1,6);
                beta0 = atan(origen(4)/origen(3)) + pi*(origen(3)<0);
                beta0 = beta0 + sum(tramos(1:i-2,1))*pi/180;
                beta0 = beta0 + 2*pi*(beta0<=-pi);
                beta0 = beta0*(abs(abs(beta0)-pi/2)>1e-3) + (pi/2)*((beta0>0)-(beta0<0))*(abs(abs(beta0)-pi/2)<=1e-3);

                if(beta0~=pi/2 && beta0~=-pi/2)
                    xdir0 = (abs(beta0)<pi/2) - (abs(beta0)>pi/2);
                    ydir0 = xdir0*tan(beta0);
                else
                    ydir0 = (beta0==pi/2) - (beta0==-pi/2);
                    xdir0 = 0;
                end
                xdir0 = xdir0 * (abs(xdir0)>1e-10);
                ydir0 = ydir0 * (abs(ydir0)>1e-10);

                x2 = tramos(i,7);
                y2 = tramos(i,8);
                xdir2 = tramos(i,3);
                ydir2 = tramos(i,4);
                beta2 = atan(ydir2/xdir2) + pi*(xdir2<0);
                beta2 = beta2*(abs(abs(beta2)-pi/2)>1e-3) + (pi/2)*((beta2>0)-(beta2<0))*(abs(abs(beta2)-pi/2)<=1e-3);

                if(mod(beta2-beta0,pi)==0)
                    beta3 = beta0 + pi/2;
                    if(beta3~=pi/2 && beta3~=-pi/2)
                        xdir3 = (abs(beta3)<pi/2) - (abs(beta3)>pi/2);
                        ydir3 = xdir3*tan(beta3);

                        if(mod(beta3,pi)==0)
                            x1 = x2;
                            y1 = y0;
                        else
                            x1 = (y0 - y2 - ydir3/xdir3*x0 + ydir2/xdir2*x2)/(ydir2/xdir2 - ydir3/xdir3);
                            y1 = ydir2/xdir2*x1 + y2 - ydir2/xdir2*x2;
                        end
                    else
                        x1 = x0;
                        y1 = origen(2);
                    end
                else
                    A = [ydir2 (-xdir2); ydir0 (-xdir0)];
                    b = [(x2*ydir2 - y2*xdir2); (x0*ydir0 - y0*xdir0)];
                    corte = zeros(1,2);

                    for t=1:2
                        Ab = [A(:,1:t-1),b,A(:,t+1:2)];
                        corte(t) = det(Ab)/det(A);
                    end

                    dist = sqrt((corte(1)-x0)^2 + (corte(2)-y0)^2);
                    if(beta2~=pi/2 && beta2~=-pi/2)
                        a = 1;
                        b = -2*corte(1);
                        c = corte(1)^2 - (dist^2/(1+(ydir2/xdir2)^2));
                        p = [a b c];
                        raices = roots(p);

                        omega = (beta2 - beta0)*(tramos(i-1,1)>0) + (beta0 - beta2)*(tramos(i-1,1)<0);
                        omega = omega + 2*pi*(omega<0);
                        if(omega > pi)
                            x1 = max(raices)*(xdir2 < 0) + min(raices)*(xdir2 > 0);
                        else
                            x1 = min(raices)*(xdir2 < 0) + max(raices)*(xdir2 > 0);
                        end
                        
                        y1 = ydir2/xdir2*x1 + y2 - ydir2/xdir2*x2;
                    else
                        x1 = x2;
                        y1 = corte(2) + dist*((beta2==pi/2) - (beta2==-pi/2));
                    end
                end
                tramos(i-1,7) = x1;
                tramos(i-1,8) = y1;

                tramos(i,5) = tramos(i-1,7);
                tramos(i,6) = tramos(i-1,8);

                long = sqrt((x2-x1)^2 + (y2-y1)^2);
                tramos(i,2) = long;

                gamma1 = beta0 + pi/2*(tramos(i-1,1)<0) - pi/2*(tramos(i-1,1)>0);
                gamma1 = gamma1 + 2*pi*(gamma1<=-pi) - 2*pi*(gamma1>pi);
                gamma1 = gamma1*(abs(abs(gamma1)-pi/2)>1e-3) + (pi/2)*((gamma1>0)-(gamma1<0))*(abs(abs(gamma1)-pi/2)<=1e-3);
                tramos(i-1,9) = gamma1;
                gamma2 = beta2 + pi/2*(tramos(i-1,1)<0) - pi/2*(tramos(i-1,1)>0);
                gamma2 = gamma2 + 2*pi*(gamma2<=-pi) - 2*pi*(gamma2>pi);
                gamma2 = gamma2*(abs(abs(gamma2)-pi/2)>1e-3) + (pi/2)*((gamma2>0)-(gamma2<0))*(abs(abs(gamma2)-pi/2)<=1e-3);
                tramos(i-1,10) = gamma2;

                if(mod(gamma1-gamma2,pi)==0)
                    xc = (x0+x1)/2;
                    yc = (y0+y1)/2;
                else
                    if(gamma1~=pi/2 && gamma1~=-pi/2)
                        xdir3 = (abs(gamma1)<pi/2) - (abs(gamma1)>pi/2);
                        ydir3 = xdir3*tan(gamma1);
                    else
                        ydir3 = (gamma1==pi/2) - (gamma1==-pi/2);
                        xdir3 = 0;
                    end

                    if(gamma2~=pi/2 && gamma2~=-pi/2)
                        xdir1 = (abs(gamma2)<pi/2) - (abs(gamma2)>pi/2);
                        ydir1 = xdir1*tan(gamma2);
                    else
                        ydir1 = (gamma2==pi/2) - (gamma2==-pi/2);
                        xdir1 = 0;
                    end

                    if(xdir3~=0 && xdir1~=0)
                        xc = (y0-y1-x0*ydir3/xdir3+x1*ydir1/xdir1)/(ydir1/xdir1-ydir3/xdir3);
                        yc = (xc-x0)*ydir3/xdir3 + y0;
                    else
                        if(xdir3==0)
                            xc = x0;
                            yc = (xc-x1)*ydir1/xdir1 + y1;
                        end
                        if(xdir1==0)
                            xc = x1;
                            yc = (xc-x0)*ydir3/xdir3 + y0;
                        end
                    end
                end

                tramos(i-1,3) = xc;
                tramos(i-1,4) = yc;

                radio = sqrt((xc-x0)^2 + (yc-y0)^2);
                tramos(i-1,2) = radio;

            otherwise % otros tramos (recta)
                tramos(i,5) = tramos(i-1,7);
                tramos(i,6) = tramos(i-1,8);
                
                alpha = atan(origen(4)/origen(3)) + pi*(origen(3)<0);
                alpha = alpha + sum(tramos(1:i,1))*pi/180;
                alpha = alpha - 2*pi*(alpha>pi) + 2*pi*(alpha<-pi);
                alpha = alpha*(abs(abs(alpha)-pi/2)>1e-3) + (pi/2)*((alpha>0)-(alpha<0))*(abs(abs(alpha)-pi/2)<=1e-3);
                
                if(alpha~=pi/2 && alpha~=-pi/2)
                    xdir = (abs(alpha)<pi/2) - (abs(alpha)>pi/2);
                    ydir = xdir*tan(alpha);
                else
                    ydir = (alpha==pi/2) - (alpha==-pi/2);
                    xdir = 0;
                end
                tramos(i,3) = xdir;
                tramos(i,4) = ydir;

                x0 = tramos(i,5);
                y0 = tramos(i,6);
                long = tramos(i,2);
                
                tramos(i,7) = long/(sqrt(1+(ydir/xdir)^2)) * ((xdir>=0) - (xdir<0)) + x0;
                tramos(i,8) = long/(sqrt(1+(xdir/ydir)^2)) * ((ydir>=0) - (ydir<0)) + y0;

        end
    else % curva
        if(tramos(i,1)==1 || tramos(i,1)==-1)
            if(sum(tramos(:,1)) - tramos(i,1) > 0)
                tramos(i,1) = 360 - (sum(tramos(:,1)) - tramos(i,1));
            else
                tramos(i,1) = -360 - (sum(tramos(:,1)) - tramos(i,1));
            end
        end
        
        switch(i)
            case 1 % primer tramo (curva)
                tramos(i,5) = origen(1);
                tramos(i,6) = origen(2);

                x0 = tramos(i,5);
                y0 = tramos(i,6);
                radio = tramos(i,2);
                alpha = atan(origen(4)/origen(3)) + pi*(origen(3)<0);
                alpha = alpha + 2*pi*(alpha<=-pi);
                beta = alpha + (tramos(i,1)>0)*pi/2 - (tramos(i,1)<0)*pi/2;
                beta = beta + 2*pi*(beta<=-pi) - 2*pi*(beta>pi);
                beta = beta*(abs(abs(beta)-pi/2)>1e-3) + (pi/2)*((beta>0)-(beta<0))*(abs(abs(beta)-pi/2)<=1e-3);

                if(beta~=pi/2 && beta~=-pi/2)
                    xdir = (beta<pi/2 && beta>-pi/2) - ~(beta<pi/2 && beta>-pi/2);
                    ydir = xdir*tan(beta);

                    a = 1;
                    b = -2*x0;
                    c = x0^2 - (radio^2/(1+(ydir/xdir)^2));
                    p = [a b c];
                    raices = roots(p);
                    xc = raices(1)*(xdir>0 && (x0-raices(1))<0 || xdir<0 && (x0-raices(1))>0) + raices(2)*(xdir>0 && (x0-raices(1))>0 || xdir<0 && (x0-raices(1))<0);
                    yc = ydir/xdir*xc + y0 - ydir/xdir*x0;
                else
                    xc = x0;
                    yc = y0 + radio*((beta==pi/2) - (beta==-pi/2));
                end

                tramos(i,3) = xc;
                tramos(i,4) = yc;

                gamma = beta + tramos(i,1)*pi/180 - pi;
                gamma = gamma + 2*pi*(gamma<=-pi) - 2*pi*(gamma>pi);
                gamma = gamma*(abs(abs(gamma)-pi/2)>1e-3) + (pi/2)*((gamma>0)-(gamma<0))*(abs(abs(gamma)-pi/2)<=1e-3);

                if(gamma~=pi/2 && gamma~=-pi/2)
                    xdir = (gamma<pi/2 && gamma>-pi/2) - ~(gamma<pi/2 && gamma>-pi/2);
                    ydir = xdir*tan(gamma);

                    a = 1;
                    b = -2*xc;
                    c = xc^2 - (radio^2/(1+(ydir/xdir)^2));
                    p = [a b c];
                    raices = roots(p);
                    xf = raices(1)*(xdir>0 && (xc-raices(1))<0 || xdir<0 && (xc-raices(1))>0) + raices(2)*(xdir>0 && (xc-raices(1))>0 || xdir<0 && (xc-raices(1))<0);
                    yf = ydir/xdir*xf + yc - ydir/xdir*xc;
                else
                    xf = xc;
                    yf = yc + radio*((gamma==pi/2) - (gamma==-pi/2));
                end

                tramos(i,7) = xf;
                tramos(i,8) = yf;

                omega = beta + pi;
                omega = omega + 2*pi*(omega<=-pi) - 2*pi*(omega>pi);
                tramos(i,9) = omega;
                tramos(i,10) = gamma;

            case m-1 % penultimo tramo (curva)
                tramos(i,5) = tramos(i-1,7);
                tramos(i,6) = tramos(i-1,8);

                alpha = atan(origen(4)/origen(3)) + pi*(origen(3)<0);
                alpha = alpha + 2*pi*(alpha<=-pi);
                beta = alpha + sum(tramos(1:i-1,1))*pi/180 + (tramos(i,1)<0)*pi/2 - (tramos(i,1)>0)*pi/2;
                beta = beta + 2*pi*(beta<=-pi) - 2*pi*(beta>pi);
                gamma = alpha + (tramos(i,1)<0)*pi/2 - (tramos(i,1)>0)*pi/2;
                gamma = gamma + 2*pi*(gamma<=-pi) - 2*pi*(gamma>pi);
                tramos(i,9) = beta;
                tramos(i,10) = gamma;

            case m % ultimo tramo (curva)
                tramos(i,5) = tramos(i-1,7);
                tramos(i,6) = tramos(i-1,8);
                tramos(i,7) = origen(1);
                tramos(i,8) = origen(2);

            otherwise % otros tramos (curva)
                tramos(i,5) = tramos(i-1,7);
                tramos(i,6) = tramos(i-1,8);
                
                x0 = tramos(i,5);
                y0 = tramos(i,6);
                radio = tramos(i,2);
                alpha = atan(origen(4)/origen(3)) + pi*(origen(3)<0);
                alpha = alpha + sum(tramos(1:i-1,1))*pi/180;
                alpha = alpha + 2*pi*(alpha<=-pi);
                beta = alpha + (tramos(i,1)>0)*pi/2 - (tramos(i,1)<0)*pi/2;
                beta = beta + 2*pi*(beta<=-pi) - 2*pi*(beta>pi);
                beta = beta*(abs(abs(beta)-pi/2)>1e-3) + (pi/2)*((beta>0)-(beta<0))*(abs(abs(beta)-pi/2)<=1e-3);

                if(beta~=pi/2 && beta~=-pi/2)
                    xdir = (beta<pi/2 && beta>-pi/2) - ~(beta<pi/2 && beta>-pi/2);
                    ydir = xdir*tan(beta);

                    a = 1;
                    b = -2*x0;
                    c = x0^2 - (radio^2/(1+(ydir/xdir)^2));
                    p = [a b c];
                    raices = roots(p);
                    xc = raices(1)*(xdir>0 && (x0-raices(1))<0 || xdir<0 && (x0-raices(1))>0) + raices(2)*(xdir>0 && (x0-raices(1))>0 || xdir<0 && (x0-raices(1))<0);
                    yc = ydir/xdir*xc + y0 - ydir/xdir*x0;
                else
                    xc = x0;
                    yc = y0 + radio*((beta==pi/2) - (beta==-pi/2));
                end

                tramos(i,3) = xc;
                tramos(i,4) = yc;

                gamma = beta + tramos(i,1)*pi/180 - pi;
                gamma = gamma + 2*pi*(gamma<=-pi) - 2*pi*(gamma>pi);
                gamma = gamma*(abs(abs(gamma)-pi/2)>1e-3) + (pi/2)*((gamma>0)-(gamma<0))*(abs(abs(gamma)-pi/2)<=1e-3);
                
                if(gamma~=pi/2 && gamma~=-pi/2)
                    xdir = (gamma<pi/2 && gamma>-pi/2) - ~(gamma<pi/2 && gamma>-pi/2);
                    ydir = xdir*tan(gamma);

                    a = 1;
                    b = -2*xc;
                    c = xc^2 - (radio^2/(1+(ydir/xdir)^2));
                    p = [a b c];
                    raices = roots(p);
                    xf = raices(1)*(xdir>0 && (xc-raices(1))<0 || xdir<0 && (xc-raices(1))>0) + raices(2)*(xdir>0 && (xc-raices(1))>0 || xdir<0 && (xc-raices(1))<0);
                    yf = ydir/xdir*xf + yc - ydir/xdir*xc;
                else
                    xf = xc;
                    yf = yc + radio*((gamma==pi/2) - (gamma==-pi/2));
                end

                tramos(i,7) = xf;
                tramos(i,8) = yf;

                omega = beta + pi;
                omega = omega + 2*pi*(omega<=-pi) - 2*pi*(omega>pi);
                tramos(i,9) = omega;
                tramos(i,10) = gamma;
        end
    end
end

%% Representar trazado central
if(representar_trazado_central)
    figure();
    title('Trazado central');
    hold on
    scatter(0, 0, 'filled');
    scatter(dim(1), 0, 'filled');
    scatter(0, dim(2), 'filled');
    scatter(dim(1), dim(2), 'filled');
    line([0 0], [0 dim(2)], 'Color','red','LineStyle','--');
    line([0 dim(1)], [0 0], 'Color','red','LineStyle','--');
    line([dim(1) dim(1)], [0 dim(2)], 'Color','red','LineStyle','--');
    line([0 dim(1)], [dim(2) dim(2)], 'Color','red','LineStyle','--');

    for i=1:m
        scatter(tramos(i,5), tramos(i,6), 'filled');
        if(tramos(i,1)==0) % recta
            line([tramos(i,5) tramos(i,7)], [tramos(i,6) tramos(i,8)]);
        else % curva
            if(tramos(i,1)>0)
                ang1 = tramos(i,9);
                ang2 = tramos(i,10);
            else
                ang1 = tramos(i,10);
                ang2 = tramos(i,9);
            end
            ang2 = ang2 + 2*pi*(ang1>ang2);
            th = ang1:pi/10000:ang2;
            xunit = tramos(i,2) * cos(th) + tramos(i,3);
            yunit = tramos(i,2) * sin(th) + tramos(i,4);
            plot(xunit, yunit);
        end
    end
    hold off
end

%% Calcular curvas y rectas que delimitan la pista
tramos_izq = tramos;
tramos_der = tramos;
tramos_izq_piano = tramos;
tramos_der_piano = tramos;

for i=1:m
    if(i==1 && tramos(1,1)==0 && tramos(m,1)==0) % recta inicial y recta final
        x0 = tramos(1,5);
        y0 = tramos(1,6);
        
        xdir = tramos(1,3);
        ydir = tramos(1,4);
        beta = atan(ydir/xdir) + pi*(xdir<0) + pi/2;
        beta = beta + 2*pi*(beta<=-pi) - 2*pi*(beta>pi);

        tramos_izq(1,5) = cos(beta)*ancho_pista/2 + x0;
        tramos_izq(1,6) = sin(beta)*ancho_pista/2 + y0;
        tramos_der(1,5) = -cos(beta)*ancho_pista/2 + x0;
        tramos_der(1,6) = -sin(beta)*ancho_pista/2 + y0;
        tramos_izq_piano(1,5) = cos(beta)*(ancho_pista/2+ancho_piano) + x0;
        tramos_izq_piano(1,6) = sin(beta)*(ancho_pista/2+ancho_piano) + y0;
        tramos_der_piano(1,5) = -cos(beta)*(ancho_pista/2+ancho_piano) + x0;
        tramos_der_piano(1,6) = -sin(beta)*(ancho_pista/2+ancho_piano) + y0;

        if(tramos(m,1)==0)
            tramos_izq(m,7) = tramos_izq(1,5);
            tramos_izq(m,8) = tramos_izq(1,6);
            tramos_der(m,7) = tramos_der(1,5);
            tramos_der(m,8) = tramos_der(1,6);
            tramos_izq_piano(m,7) = tramos_izq_piano(1,5);
            tramos_izq_piano(m,8) = tramos_izq_piano(1,6);
            tramos_der_piano(m,7) = tramos_der_piano(1,5);
            tramos_der_piano(m,8) = tramos_der_piano(1,6);
        end
    else
        if(tramos(i,1)~=0) % curva
            tramos_izq(i,2) = tramos(i,2) - ancho_pista/2*((tramos(i,1)>0) - (tramos(i,1)<0));
            tramos_der(i,2) = tramos(i,2) + ancho_pista/2*((tramos(i,1)>0) - (tramos(i,1)<0));
            tramos_izq_piano(i,2) = tramos(i,2) - (ancho_pista/2 + ancho_piano)*((tramos(i,1)>0) - (tramos(i,1)<0));
            tramos_der_piano(i,2) = tramos(i,2) + (ancho_pista/2 + ancho_piano)*((tramos(i,1)>0) - (tramos(i,1)<0));
            
            x0 = tramos(i,5);
            y0 = tramos(i,6);
            alpha = tramos(i,9);
            beta = alpha + pi*(tramos(i,1)>0);
            beta = beta + 2*pi*(beta<=-pi) - 2*pi*(beta>pi);

            tramos_izq(i,5) = cos(beta)*ancho_pista/2 + x0;
            tramos_izq(i,6) = sin(beta)*ancho_pista/2 + y0;
            tramos_der(i,5) = -cos(beta)*ancho_pista/2 + x0;
            tramos_der(i,6) = -sin(beta)*ancho_pista/2 + y0;
            tramos_izq_piano(i,5) = cos(beta)*(ancho_pista/2+ancho_piano) + x0;
            tramos_izq_piano(i,6) = sin(beta)*(ancho_pista/2+ancho_piano) + y0;
            tramos_der_piano(i,5) = -cos(beta)*(ancho_pista/2+ancho_piano) + x0;
            tramos_der_piano(i,6) = -sin(beta)*(ancho_pista/2+ancho_piano) + y0;

            j = (i-1)*(i>1) + m*(i==1);
            if(tramos(j,1)==0)
                tramos_izq(j,7) = tramos_izq(i,5);
                tramos_izq(j,8) = tramos_izq(i,6);
                tramos_der(j,7) = tramos_der(i,5);
                tramos_der(j,8) = tramos_der(i,6);
                tramos_izq_piano(j,7) = tramos_izq_piano(i,5);
                tramos_izq_piano(j,8) = tramos_izq_piano(i,6);
                tramos_der_piano(j,7) = tramos_der_piano(i,5);
                tramos_der_piano(j,8) = tramos_der_piano(i,6);
            end

            xf = tramos(i,7);
            yf = tramos(i,8);
            alpha = tramos(i,10);
            beta = alpha + pi*(tramos(i,1)>0);
            beta = beta + 2*pi*(beta<=-pi) - 2*pi*(beta>pi);

            tramos_izq(i,7) = cos(beta)*ancho_pista/2 + xf;
            tramos_izq(i,8) = sin(beta)*ancho_pista/2 + yf;
            tramos_der(i,7) = -cos(beta)*ancho_pista/2 + xf;
            tramos_der(i,8) = -sin(beta)*ancho_pista/2 + yf;
            tramos_izq_piano(i,7) = cos(beta)*(ancho_pista/2+ancho_piano) + xf;
            tramos_izq_piano(i,8) = sin(beta)*(ancho_pista/2+ancho_piano) + yf;
            tramos_der_piano(i,7) = -cos(beta)*(ancho_pista/2+ancho_piano) + xf;
            tramos_der_piano(i,8) = -sin(beta)*(ancho_pista/2+ancho_piano) + yf;

            j = (i+1)*(i<m) + (i==m);
            if(tramos(j,1)==0)
                tramos_izq(j,5) = tramos_izq(i,7);
                tramos_izq(j,6) = tramos_izq(i,8);
                tramos_der(j,5) = tramos_der(i,7);
                tramos_der(j,6) = tramos_der(i,8);
                tramos_izq_piano(j,5) = tramos_izq_piano(i,7);
                tramos_izq_piano(j,6) = tramos_izq_piano(i,8);
                tramos_der_piano(j,5) = tramos_der_piano(i,7);
                tramos_der_piano(j,6) = tramos_der_piano(i,8);
            end
        end
    end
end

%% Representar trazado limite del circuito
if(representar_trazado_limite)
    figure();
    title('Trazado limite');
    hold on
    scatter(0, 0, 'filled');
    scatter(dim(1), 0, 'filled');
    scatter(0, dim(2), 'filled');
    scatter(dim(1), dim(2), 'filled');
    line([0 0], [0 dim(2)], 'Color','red','LineStyle','--');
    line([0 dim(1)], [0 0], 'Color','red','LineStyle','--');
    line([dim(1) dim(1)], [0 dim(2)], 'Color','red','LineStyle','--');
    line([0 dim(1)], [dim(2) dim(2)], 'Color','red','LineStyle','--');

    for i=1:m
        scatter(tramos_izq(i,5), tramos_izq(i,6), 'filled');
        scatter(tramos_der(i,5), tramos_der(i,6), 'filled');
        scatter(tramos_izq_piano(i,5), tramos_izq_piano(i,6), 'filled');
        scatter(tramos_der_piano(i,5), tramos_der_piano(i,6), 'filled');

        if(tramos(i,1)==0) % recta
            line([tramos_izq(i,5) tramos_izq(i,7)], [tramos_izq(i,6) tramos_izq(i,8)]);
            line([tramos_der(i,5) tramos_der(i,7)], [tramos_der(i,6) tramos_der(i,8)]);
            line([tramos_izq_piano(i,5) tramos_izq_piano(i,7)], [tramos_izq_piano(i,6) tramos_izq_piano(i,8)]);
            line([tramos_der_piano(i,5) tramos_der_piano(i,7)], [tramos_der_piano(i,6) tramos_der_piano(i,8)]);
            
        else % curva
            if(tramos(i,1)>0)
                ang1 = tramos(i,9);
                ang2 = tramos(i,10);
            else
                ang1 = tramos(i,10);
                ang2 = tramos(i,9);
            end
            ang2 = ang2 + 2*pi*(ang1>ang2);
            th = ang1:pi/1000:ang2;
            
            xunit = tramos_izq(i,2) * cos(th) + tramos_izq(i,3);
            yunit = tramos_izq(i,2) * sin(th) + tramos_izq(i,4);
            plot(xunit, yunit);
            xunit = tramos_der(i,2) * cos(th) + tramos_der(i,3);
            yunit = tramos_der(i,2) * sin(th) + tramos_der(i,4);
            plot(xunit, yunit);
            xunit = tramos_izq_piano(i,2) * cos(th) + tramos_izq_piano(i,3);
            yunit = tramos_izq_piano(i,2) * sin(th) + tramos_izq_piano(i,4);
            plot(xunit, yunit);
            xunit = tramos_der_piano(i,2) * cos(th) + tramos_der_piano(i,3);
            yunit = tramos_der_piano(i,2) * sin(th) + tramos_der_piano(i,4);
            plot(xunit, yunit);
        end
    end
    hold off
end

%% GENERAR CIRCUITO
if(generar_circuito)
%% Crear matriz de la imagen
circuito = uint8(zeros(dim(2)/mm_pix, dim(1)/mm_pix, 3));
R = uint8(fondo(1)*ones(dim(2)/mm_pix, dim(1)/mm_pix));
G = uint8(fondo(2)*ones(dim(2)/mm_pix, dim(1)/mm_pix));
B = uint8(fondo(3)*ones(dim(2)/mm_pix, dim(1)/mm_pix));

%% Pintar pista con degradado
lineas = 2*ancho_pista/mm_pix;

for i=1:m
    if(tramos(i,1)==0) % recta
        xi0 = tramos_izq(i,5);
        yi0 = tramos_izq(i,6);
        xd0 = tramos_der(i,5);
        yd0 = tramos_der(i,6);
        xif = tramos_izq(i,7);
        yif = tramos_izq(i,8);
        xdf = tramos_der(i,7);
        ydf = tramos_der(i,8);
        xdir = tramos(i,3);
        ydir = tramos(i,4);
        modulo = tramos(i,2);
        
        for j=1:lineas
            x0 = ((lineas-j)*xi0 + (j-1)*xd0)/(lineas-1);
            y0 = ((lineas-j)*yi0 + (j-1)*yd0)/(lineas-1);
            xf = ((lineas-j)*xif + (j-1)*xdf)/(lineas-1);
            yf = ((lineas-j)*yif + (j-1)*ydf)/(lineas-1);
            
            if(abs(x0-xf)>abs(y0-yf))
                xunit = x0*(x0<xf) + xf*(x0>=xf):max(mm_pix/2,mm_pix/2*(xf-x0)/modulo):xf*(x0<xf) + x0*(x0>=xf);
                yunit = (xunit-x0)*ydir/xdir + y0;
            else
                yunit = y0*(y0<yf) + yf*(y0>=yf):max(mm_pix/2,mm_pix/2*(yf-y0)/modulo):yf*(y0<yf) + y0*(y0>=yf);
                xunit = (yunit-y0)*xdir/ydir + x0;
            end

            if(abs(xd0-xi0)>abs(yd0-yi0))
                color = (2*(color_exterior_pista-color_interior_pista)/(xd0-xi0)*(x0-xd0)+color_exterior_pista)*(x0<(xi0+xd0)/2 && xd0<xi0 || x0>=(xi0+xd0)/2 && xd0>xi0) + (2*(color_interior_pista-color_exterior_pista)/(xd0-xi0)*(x0-xi0)+color_exterior_pista)*(x0>=(xi0+xd0)/2 && xd0<xi0 || x0<(xi0+xd0)/2 && xd0>xi0);
            else
                color = (2*(color_exterior_pista-color_interior_pista)/(yd0-yi0)*(y0-yd0)+color_exterior_pista)*(y0<(yi0+yd0)/2 && yd0<yi0 || y0>=(yi0+yd0)/2 && yd0>yi0) + (2*(color_interior_pista-color_exterior_pista)/(yd0-yi0)*(y0-yi0)+color_exterior_pista)*(y0>=(yi0+yd0)/2 && yd0<yi0 || y0<(yi0+yd0)/2 && yd0>yi0);
            end
            
            for q=1:size(xunit,2)
                R(size(circuito,1)-floor(yunit(q)/mm_pix),floor(xunit(q)/mm_pix)) = color;
                G(size(circuito,1)-floor(yunit(q)/mm_pix),floor(xunit(q)/mm_pix)) = color;
                B(size(circuito,1)-floor(yunit(q)/mm_pix),floor(xunit(q)/mm_pix)) = color;
            end
        end
            
    else % curva
        ri = tramos_izq(i,2);
        rd = tramos_der(i,2);
        xc = tramos(i,3);
        yc = tramos(i,4);
        modulo = pi*tramos(i,2);
        
        if(tramos(i,1)>0)
            ang1 = tramos(i,9);
            ang2 = tramos(i,10);
        else
            ang1 = tramos(i,10);
            ang2 = tramos(i,9);
        end
        ang2 = ang2 + 2*pi*(ang1>ang2);

        for j=1:lineas
            r = ((lineas-j)*(ri*(ri<rd) + rd*(ri>rd)) + (j-1)*(rd*(ri<rd) + ri*(ri>rd)))/(lineas-1);
            
            th = ang1:pi*mm_pix/2/modulo:ang2;
            xunit = r * cos(th) + xc;
            yunit = r * sin(th) + yc;

            color = (2*(color_exterior_pista-color_interior_pista)/(rd-ri)*(r-rd)+color_exterior_pista)*(r<(ri+rd)/2 && rd<ri || r>=(ri+rd)/2 && rd>ri) + (2*(color_interior_pista-color_exterior_pista)/(rd-ri)*(r-ri)+color_exterior_pista)*(r>=(ri+rd)/2 && rd<ri || r<(ri+rd)/2 && rd>ri);

            for q=1:size(xunit,2)
                R(size(circuito,1)-floor(yunit(q)/mm_pix),floor(xunit(q)/mm_pix)) = color;
                G(size(circuito,1)-floor(yunit(q)/mm_pix),floor(xunit(q)/mm_pix)) = color;
                B(size(circuito,1)-floor(yunit(q)/mm_pix),floor(xunit(q)/mm_pix)) = color;
            end
        end
    end
end

%% Pintar pianos
lineas = 2*ancho_piano/mm_pix;
dist_anterior_izq = 0;
dist_anterior_der = 0;
long_izq = 0;
long_der = 0;
for i=1:m
    if(tramos(i,1)==0)
        long_izq = long_izq + tramos(i,2);
        long_der = long_der + tramos(i,2);
    else
        long_izq = long_izq + tramos_izq(i,2)*(abs(tramos(i,1)))*pi/180;
        long_der = long_der + tramos_der(i,2)*(abs(tramos(i,1)))*pi/180;
    end
end
num_tramos_piano_izq = round(long_izq/(2*ancho_piano));
num_tramos_piano_der = round(long_der/(2*ancho_piano));
long_piano_izq = long_izq/(2*num_tramos_piano_izq);
long_piano_der = long_der/(2*num_tramos_piano_der);

for i=1:m
    if(tramos(i,1)==0) % recta
        xi0 = tramos_izq(i,5);
        yi0 = tramos_izq(i,6);
        xpi0 = tramos_izq_piano(i,5);
        ypi0 = tramos_izq_piano(i,6);
        xd0 = tramos_der(i,5);
        yd0 = tramos_der(i,6);
        xpd0 = tramos_der_piano(i,5);
        ypd0 = tramos_der_piano(i,6);
        xif = tramos_izq(i,7);
        yif = tramos_izq(i,8);
        xpif = tramos_izq_piano(i,7);
        ypif = tramos_izq_piano(i,8);
        xdf = tramos_der(i,7);
        ydf = tramos_der(i,8);
        xpdf = tramos_der_piano(i,7);
        ypdf = tramos_der_piano(i,8);
        xdir = tramos(i,3);
        ydir = tramos(i,4);
        modulo = tramos(i,2);
        modulo_izq = modulo;
        modulo_der = modulo;
        
        for j=1:lineas
            x0 = ((lineas-j)*xpi0 + (j-1)*xi0)/(lineas-1);
            y0 = ((lineas-j)*ypi0 + (j-1)*yi0)/(lineas-1);
            xf = ((lineas-j)*xpif + (j-1)*xif)/(lineas-1);
            yf = ((lineas-j)*ypif + (j-1)*yif)/(lineas-1);
            
            if(abs(x0-xf)>abs(y0-yf))
                xunit = x0*(x0<xf) + xf*(x0>=xf):max(mm_pix/2,mm_pix/2*(xf-x0)/modulo):xf*(x0<xf) + x0*(x0>=xf);
                yunit = (xunit-x0)*ydir/xdir + y0;
            else
                yunit = y0*(y0<yf) + yf*(y0>=yf):max(mm_pix/2,mm_pix/2*(yf-y0)/modulo):yf*(y0<yf) + y0*(y0>=yf);
                if(x0~=xf)
                    xunit = (yunit-y0)*xdir/ydir + x0;
                else
                    xunit = ones(size(yunit))*x0;
                end
            end

            for q=1:size(xunit,2)
                t = q*(tramos(i,5)<=tramos(i,7)) + (size(xunit,2)-q+1)*(tramos(i,5)>tramos(i,7));
                
                dist_tramo_actual_izq = sqrt((x0-xunit(t))^2 + (y0-yunit(t))^2);
                dist = dist_tramo_actual_izq + dist_anterior_izq;

                if((dist/(2*long_piano_izq)-floor(dist/(2*long_piano_izq)))<0.5)
                    color = color_piano_1;
                else
                    color = color_piano_2;
                end
                
                R(size(circuito,1)-floor(yunit(t)/mm_pix),floor(xunit(t)/mm_pix)) = color(1);
                G(size(circuito,1)-floor(yunit(t)/mm_pix),floor(xunit(t)/mm_pix)) = color(2);
                B(size(circuito,1)-floor(yunit(t)/mm_pix),floor(xunit(t)/mm_pix)) = color(3);
            end
            
            x0 = ((lineas-j)*xd0 + (j-1)*xpd0)/(lineas-1);
            y0 = ((lineas-j)*yd0 + (j-1)*ypd0)/(lineas-1);
            xf = ((lineas-j)*xdf + (j-1)*xpdf)/(lineas-1);
            yf = ((lineas-j)*ydf + (j-1)*ypdf)/(lineas-1);
            
            if(abs(x0-xf)>abs(y0-yf))
                xunit = x0*(x0<xf) + xf*(x0>=xf):max(mm_pix/2,mm_pix/2*(xf-x0)/modulo):xf*(x0<xf) + x0*(x0>=xf);
                yunit = (xunit-x0)*ydir/xdir + y0;
            else
                yunit = y0*(y0<yf) + yf*(y0>=yf):max(mm_pix/2,mm_pix/2*(yf-y0)/modulo):yf*(y0<yf) + y0*(y0>=yf);
                if(x0~=xf)
                    xunit = (yunit-y0)*xdir/ydir + x0;
                else
                    xunit = ones(size(yunit))*x0;
                end
            end

            for q=1:size(xunit,2)
                t = q*(tramos(i,5)<=tramos(i,7)) + (size(xunit,2)-q+1)*(tramos(i,5)>tramos(i,7));
                
                dist_tramo_actual_der = sqrt((x0-xunit(t))^2 + (y0-yunit(t))^2);
                dist = dist_tramo_actual_der + dist_anterior_der;

                if((dist/(2*long_piano_der)-floor(dist/(2*long_piano_der)))<0.5)
                    color = color_piano_1;
                else
                    color = color_piano_2;
                end
                
                R(size(circuito,1)-floor(yunit(t)/mm_pix),floor(xunit(t)/mm_pix)) = color(1);
                G(size(circuito,1)-floor(yunit(t)/mm_pix),floor(xunit(t)/mm_pix)) = color(2);
                B(size(circuito,1)-floor(yunit(t)/mm_pix),floor(xunit(t)/mm_pix)) = color(3);
            end
        end
            
    else % curva
        ri = tramos_izq(i,2);
        rpi = tramos_izq_piano(i,2);
        rd = tramos_der(i,2);
        rpd = tramos_der_piano(i,2);
        xc = tramos(i,3);
        yc = tramos(i,4);
        modulo = abs(tramos(i,2)*tramos(i,1)*pi/180);
        modulo_izq = abs(tramos_izq(i,2)*tramos_izq(i,1)*pi/180);
        modulo_der = abs(tramos_der(i,2)*tramos_der(i,1)*pi/180);
        
        if(tramos(i,1)>0)
            ang1 = tramos(i,9);
            ang2 = tramos(i,10);
        else
            ang1 = tramos(i,10);
            ang2 = tramos(i,9);
        end
        ang2 = ang2 + 2*pi*(ang1>ang2);
        
        for j=1:lineas
            r = ((lineas-j)*rpi + (j-1)*ri)/(lineas-1);
            
            th = ang1:mm_pix/4/modulo:ang2;
            xunit = r * cos(th) + xc;
            yunit = r * sin(th) + yc;

            for q=1:size(xunit,2)
                t = q*(tramos(i,1)>0) + (size(xunit,2)-q+1)*(tramos(i,1)<0);

                dist_tramo_actual_izq = ((th(t)-ang1)*(tramos(i,1)>0) + (ang2-th(t))*(tramos(i,1)<0))*ri;
                dist = dist_tramo_actual_izq + dist_anterior_izq;
                
                if((dist/(2*long_piano_izq)-floor(dist/(2*long_piano_izq)))<0.5)
                    color = color_piano_1;
                else
                    color = color_piano_2;
                end
                
                R(size(circuito,1)-floor(yunit(t)/mm_pix),floor(xunit(t)/mm_pix)) = color(1);
                G(size(circuito,1)-floor(yunit(t)/mm_pix),floor(xunit(t)/mm_pix)) = color(2);
                B(size(circuito,1)-floor(yunit(t)/mm_pix),floor(xunit(t)/mm_pix)) = color(3);
            end

            r = ((lineas-j)*rd + (j-1)*rpd)/(lineas-1);
            
            th = ang1:mm_pix/4/modulo:ang2;
            xunit = r * cos(th) + xc;
            yunit = r * sin(th) + yc;

            for q=1:size(xunit,2)
                t = q*(tramos(i,1)>0) + (size(xunit,2)-q+1)*(tramos(i,1)<0);

                dist_tramo_actual_der = ((th(t)-ang1)*(tramos(i,1)>0) + (ang2-th(t))*(tramos(i,1)<0))*rd;
                dist = dist_tramo_actual_der + dist_anterior_der;

                if((dist/(2*long_piano_der)-floor(dist/(2*long_piano_der)))<0.5)
                    color = color_piano_1;
                else
                    color = color_piano_2;
                end

                R(size(circuito,1)-floor(yunit(t)/mm_pix),floor(xunit(t)/mm_pix)) = color(1);
                G(size(circuito,1)-floor(yunit(t)/mm_pix),floor(xunit(t)/mm_pix)) = color(2);
                B(size(circuito,1)-floor(yunit(t)/mm_pix),floor(xunit(t)/mm_pix)) = color(3);
            end
        end
    end
    
    dist_anterior_izq = dist_anterior_izq + modulo_izq;
    dist_anterior_der = dist_anterior_der + modulo_der;
end

%% Pintar marcas de salida y l�nea de meta
i = 1;
while (tramos(i,1)~=0)
    i = i+1;
end
recta_meta = i;

xi0 = tramos_izq(i,5);
yi0 = tramos_izq(i,6);
xpi0 = tramos_izq_piano(i,5);
ypi0 = tramos_izq_piano(i,6);
xd0 = tramos_der(i,5);
yd0 = tramos_der(i,6);
xpd0 = tramos_der_piano(i,5);
ypd0 = tramos_der_piano(i,6);
xif = tramos_izq(i,7);
yif = tramos_izq(i,8);
xpif = tramos_izq_piano(i,7);
ypif = tramos_izq_piano(i,8);
xdf = tramos_der(i,7);
ydf = tramos_der(i,8);
xpdf = tramos_der_piano(i,7);
ypdf = tramos_der_piano(i,8);
xdir = tramos(i,3);
ydir = tramos(i,4);
modulo = tramos(i,2);
lineas = 2*modulo/mm_pix;
vuelta = sum(tramos(:,1));

xei0 = xpi0 - (xi0 - xpi0);
yei0 = ypi0 - (yi0 - ypi0);
xed0 = xpd0 - (xd0 - xpd0);
yed0 = ypd0 - (yd0 - ypd0);
xeif = xpif - (xif - xpif);
yeif = ypif - (yif - ypif);
xedf = xpdf - (xdf - xpdf);
yedf = ypdf - (ydf - ypdf);

for j=1:lineas
    x0 = ((lineas-j)*xpi0 + (j-1)*xpif)/(lineas-1);
    y0 = ((lineas-j)*ypi0 + (j-1)*ypif)/(lineas-1);
    xf = ((lineas-j)*xei0 + (j-1)*xeif)/(lineas-1);
    yf = ((lineas-j)*yei0 + (j-1)*yeif)/(lineas-1);

    if(abs(x0-xf)>abs(y0-yf))
        xunit = x0*(x0<xf) + xf*(x0>=xf):max(mm_pix/2,mm_pix/2*(xf-x0)/modulo):xf*(x0<xf) + x0*(x0>=xf);
        yunit = -(xunit-x0)*xdir/ydir + y0;
    else
        yunit = y0*(y0<yf) + yf*(y0>=yf):max(mm_pix/2,mm_pix/2*(yf-y0)/modulo):yf*(y0<yf) + y0*(y0>=yf);
        if(x0~=xf)
            xunit = -(yunit-y0)*ydir/xdir + x0;
        else
            xunit = ones(size(yunit))*x0;
        end
    end

    dist = sqrt((x0-xpi0)^2 + (y0-ypi0)^2);

    color = fondo;
    marca = 1;
    while(isequal(color,negro)==0 && marca<=marcas_salida)
        if(vuelta>0) % sentido antihorario
            if(mod(marca,2)==1 && dist>(100+(marcas_salida-marca)*separacion_salida) && dist<(100+(marcas_salida-marca)*separacion_salida+20))
                color = negro;
            end
        else % sentido horario
            if(mod(marca,2)==0 && dist>(100+(marcas_salida-marca)*separacion_salida) && dist<(100+(marcas_salida-marca)*separacion_salida+20))
                color = negro;
            end
        end
        marca = marca + 1;
    end
    
    if(isequal(color,negro)==0)
        if(dist>(100+(marcas_salida-1)*separacion_salida+200) && dist<=(100+(marcas_salida-1)*separacion_salida+200+sqrt((x0-xf)^2+(y0-yf)^2)))
            for q=1:size(xunit,2)
                if(q>size(xunit,2)/3 && q<=size(xunit,2)*2/3)
                    if(dist>(100+(marcas_salida-1)*separacion_salida+200+sqrt((x0-xf)^2+(y0-yf)^2)/3) && dist<=(100+(marcas_salida-1)*separacion_salida+200+sqrt((x0-xf)^2+(y0-yf)^2)*2/3))
                        color = negro;
                    else
                        color = blanco;
                    end
                else
                    if(dist>(100+(marcas_salida-1)*separacion_salida+200+sqrt((x0-xf)^2+(y0-yf)^2)/3) && dist<=(100+(marcas_salida-1)*separacion_salida+200+sqrt((x0-xf)^2+(y0-yf)^2)*2/3))
                        color = blanco;
                    else
                        color = negro;
                    end
                end
                
                R(size(circuito,1)-floor(yunit(q)/mm_pix),floor(xunit(q)/mm_pix)) = color(1);
                G(size(circuito,1)-floor(yunit(q)/mm_pix),floor(xunit(q)/mm_pix)) = color(2);
                B(size(circuito,1)-floor(yunit(q)/mm_pix),floor(xunit(q)/mm_pix)) = color(3);
            end
        end
    else
        for q=1:size(xunit,2)
            R(size(circuito,1)-floor(yunit(q)/mm_pix),floor(xunit(q)/mm_pix)) = color(1);
            G(size(circuito,1)-floor(yunit(q)/mm_pix),floor(xunit(q)/mm_pix)) = color(2);
            B(size(circuito,1)-floor(yunit(q)/mm_pix),floor(xunit(q)/mm_pix)) = color(3);
        end
    end

    x0 = ((lineas-j)*xed0 + (j-1)*xedf)/(lineas-1);
    y0 = ((lineas-j)*yed0 + (j-1)*yedf)/(lineas-1);
    xf = ((lineas-j)*xpd0 + (j-1)*xpdf)/(lineas-1);
    yf = ((lineas-j)*ypd0 + (j-1)*ypdf)/(lineas-1);

    if(abs(x0-xf)>abs(y0-yf))
        xunit = x0*(x0<xf) + xf*(x0>=xf):max(mm_pix/2,mm_pix/2*(xf-x0)/modulo):xf*(x0<xf) + x0*(x0>=xf);
        yunit = -(xunit-x0)*xdir/ydir + y0;
    else
        yunit = y0*(y0<yf) + yf*(y0>=yf):max(mm_pix/2,mm_pix/2*(yf-y0)/modulo):yf*(y0<yf) + y0*(y0>=yf);
        if(x0~=xf)
            xunit = -(yunit-y0)*ydir/xdir + x0;
        else
            xunit = ones(size(yunit))*x0;
        end
    end

    dist = sqrt((x0-xed0)^2 + (y0-yed0)^2);

    color = fondo;
    marca = 1;
    while(isequal(color,negro)==0 && marca<=marcas_salida)
        if(vuelta>0) % sentido antihorario
            if(mod(marca,2)==0 && dist>(100+(marcas_salida-marca)*separacion_salida) && dist<(100+(marcas_salida-marca)*separacion_salida+20))
                color = negro;
            end
        else % sentido horario
            if(mod(marca,2)==1 && dist>(100+(marcas_salida-marca)*separacion_salida) && dist<(100+(marcas_salida-marca)*separacion_salida+20))
                color = negro;
            end
        end
        marca = marca + 1;
    end
    
    if(isequal(color,negro)==0)
        if(dist>(100+(marcas_salida-1)*separacion_salida+200) && dist<(100+(marcas_salida-1)*separacion_salida+200+sqrt((x0-xf)^2+(y0-yf)^2)))
            for q=1:size(xunit,2)
                if(q>size(xunit,2)/3 && q<=size(xunit,2)*2/3)
                    if(dist>(100+(marcas_salida-1)*separacion_salida+200+sqrt((x0-xf)^2+(y0-yf)^2)/3) && dist<=(100+(marcas_salida-1)*separacion_salida+200+sqrt((x0-xf)^2+(y0-yf)^2)*2/3))
                        color = negro;
                    else
                        color = blanco;
                    end
                else
                    if(dist>(100+(marcas_salida-1)*separacion_salida+200+sqrt((x0-xf)^2+(y0-yf)^2)/3) && dist<=(100+(marcas_salida-1)*separacion_salida+200+sqrt((x0-xf)^2+(y0-yf)^2)*2/3))
                        color = blanco;
                    else
                        color = negro;
                    end
                end

                R(size(circuito,1)-floor(yunit(q)/mm_pix),floor(xunit(q)/mm_pix)) = color(1);
                G(size(circuito,1)-floor(yunit(q)/mm_pix),floor(xunit(q)/mm_pix)) = color(2);
                B(size(circuito,1)-floor(yunit(q)/mm_pix),floor(xunit(q)/mm_pix)) = color(3);
            end
        end
    else
        for q=1:size(xunit,2)
            R(size(circuito,1)-floor(yunit(q)/mm_pix),floor(xunit(q)/mm_pix)) = color(1);
            G(size(circuito,1)-floor(yunit(q)/mm_pix),floor(xunit(q)/mm_pix)) = color(2);
            B(size(circuito,1)-floor(yunit(q)/mm_pix),floor(xunit(q)/mm_pix)) = color(3);
        end
    end
end

%% Exportar circuito
circuito(:,:,1) = uint8(R);
circuito(:,:,2) = uint8(G);
circuito(:,:,3) = uint8(B);

imwrite(circuito,'circuito.bmp');

%% Mostrar circuito
if(mostrar_circuito)
    figure();
    title('Circuito');
    imshow(circuito);
end

%% END GENERAR CIRCUITO
end

time = clock;
%fprintf('\nfinished at %d:%d:%2.0f\n', time(4), time(5), time(6));
%input('Pulsa enter para finalizar el programa', 's');
