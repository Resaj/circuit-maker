% Coordenadas para el circuito MGW 2015

function [dim_cto origen_cto tramos_cto] = coord_mge2017()
%% Dimensiones del circuito
X_cto = 6500; % mm
Y_cto = 2000; % mm
dim_cto = [X_cto Y_cto];

%% Punto y direcci�n en origen
x0_pos = 1130;
y0_pos = 400;
x0_dir = 1;
y0_dir = 0;
origen_cto = [x0_pos y0_pos x0_dir y0_dir];

%% Trazado [tipo, longitud]
    % tipo:
        % 0 = recta
        % -ang = �ngulo de curva a izquierda en grados
        % ang = �ngulo de curva a derecha en grados
    % longitud:
        % distancia en mm para recta
        % radio en mm para curva

numero_de_tramos = 6;
tramos_cto = zeros(numero_de_tramos,2);

tramos_cto(1,:) = [0 4240];
tramos_cto(2,:) = [195 600];
tramos_cto(3,:) = [0 1873];
tramos_cto(4,:) = [-30 600];
tramos_cto(5,:) = [0 0];
tramos_cto(6,:) = [1 0];
% Los �ltimos dos tramos se dejan con la longitud/radio a 0 y con el �ngulo
% de curva a 1 o -1 para que el programa los calcule de manera �ptima para
% cerrar el circuito.