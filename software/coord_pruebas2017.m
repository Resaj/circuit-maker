% Coordenadas para el circuito MGW 2015

function [dim_cto origen_cto tramos_cto] = coord_pruebas2017()
%% Dimensiones del circuito
X_cto = 2840; % mm
Y_cto = 1600; % mm
dim_cto = [X_cto Y_cto];

%% Punto y dirección en origen
x0_pos = 770;
y0_pos = 320;
x0_dir = 1;
y0_dir = 0;
origen_cto = [x0_pos y0_pos x0_dir y0_dir];

%% Trazado [tipo, longitud]
    % tipo:
        % 0 = recta
        % -ang = ángulo de curva a izquierda en grados
        % ang = ángulo de curva a derecha en grados
    % longitud:
        % distancia en mm para recta
        % radio en mm para curva

numero_de_tramos = 6;
tramos_cto = zeros(numero_de_tramos,2);

tramos_cto(1,:) = [0 1300];
tramos_cto(2,:) = [210 500];
tramos_cto(3,:) = [0 181.85];
tramos_cto(4,:) = [-60 485];
tramos_cto(5,:) = [0 0];
tramos_cto(6,:) = [1 0];
% Los últimos dos tramos se dejan con la longitud/radio a 0 y con el ángulo
% de curva a 1 o -1 para que el programa los calcule de manera óptima para
% cerrar el circuito.