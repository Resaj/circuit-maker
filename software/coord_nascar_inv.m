% Coordenadas para el circuito nascar invertido

function [dim_cto origen_cto tramos_cto] = coord_nascar_inv()
%% Dimensiones del circuito
X_cto = 3000; % mm
Y_cto = 2000; % mm
dim_cto = [X_cto Y_cto];

%% Punto y direcci�n en origen
x0_pos = 1070;
y0_pos = 430;
x0_dir = -1;
y0_dir = 0;
origen_cto = [x0_pos y0_pos x0_dir y0_dir];

%% Trazado [tipo, longitud]
    % tipo:
        % 0 = recta
        % ang = �ngulo de curva a izquierda en grados
        % -ang = �ngulo de curva a derecha en grados
    % longitud:
        % distancia en mm para recta
        % radio en mm para curva

numero_de_tramos = 4;
tramos_cto = zeros(numero_de_tramos,2);

tramos_cto(1,:) = [-180 570];
tramos_cto(2,:) = [0 860];
tramos_cto(3,:) = [-1 0];
tramos_cto(4,:) = [0 0];
% Los �ltimos dos tramos se dejan con la longitud/radio a 0 y con el �ngulo
% de curva a 1 o -1 para que el programa los calcule de manera �ptima para
% cerrar el circuito.