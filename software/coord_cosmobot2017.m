% Coordenadas para el circuito de Cosmobot 2017

function [dim_cto origen_cto tramos_cto] = coord_cosmobot2017()
%% Dimensiones del circuito
X_cto = 7000; % mm
Y_cto = 2000; % mm
dim_cto = [X_cto Y_cto];

%% Punto y direcci�n en origen
x0_pos = 1000;
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

numero_de_tramos = 10;
tramos_cto = zeros(numero_de_tramos,2);

tramos_cto(1,:) = [0 5000];
tramos_cto(2,:) = [190 600];
tramos_cto(3,:) = [0 2400];
tramos_cto(4,:) = [-30 570];
tramos_cto(5,:) = [0 1100];
tramos_cto(6,:) = [20 570];
tramos_cto(7,:) = [0 1010];
tramos_cto(8,:) = [90 570];
tramos_cto(9,:) = [0 0];
tramos_cto(10,:) = [1 0];
% Los �ltimos dos tramos se dejan con la longitud/radio a 0 y con el �ngulo
% de curva a 1 o -1 para que el programa los calcule de manera �ptima para
% cerrar el circuito.