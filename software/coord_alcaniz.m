% Coordenadas para el circuito de Alcañiz

function [dim_cto origen_cto tramos_cto] = coord_alcaniz()
%% Dimensiones del circuito
X_cto = 7000; % mm
Y_cto = 2000; % mm
dim_cto = [X_cto Y_cto];

%% Punto y direccion en origen
x0_pos = 1590;
y0_pos = 320;
x0_dir = 1;
y0_dir = 0;
origen_cto = [x0_pos y0_pos x0_dir y0_dir];

%% Trazado [tipo, longitud]
    % tipo:
        % 0 = recta
        % -ang = angulo de curva a izquierda en grados
        % ang = angulo de curva a derecha en grados
    % longitud:
        % distancia en mm para recta
        % radio en mm para curva

numero_de_tramos = 14;
tramos_cto = zeros(numero_de_tramos,2);

tramos_cto(1,:) = [0 4530];
tramos_cto(2,:) = [165 575];
tramos_cto(3,:) = [0 1000];
tramos_cto(4,:) = [90 480];
tramos_cto(5,:) = [-75 600];
tramos_cto(6,:) = [-60 480];
tramos_cto(7,:) = [0 340];
tramos_cto(8,:) = [45 520];
tramos_cto(9,:) = [95 520];
tramos_cto(10,:) = [-125 520];
tramos_cto(11,:) = [0 800];
tramos_cto(12,:) = [170 500];
tramos_cto(13,:) = [1 0];
tramos_cto(14,:) = [0 0];
% Los ultimos dos tramos se dejan con la longitud/radio a 0 y con el angulo
% de curva a 1 o -1 para que el programa los calcule de manera optima para
% cerrar el circuito.