%--------------------------------------------------------------------------
% TRABAJO PRÁCTICO Nº2
% PROCESAMIENTO DE SEÑALES DE RADAR
% MAESTRÍA EN RADARES E INSTRUMENTACIÓN UNC-IUA
%--------------------------------------------------------------------------
clc;
clear all;
%--------------------------------------------------------------------------
% Declaración de variables
%--------------------------------------------------------------------------
fs= 10e6;           % frec de sampleo
PRF=200;            % PRF 200 hz
blindRange=2000;    % Zona ciega, rango ciego
L=54;               % fast time meaurements
M=2048;             % slow time 
c=3e8;              % velocidad de la onda en ms


X=2000+c*(1:L)/fs;
Y=1:M;


%--------------------------------------------------------------------------
% 1) Apertura del archivo y graficar VV
%--------------------------------------------------------------------------
data=load('procNov11stare0.mat');
Z=rot90(data.vv,3);

figure(1);
phIm=pcolor(Y, X ,abs(Z));
set(phIm,'edgecolor','none');
title('|(procNov11stare0.mat)| VV');
ylabel('Rango [m]');
xlabel('Nº PRF');
colorbar;

%--------------------------------------------------------------------------
% 2) calcular T para un pulso
%--------------------------------------------------------------------------

%se le solicita al usuario que ingrese una Pfa deseada
pfa=-1;
while (pfa<0)
    pfa=input('Ingrese la Pfa deseada (10^-6): ');
end;

%se le solicita al usuario que ingrese el tamaño de la ventana de
%referencia
ref_win=-1;
while (ref_win<0)
    ref_win=input('Ingrese la cantidad de muestras que debe poseer la ventana de referencia:');
end

if (mod(ref_win,2) == 0),
    ref_win= ref_win+1;
end;

% definir que pulso procesar
n_pulso=-1;
while(n_pulso<1)
    n_pulso=input('Ingrese el pulso que desea procesar:');
end

%Definimos e inicializamos
T=zeros(L,M);

%actualizamos los punteros a la celdas de referencia
register=zeros(ref_win);      % cell + 2 of vecinity
left_window=1:(ref_win/2-2.5);       %  
right_window=(ref_win/2+3.5):ref_win;    % 
cut=ref_win/2+0.5;                  % cell under test

C=2;            % parametro de forma

%bucle del CFAR

for l=1:L
    % Se obtiene la intensidad, se corre un lugar el resitro, y se ingresa
    % el nuevo valor
    Pxx=Z(l,n_pulso).*conj(Z(l,n_pulso))/(M*L);       % Intensidad
    register = circshift(register,1);       % Se corre todo un reistro ('clk')
    register(1)=Pxx;                        % se guarda

    %alfa
    raizMpfa=nthroot(pfa,(ref_win-4));
    alfa=((1-raizMpfa)/(raizMpfa/(ref_win-4)))^(1/2);
    % parámetro de escala
    B=((1/(ref_win-4)).*(sum(register(left_window).^C)+sum(register(cut).^C)+sum(register(right_window).^C)))^(1/C); % 6

    % threshold of the form
    T(l,n_pulso)=alfa.*B; % 7
    %T(l,m)=((pfa^(-1/M)-1)*(sum(register(left_window).^C)+sum(register(cut).^C)+sum(register(right_window).^C)))^(1/C);      % 18
end

%--------------------------------------------------------------------------
% 3) Graficar la intensidad de los datos correspondientes al pulso 180 de la polarización VV.
% Utilizar la rutina anterior para graficar en la misma figura el umbral para M igual a 16 y
% 32, y para P F A igual a 10 −2 y 10 −3 . Repetir para el pulso 155 de VV. Comentar respecto
% del comportamiento del umbral en función de M y P F A .
%--------------------------------------------------------------------------

%en Z: la señal recibida, T el umbral
if(ref_win==17)
    ref_win=16;
elseif(ref_win==33)
    ref_win=32;
end
figure(3);
plot(X,Z(:,n_pulso).*conj(Z(:,n_pulso))/(M*L),X, T(:,n_pulso));
title(strcat('Comparación Pulso 180 Intensidad VV y T - M= ',num2str(ref_win),' y Pfa= ', num2str(pfa)));
xlabel('Rabngo [m]'),ylabel('Intensidad');

%--------------------------------------------------------------------------
% 4) Generar un mapa de detección a partir de los datos VV usando nuevamente la rutina
% anterior para procesar todos los pulsos. Asignar valor uno cuando se detecta un objetivo
% (el dato es mayor que el umbral) y cero cuando no hay detección. Calcular la probabilidad
% de falsa alarma empı́rica y comparar con la P F A utilizada para establecer el umbral de
% detección. Comentar brevemente.
%--------------------------------------------------------------------------

%Definimos e inicializamos
detected=zeros(L,M);


%bucle del CFAR
for m=1:M
    for l=1:L
        % Se obtiene la intensidad, se corre un lugar el resitro, y se ingresa
        % el nuevo valor
        Pxx=Z(l,m).*conj(Z(l,m))/(M*L);       % Intensidad
        register = circshift(register,1);       % Se corre todo un reistro ('clk')
        register(1)=Pxx;                        % se guarda

        raizMpfa=nthroot(pfa,(ref_win-4));
        alfa=((1-raizMpfa)/(raizMpfa/(ref_win-4)))^(1/2);
        % parámetro de escala
        B=((1/(ref_win-4)).*(sum(register(left_window).^C)+sum(register(cut).^C)+sum(register(right_window).^C)))^(1/C); % 6

        % threshold of the form
        T(l,m)=alfa.*B; % 7
        %T(l,m)=((pfa^(-1/M)-1)*(sum(register(left_window).^C)+sum(register(cut).^C)+sum(register(right_window).^C)))^(1/C);      % 18

        %Detector
        if T(l,m) < register(cut)
            detected(l,m)=1;
        elseif T(l,m) > register(cut)
            detected(l,m)=0;
        end
    end
end
figure(4);
title('Detection VV');
phIm=image(Y, X ,detected);
set(phIm,'edgecolor','none');
caxis([0,1]);
ylabel('Rango [m]');
xlabel('Nº PRF');
colorbar;


