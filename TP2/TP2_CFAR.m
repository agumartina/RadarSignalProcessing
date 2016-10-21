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

%garantizamos que siemrpe sea par
if (mod(ref_win,2) ~= 0),
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
vecinity=2;                   % definimos cuantas celdas de vecindad
long_register=ref_win+2*vecinity+1; % cell + 2 of vecinity +1 cut

register=zeros(long_register);      
left_window=1:(ref_win/2);       %  
right_window=(long_register-ref_win/2+1):long_register;    % 
cut=ref_win/2+vecinity+1;                  % cell under test

C=2;            % parametro de forma

%%%%%% Bucle del CFAR %%%%%%

%llenado del llenado del registro
for l=1:(long_register-1)
    Pxx=abs(Z(l,n_pulso));                  % Intensidad
    register = circshift(register,1);       % Se corre todo un regitro ('clk')
    register(1)=Pxx;                        % se guarda
    
    %llenamos los primeros  
end

% calculamos el alfa para la pfa seleccionada
raizMpfa=nthroot(pfa,(ref_win));
alfa=((1-raizMpfa)/(raizMpfa/(ref_win-4)))^(1/2);

for l=long_register:L
    % Se obtiene la intensidad, se corre un lugar el resitro, y se ingresa
    % el nuevo valor
    Pxx=abs(Z(l,n_pulso));                  % Intensidad
    register = circshift(register,1);       % Se corre todo un regitro ('clk')
    register(1)=Pxx;                        % se guarda
    
    % parámetro de escala
    B=((1/(ref_win-4)).*(sum(register(left_window).^C)+sum(register(right_window).^C)))^(1/C); % 6

    % el T correspondiente a la posición de CUT 
    T(l-(ref_win/2+vecinity+1),n_pulso)=alfa.*B;	% 7
    %T(l,m)=((pfa^(-1/M)-1)*(sum(register(left_window).^C)+sum(register(cut).^C)+sum(register(right_window).^C)))^(1/C);      % 18
end

%--------------------------------------------------------------------------
% 3) Graficar la intensidad de los datos correspondientes al pulso 180 de la polarización VV.
% Utilizar la rutina anterior para graficar en la misma figura el umbral para M igual a 16 y
% 32, y para P F A igual a 10 −2 y 10 −3 . Repetir para el pulso 155 de VV. Comentar respecto
% del comportamiento del umbral en función de M y P F A .
%--------------------------------------------------------------------------

figure(3);
plot(X,abs(Z(:,n_pulso)),X, T(:,n_pulso));
title(strcat('Comparación Pulso ', num2str(n_pulso) , ' Intensidad VV - M= ',num2str(ref_win),' y Pfa= ', num2str(pfa)));
xlabel('Rabngo [m]'),ylabel('Intensidad');
legend('Intensidad Rx','Umbral T');

%--------------------------------------------------------------------------
% 4) Generar un mapa de detección a partir de los datos VV usando nuevamente la rutina
% anterior para procesar todos los pulsos. Asignar valor uno cuando se detecta un objetivo
% (el dato es mayor que el umbral) y cero cuando no hay detección. Calcular la probabilidad
% de falsa alarma empı́rica y comparar con la Pfa utilizada para establecer el umbral de
% detección. Comentar brevemente.
%--------------------------------------------------------------------------

%Definimos e inicializamos
detected=zeros(L,M);
pfa_empirica=zeros(M);
register=zeros(long_register); 

% calculamos el alfa para la pfa dada
raizMpfa=nthroot(pfa,(ref_win));
alfa=((1-raizMpfa)/(raizMpfa/(ref_win-4)))^(1/2);

%bucle del CFAR
for m=1:M
    for l=1:L
        % Se obtiene la intensidad, se corre un lugar el registro, y se ingresa
        % el nuevo valor
        Pxx=abs(Z(l,m));       % Intensidad
        register = circshift(register,1);       % Se corre todo un reistro ('clk')
        register(1)=Pxx;                        % se guarda
        
        % una vez que se llene el registro
        if(not((m==1)&&(l<long_register)))
            % parámetro de escala
            B=((1/(ref_win-4)).*(sum(register(left_window).^C)+sum(register(right_window).^C)))^(1/C); % 6

            % Posicionamos los punteros al CUT en la matriz T 
            posCut=l-(ref_win/2+vecinity+1);
            pulso=m;
            if (posCut<=0)                      %si fue negativo, es porque es del pulso anterior, si fue 0 es el último
                posCut=L+posCut;
                pulso=m-1;          
            end
            T(posCut,pulso)=alfa.*B;	% 7

            %Detector
            if T(posCut,pulso) < register(cut)
                detected(posCut,pulso)=1;
            elseif T(posCut,pulso) > register(cut)
                detected(posCut,pulso)=0;
            end
        end
    end
    % pfa_empirica, suma todos los match en cada ray
    % si detecta un match en el rango 2700 m (l=23 o 24), descuenta uno
    pfa_empirica(m)=sum(detected(:,m));
    if((detected(23,m)==1)||(detected(24,m)==1))
        pfa_empirica(m)=pfa_empirica(m)-1;
    end
    pfa_empirica(m)=pfa_empirica(m)/L;
end
%calculamos de las M mediciones
pfa_empirica_total=sum(pfa_empirica(:))/M;

%se reordena la matriz matchs

figure(4);

imagesc(Y, fliplr(X) ,detected);
title('Detection VV'),
ylabel('Rango [m]'),
xlabel('Nº PRF'),
colorbar,
caxis([0,1]);

figure(2);
plot(X,abs(Z(:,n_pulso)),X, T(:,n_pulso));
title(strcat('Comparación Pulso ', num2str(n_pulso) ,' Intensidad VV - M= ',num2str(ref_win),' y Pfa= ', num2str(pfa)));
xlabel('Rabngo [m]'),ylabel('Intensidad');
legend('Intensidad Rx','Umbral T');

% probabilidad de falsa añarma empírica
% false alarm rate = false targets per PRT / nº of rangecells
disp(strcat('El Pfa calculado es: ',num2str(pfa_empirica_total), ' y el Pfa establecido: ', num2str(pfa),' con un error de ', num2str(pfa_empirica_total-pfa) ));

% 5) Construir una rutina para implementar el detector CFAR que se describe en la Sección 3.1
% de [1], donde el parámetro de forma es desconocido. Generar un mapa de detección a partir
% de los datos VV para M igual a 16 y 32, y para P F A igual a 10 −2 y 10 −3 (usar la Figura
% 2 de [1] para determinar el umbral). Comparar con los resultados del inciso anterior.

%definición de variables del bucle
detected2=zeros(L,M);
register=zeros(ref_win);      % reiniciamos el registro del CFAR


%bucle del CFAR C y B desconocidos
for m=1:M
    for l=1:L
        % Se obtiene la intensidad, se corre un lugar el resitro, y se ingresa
        % el nuevo valor
        Pxx=abs(Z(l,m));                        % Intensidad
        register = circshift(register,1);       % Se corre todo un reistro ('clk')
        register(1)=Pxx;                        % se guarda

        % Método iterativo para encontrar parámetro de forma C
        testC=0;
        error=99;
        invC=0;
        while(error>=0.01)
            
            invC=(sum((register(left_window).^testC).*log10(register(left_window)))+sum((register(cut).^testC).*log10(register(cut)))+sum((register(right_window).^testC).*log10(register(right_window))) /...
                (sum(register(left_window).^testC)+sum(register(cut).^testC)+sum(register(right_window).^testC)))-...
                ((1/(ref_win-4)).*(sum(log10(register(left_window)))+log10(register(cut))+sum(log10(register(right_window))))); %ec 29
            
            error=testC-1/invC;
            if (error>0)
                testC=testC+0.01;
            elseif (error <0)
                testC=testC-0.01;
            end
        end    
        C=testC;
        
        % parámetro de escala
        B=((1/(ref_win-4)).*(sum(register(left_window).^C)+sum(register(cut).^C)+sum(register(right_window).^C)))^(1/C); % 6
        
        % alfa
        
        
        % threshold of the form
        T(l,m)=B*alfa^(1/c); % 7
        %T(l,m)=((pfa^(-1/M)-1)*(sum(register(left_window).^C)+sum(register(cut).^C)+sum(register(right_window).^C)))^(1/C);      % 18

        %Detector
        if T(l,m) < register(cut)
            detected2(l,m)=1;
        elseif T(l,m) > register(cut)
            detected2(l,m)=0;
        end
        
        % pfa_empirica, suma todos los match en cada ray
        % si detecta un match en el rango 2700 m, descuenta uno
        pfa_empirica(m)=sum(detected2(:,m));
        if((detected2(23,m)==1)||(detected2(24,m)==1))
            pfa_empirica(m)=pfa_empirica(m)-1;
        end
        pfa_empirica(m)=pfa_empirica(m)/L;
    end
end

