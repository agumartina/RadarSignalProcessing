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
f0= 9.39e9;          % frec del radar
PRF=200;            % PRF 200 hz
blindRange=2000;    % Zona ciega, rango ciego
L=54;               % fast time meaurements
M=2048;             % slow time 
c=3e8;              % velocidad de la onda en ms
lambda=c/f0;

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

register=zeros(1,long_register);      
left_window=1:(ref_win/2);       %  
right_window=(long_register-ref_win/2+1):long_register;    % 
cut=ref_win/2+vecinity+1;                  % cell under test

C=2;            % parametro de forma

%%%%%% Bucle del CFAR %%%%%%

%llenado del llenado del registro
for l=1:(long_register-1)
    Pxx=abs(Z(l,n_pulso));                  % Intensidad
    register = circshift(register,1,2);       % Se corre todo un regitro ('clk')
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
    register = circshift(register,1,2);       % Se corre todo un regitro ('clk')
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
pfa_empirica=zeros(1,M);
%register=zeros(long_register); 

% calculamos el alfa para la pfa dada
raizMpfa=nthroot(pfa,(ref_win));
alfa=((1-raizMpfa)/(raizMpfa/(ref_win-4)))^(1/2);



%bucle del CFAR
for m=1:M
    for l=1:L
        % Se obtiene la intensidad, se corre un lugar el registro, y se ingresa
        % el nuevo valor
        Pxx=abs(Z(l,m));       % Intensidad
       
        register = circshift(register,1,2);       % Se corre todo un reistro ('clk')
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
pfa_empirica=zeros(1,M);

figure(4);
imagesc(Y, X ,flip(detected));
title('Detection VV'),
ylabel('Rango [m]'),
xlabel('Nº PRF'),
colorbar,
caxis([0,1]);

figure(5),
plot(X,abs(Z(:,n_pulso)),X, T(:,n_pulso)),
title(strcat('Comparación Pulso ', num2str(n_pulso) ,' Intensidad VV - M= ',num2str(ref_win),' y Pfa= ', num2str(pfa))),
xlabel('Rabngo [m]'),ylabel('Intensidad'),
legend('Intensidad Rx','Umbral T');

% probabilidad de falsa añarma empírica
% false alarm rate = false targets per PRT / nº of rangecells
disp(strcat('El Pfa calculado es: ',num2str(pfa_empirica_total), ' y el Pfa establecido: ', num2str(pfa),' con un error de ', num2str(pfa_empirica_total-pfa) ));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARTE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 5) Construir una rutina para implementar el detector CFAR que se describe en la Sección 3.1
% % de [1], donde el parámetro de forma es desconocido. Generar un mapa de detección a partir
% % de los datos VV para M igual a 16 y 32, y para P F A igual a 10 −2 y 10 −3 (usar la Figura
% % 2 de [1] para determinar el umbral). Comparar con los resultados del inciso anterior.
% 
% %  B (se toman M/2 valores a cada lado de la CUT,
% ifftshift
% %definición de variables del bucle
% detected2=zeros(L,M);
% pfa_empirica2=zeros(1,M);
% register2=zeros(1,long_register);  
% T=zeros(L,M);                       %reseteamos la matriz
% 
% %calculo de alfa
% done=1;
% while(done==1)
%     switch pfa
%         case 0.01
%             alfa=-log10(pfa);
%             disp('0.01');
%             done=0;
%         case 0.001 %de la figura 2 
%             disp('0.001');
%             if(ref_win==16)
%                 alfa=17;
%             elseif(ref_win==32)
%                 alfa=11;
%             end
%             done=0;
%         otherwise
%             disp('Debe ingressar como Pfa 10^-2 o 10^-3');
%             pfa=input('Ingrese la Pfa deseada: ');
%     end
% end
% 
% %bucle del CFAR C y B desconocidos
% for m=1:M
%     for l=1:L
%         % Se obtiene la intensidad, se corre un lugar el registro, y se ingresa
%         % el nuevo valor
%         Pxx=abs(Z(l,m));       % Intensidad
% 
%         register2 = circshift(register2,1,2);       % Se corre todo un reistro ('clk')
%         register2(1)=Pxx;                        % se guarda
%         
%         % una vez que se llene el registro
%         if(not((m==1)&&(l<=long_register)))
%             % parámetro de escala y forma
%             Xmerged=[register2(left_window),register2(right_window)]; 
%             parmhat = wblfit(Xmerged);
%             B=parmhat(1);
%             c=parmhat(2);
%             
%             % Posicionamos los punteros al CUT en la matriz T 
%             posCut=l-(ref_win/2+vecinity+1);
%             pulso=m;
%             if (posCut<=0)                      %si fue negativo, es porque es del pulso anterior, si fue 0 es el último
%                 posCut=L+posCut;
%                 pulso=m-1;          
%             end
%             % threshold of the form
%             T(posCut,pulso)=B*alfa^(1/c); % 32
% 
%             %Decision 
%             if T(posCut,pulso) < register2(cut)
%                 detected2(posCut,pulso)=1;
%             elseif T(posCut,pulso) > register2(cut)
%                 detected2(posCut,pulso)=0;
%             end
%         end
%     end
%     % pfa_empirica, suma todos los match en cada ray
%     % si detecta un match en el rango 2700 m (l=23 o 24), descuenta uno
%     pfa_empirica2(m)=sum(detected2(:,m));
%     if((detected2(23,m)==1)||(detected2(24,m)==1))
%         pfa_empirica2(m)=pfa_empirica2(m)-1;
%     end
%     pfa_empirica2(m)=pfa_empirica2(m)/L;
% end
% 
% %calculamos de las M mediciones
% pfa_empirica_total2=sum(pfa_empirica2(:))/M;
% 
% %se reordena la matriz matchs
% 
% figure(6);
% 
% imagesc(Y, fliplr(X) ,detected2);
% title('Detection VV ML CFAR Webull'),
% ylabel('Rango [m]'),
% xlabel('Nº PRF'),
% colorbar,
% caxis([0,1]);
% 
% figure(7);
% plot(X,abs(Z(:,n_pulso)),X, T(:,n_pulso));
% title(strcat('Comparación Pulso ', num2str(n_pulso) ,' Intensidad VV - M= ',num2str(ref_win),' y Pfa= ', num2str(pfa)));
% xlabel('Rabngo [m]'),ylabel('Intensidad');
% legend('Intensidad Rx','Umbral T');
% 
% % probabilidad de falsa añarma empírica
% % false alarm rate = false targets per PRT / nº of rangecells
% disp(strcat('El Pfa calculado es: ',num2str(pfa_empirica_total2), ' y el Pfa establecido: ', num2str(pfa),' con un error de ', num2str(pfa_empirica_total2-pfa) ));

% 6 Generar y graficar una imagen de rango/Doppler usando los primeros 256 pulsos de la
% matriz VV. Utilizar los comandos fft y fftshift. Repetir para las siguientes dos ventanas de
% 256 pulsos y observar el cambio de Doppler para el objetivo y el clutter.

% declaración e inicialización de variables
f=linspace(-PRF/2,PRF/2,256); %Frequency Vector
vd=linspace(-PRF/2,PRF/2,256)*lambda/2; % vector velocidad

% Inicializamos fariables
fd = zeros(8,L,256);

for i=1:8  
    %limites del DWELL
    iup=i*256;
    idown=iup-256+1;
    
    % Convertimos las Lp muestras de 128 DWELL
    for j=1:L     
        fd(i,j,:)=abs(fftshift(fft(ifftshift(abs(Z(j,idown:iup))))));
    end
    
    % se hace la fft y se centra
    tmpf=fftshift(fft(ifftshift(abs(Z(25,idown:iup)))));
       
    % normalizamos
    normA = tmpf - min(tmpf);
    normA = normA ./ max(normA(:));
    figure(9),
    plot(f/PRF,abs(normA)),
    title('Espectro de todo rango')
    xlabel('f/PRF'),
    ylabel('|Amplitud Normalizada|'),
    drawnow
     
    % 2D map using view
    ftmp=squeeze(fd(i,:,:));
    figure(8),
    imagesc(vd, X,ftmp);
    title('Grafico Velocidad/Rango Sin Filtro'),
    ylabel('Rango [m]'),
    xlabel('Vd'),
    colorbar,
    caxis([0,10]);
    drawnow
    
    pause(0.5);
end

% 7. Utilizar la rutina del detector de la Sección 2.1 de [1] para determinar la presencia de
% objetivos en cada banda de frecuencia Doppler. Graficar los resultados para M = 32 y
% Pfa = 10^−3 .

%Definimos e inicializamos
detected3=zeros(L,256);
pfa_empirica3=zeros(1,256);
register3=zeros(1,long_register);

% calculamos el alfa para la pfa dada
raizMpfa=nthroot(pfa,(ref_win));
alfa=((1-raizMpfa)/(raizMpfa/(ref_win-4)))^(1/2);

%bucle del CFAR
for m=1:8

    %limites del DWELL
    iup=m*256;
    idown=iup-256+1;
    
    
    for k=1:256 
        for l=1:L
            % Se obtiene la intensidad, se corre un lugar el registro, y se ingresa
            % el nuevo valor
            Pxx=abs(fd(m,l,k));       % Intensidad

            register3 = circshift(register3,1,2);       % Se corre todo un reistro ('clk')
            register3(1)=Pxx;                           % se guarda

            % una vez que se llene el registro
            if(not((k==1)&&(l<long_register)))
                % parámetro de escala
                B=((1/(ref_win-4)).*(sum(register3(left_window).^C)+sum(register3(right_window).^C)))^(1/C); % 6

                % Posicionamos los punteros al CUT en la matriz T 
                posCut=l-(ref_win/2+vecinity+1);
                pulso=idown+k;
                if (posCut<=0)                      %si fue negativo, es porque es del pulso anterior, si fue 0 es el último
                    posCut=L+posCut;
                    pulso=idown+k-1;          
                end
                T(posCut,pulso)=alfa.*B;	% 7

                %Detector
                if T(posCut,pulso) < register3(cut)
                    detected3(posCut,pulso)=1;
                elseif T(posCut,pulso) > register3(cut)
                    detected3(posCut,pulso)=0;
                end
            end
        end
    end
    
    % pfa_empirica, suma todos los match en cada ray
    % si detecta un match en el rango 2700 m (l=23 o 24), descuenta uno
    pfa_empirica3(k)=sum(detected3(:,k));
    if((detected3(23,k)==1)||(detected3(24,k)==1))
        pfa_empirica3(k)=pfa_empirica3(k)-1;
    end
    pfa_empirica3(k)=pfa_empirica3(k)/L;
end

%calculamos de las M mediciones
pfa_empirica_total3=sum(pfa_empirica3(:))/M;

figure(10);
imagesc(vd/PRF, X,detected3);
title('Detection VV Doppler'),
ylabel('Rango [m]'),
xlabel('Vd/PRF'),
colorbar,
caxis([0,1]);

figure(11),
plot(X,abs(squeeze(fd(1,:,n_pulso))),X, T(:,n_pulso)),
title(strcat('Comparación T Doppler Pulso ', num2str(n_pulso) ,' Intensidad VV - M= ',num2str(ref_win),' y Pfa= ', num2str(pfa))),
xlabel('Rabngo [m]'),ylabel('Intensidad'),
legend('Intensidad Rx','Umbral T');

% probabilidad de falsa añarma empírica
% false alarm rate = false targets per PRT / nº of rangecells
disp(strcat('El Pfa para Doppler calculado es: ',num2str(pfa_empirica_total3), ' y el Pfa establecido: ', num2str(pfa),' con un error de ', num2str(pfa_empirica_total3-pfa) ));