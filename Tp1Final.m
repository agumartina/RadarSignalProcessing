% ------------------------------------------------------------------------
% Notes:
% The script shows how to read the data file and organize the data in radar
% dwells of M pulses.
% ------------------------------------------------------------------------
clc;
clear all;
% ------------------------------------------------------------------------
% Signal parameters
% ------------------------------------------------------------------------
c = 3e8;        % Signal propagation velocity [m/s]
f0 = 0.9e9;     % Carrier freq of the transmitted signal [Hz]
fs = 12.5e6;    % DAC rate samples [samples/s] 
M = 128;        % Pulses per dwell
TauP = 0.001;   % Pulse length [s]
PRI = 1*TauP;   % Slow time sampling interval or PRI [s]
PRF=1/PRI;      % Pulse repetition frequency
beta = 2.5e6 ;  % Chirp signal band-width [Hz]
L = round(fs*PRI); % Number of range cells
T=1/fs;         % período de sampleo hz

Lp = 2^(nextpow2(L)+1);
fsp=round(Lp/PRI);
Tp=1/fsp;
t_t  = (0:L-1)*T; 
t_tp = (0:Lp-1)*Tp;
% t_f = (-(Lp-1)/2:(Lp-1)/2)*Tp;       %tiempo para graficar en f
% 
% f= fs*(-beta/2+152:beta/2-152-1)/(beta);
% fp= linspace(0,beta,Lp)*fs/beta;
% fp_centred=linspace(-beta/2,beta/2,Lp);
fp_centred = fs*(-Lp/2:Lp/2-1)/Lp; %Frequency Vector
fp_centred_f0 =PRF*(-(Lp/2):(Lp/2-1))/Lp; %Frequency Vector
fp_centred_128 =-(64):(64-1); %Frequency Vector



% ------------------------------------------------------------------------
% Data files
% ------------------------------------------------------------------------
FileName = 'MartinWalking_Fs=12.5MHz_B=2.5MHz_Fc=900MHz_I16.bin' ;
FileId = fopen(FileName);
A = fread(FileId,[1 2],'uint32','b');
fseek(FileId,0,'bof') ;

% ------------------------------------------------------------------------
% Movie
% ------------------------------------------------------------------------
% h10 = figure(10); 
% FrameRate = 1/(PRI*M);
% set(gca,'nextplot','replace','Visible','off');
% mov=VideoWriter('Movie.mp4');

%--------------------------------------------------------------------------
% Definición de variables para procesamiento
%--------------------------------------------------------------------------
% siendo L el largo de la sequencia de datos a transformar, usar
Signal = zeros(L,M,81); %reservamos memoria para la matriz de datos
SFiltrada = zeros(Lp,M,81); 

% terting para 1 p/ primer iteración
nDWEL=100;
nk=10;

% ------------------------------------------------------------------------
% Chirp - Adaptive Filter
% ------------------------------------------------------------------------

% calculamos nuevamente 4.82
theta_t=(pi*beta*(t_t.^2))./TauP;
x_t=exp(1i*theta_t);

%de la ecuación 4.55 del libro, tenemos
hp=flip(x_t');
Hp_f=fft(hp,Lp);     % transformamos y corremo

% ------------------------------------------------------------------------
% Processing
% ------------------------------------------------------------------------
k=0;
yExtra = [] ;

while (~feof(FileId) ) 
    k = k + 1;
    
    % Read data: 
    %   - batch of 2M samples (complex numbers)
    %   - fill vectors of size LM
    yTemp = [] ;
    if ((L*M-size(yExtra,1)) > 0 )
        for i=1:(L*M-size(yExtra,1))/A(2)+1
            A = fread(FileId,[1 2],'uint32','b');
            if ( ~isempty(A) )
                recv = fread(FileId,fliplr(A),'int16','b');
                yTemp((i-1)*A(2)+1:(i)*A(2),1) = recv(:,1) + 1i*recv(:,2);  
            end
        end
    end
    
    yTemp = [yExtra ; yTemp] ;
    if ( length(yTemp) > L*M )
        y = yTemp(1:L*M,1) ; % Throw out the extra samples
        yExtra = yTemp(L*M+1:end,:);
    else
        y = yTemp ;
        yExtra =[] ;
    end
    
% Movie
    if ( length(y) == L*M )
        
        Y=reshape(y,L,M); % input, filas_output, columnas_output
        Signal(:,:,k)=Y;

% -----------------------------------------------------------------
%       figure(h10);
%       subplot(2,1,1),
%       pcolor((1:M),(1:L),real(Y)), shading flat, 
%       xlabel('Pulses per dwell'); 
%       ylabel('Number of range cells');
%       subplot(2,1,2),
%       pcolor((1:M)*PRI,(1:L)*c/2/fs/1000,real(Out)), shading flat, 
%       xlabel('Pulses per dwell'); 
%       ylabel('Number of range cells');
        
%       pcolor((1:M)*PRI, (1:L)*c/2/fs/1000, real(Out) ), shading flat;
%       FramK = getframe(h10);
%       mov = addframe(mov,FramK);
    end
end

%mov = close(mov);
fclose(FileId) ;


%--------------------------------------------------------------------------
% Funciones Window
%--------------------------------------------------------------------------

disp('Seleccione la función windos que desee implementar:');
disp(sprintf ('1) Rectangular\n2) Flat top\n3) Hann\n4) Hamming\nSin window otro número'));
cw=input('Ingrese Opcion: ');
x=1:M;
switch cw
    case 1 %Rectangular
        w=1-abs((x-0.5*(M-1))/(0.5*M));
    case 2 %Flat top 
        w=1-1.93*cos(2*pi*x/(M-1))+1.29*cos(4*pi*x/(M-1))-0.388*cos(6*pi*x/(M-1))+0.028*cos(8*pi*x/(M-1));
    case 3 %Hann
        w=0.5*(1-cos(2*pi*x/(M-1)));
    case 4 %Hamming
        w=25/46-(21/46)*cos(2*pi*x/(M-1));
    otherwise
        w=1;
end

%--------------------------------------------------------------------------
% Procesamiento
%--------------------------------------------------------------------------
% IMPLEMENTACION DEL FILTRO
% La siguiente rutina recorre todo el cubo de datos y procesa la señal
% recibida con el filtro Chirp y la almacena en SFiltrada
for nk=1:81
        for nDWEL=1:M
            %Tomamos el primer ray y Aplicamos FFT
            ray_t=Signal(:,nDWEL,nk);
            ray_f=fft(ray_t, Lp);

            %Aplicamos el Filtro en el dominio de F y lo volvemos a
            % transformar en t
            out_f=ray_f.*Hp_f;
            out=ifft(out_f, Lp);
            
            SFiltrada(:,nDWEL,nk)=out;
               
        end
end


% TEST DOPPLER SIN FILTRO %
plano = zeros(81,4*M);
promedio = zeros (81, M);
for nk=1:81
        %buscamos un maximo por iteración para un dwell cualquiera
        [ValMax, Lmax]=max(SFiltrada(:,10,nk));

        % aplicamos el window
      
        win=SFiltrada(Lmax,:,nk).*w;
            
        % se hace la fft y se centra
        tmpf=fftshift(fft(win,Lp));
        
        %tmpf.*conj(tmpf)/(L*L); %computing power with proper scaling
        
%grafico
        Z=abs(tmpf.*conj(tmpf)/(L*L));
        
        figure(33),
        plot(fp_centred_f0(-TauP/10:TauP/10)/TauP,Z),
        %ylim([0,1e9]),
        xlabel('Frequency'),
        ylabel('Amplitud'),
        shading interp;
        drawnow
        pause(1);
end

% PROCESAMIENTO DOPPLER
% La siguiente rutina recorre todo el cubo de datos y procesa la señal
% recibida con el filtro Chirp y la almacena en SFiltrada

%TEST
% % %malloc de las variables a utilizar
% plano = zeros(81,2*M);
% promedio = zeros (81, M);
% 
% % implementación del window
% %window function
% x=1:128;
% 
% w=sin(pi*x/(128-1)).^2;         %HANN
% w=1;
% for nk=1:81
%         %calculamos el promedio del pulso
%         tmp=0;
%         for Li=1:Lp   
%             %se lee una fila en slow time
%             tmp=tmp+SFiltrada(Li,:,nk);
%         end
%         
%         promedio(nk,:)=tmp/Lp;
%         
%         % aplicamos el window
%         win=promedio(nk,:)*w;
%             
%         % se hace la fft y se centra
%         tmpf=fftshift(fft(win,2*M));
%         
%         plano(nk,:)=tmpf.*conj(tmpf)/(Lp*Lp); %computing power with proper scaling
%         
%         Z=plano(nk,:);
%         
%         figure(33),
%         plot(fd,Z),
%         ylim([-10,50]),
%         xlabel('Frequency'),
%         ylabel('Range'),
%         shading interp;
%         drawnow
% end

%--------------------------------------------------------------------------
% Grafuicamos
%--------------------------------------------------------------------------
% 
% %graph f salida
% out_f_c=fftshift(out_f);
% %out_DB=20*log10(out_f_c/max(out_f_c));
% %Pxx=out_f_c.*conj(out_f_c)/(Lp*Lp); %computing power with proper scaling
% figure(1),
% plot(fp_centred/beta,abs(out_f_c)/Lp),
% title('espectro salida');
% 
% %graph f rx
% ray_f_c=fftshift(ray_f);
% %ray_DB=20*log10(ray_f_c/max(ray_f_c));
% %PxxR=ray_f_c.*conj(ray_f_c)/(Lp*Lp); %computing power with proper scaling
% figure(5),
% plot(fp_centred/beta, abs(ray_f_c)/Lp),
% title('espectro rx');
% 
% %graph f chirp
% Hp_f_c=fftshift(Hp_f);
% %chi_DB=20*log10(Hp_f_c/max(Hp_f_c));
% %PxxH=Hp_f_c.*conj(Hp_f_c)/(Lp*Lp); %computing power with proper scaling
% figure(6),
% plot(fp_centred/beta, abs(Hp_f_c)/Lp),
% % ylim([-10,1]),
% % xlim([-0.6, 0.6]),
% title('espectro chirp');
% 
% %-------------- tiempo ----------------
% %graph tiempo
% figure(2),
% plot(t_tp/TauP, transpose(out)),
% title('Señal Filtrada'),
% ylabel('Amplitud'),
% xlabel('t/TauP');
% 
% %graph tiempo
% figure(3),
% plot(t_t/TauP, fliplr(x_t)),
% title('Señal Chirp'),
% ylabel('Amplitud'),
% xlabel('t/TauP');
% 
% %graph tiempo
% figure(4),
% plot(t_t/TauP, ray_t),
% title('RX'),
% ylabel('Amplitud'),
% xlabel('t/TauP');