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
lambda=c/f0;


Lp = 2^(nextpow2(L)+1);
fsp=round(Lp/PRI);
Tp=1/fsp;
t_t  = (0:L-1)*T; 

fp_centred_PRF =linspace(-PRF/2,PRF/2,Lp); %Frequency Vector
vd=linspace(-PRF/2,PRF/2,4*M)*lambda/2;
range=(Tp*(0:Lp-1))*c/2;


% ------------------------------------------------------------------------
% Data files
% ------------------------------------------------------------------------
FileName = 'MartinWalking_Fs=12.5MHz_B=2.5MHz_Fc=900MHz_I16.bin' ;
FileId = fopen(FileName);
A = fread(FileId,[1 2],'uint32','b');
fseek(FileId,0,'bof') ;

%--------------------------------------------------------------------------
% Definición de variables para procesamiento
%--------------------------------------------------------------------------
% siendo L el largo de la sequencia de datos a transformar, usar
Signal = zeros(L,M,81); %reservamos memoria para la matriz de datos
SFiltrada = zeros(Lp,M,81); 

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
%       subplot(2,1,1),fftshift
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
    case 2% Velocidad Doppler
        for i=1:Lp
            fd(i,:)=abs(fft(SFiltrada(i,:,nk),L));
        end %Flat top 
        w=1-1.93*cos(2*pi*x/(M-1))+1.29*cos(4*pi*x/(M-1))-0.388*cos(6*pi*x/(M-1))+0.028*cos(8*pi*x/(M-1));
    case 3 %Hann
        w=0.5*(1-cos(2*pi*x/(M-1)));
    case 4 %Hamming
        w=25/46-(21/46)*cos(2*pi*x/(M-1));
    otherwise
        w=1;
end

% ------------------------------------------------------------------------
% Movie
% ------------------------------------------------------------------------
    
    set(gca,'nextplot','replace','Visible','off');
    mov=VideoWriter('Movie.mp4', 'Uncompressed AVI');
    mov.FrameRate=1/(PRI*M);
    open(mov);
    
%--------------------------------------------------------------------------
% Procesamiento
%--------------------------------------------------------------------------
% IMPLEMENTACION DEL FILTRO
% La siguiente rutina recorre todo el cubo de datos y procesa la señal
% recibida con el filtro Chirp y la almacena en SFiltrada

% Inicializamos fariables
fd = zeros(Lp, 4*M);

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
        
        % DOPPLER
        % Realixamos la FFT a M pulsos, en el Li que encontramos un máximo

        % buscamos un maximo por iteración para un dwell cualquiera
        [ValMax, Lmax]=max(SFiltrada(:,10,nk));

        % aplicamos el window
        win=SFiltrada(Lmax,:,nk).*w;
        
        % Convertimos las Lp muestras de 128 DWEL
        for i=1:Lp
          fd(i,:)=abs(fftshift(fft(SFiltrada(i,:,nk),4*M)));
        end    
        
        % se hace la fft y se centra
        tmpf=fftshift(fft(win,Lp));
               
        % normalizamos
        normA = tmpf - min(tmpf(:));
        normA = normA ./ max(normA(:));
         
        % graficamos
        %figure(33),
        h10 = figure(10);
        subplot(2,1,1),
        plot(fp_centred_PRF,abs(normA)),
        
        xlim([-150,150]),
        xlabel('Frecuencia'),
        ylabel('|Amplitud Normalizada|'),
        drawnow
                
        % RANGO / VELOCIDAD
        % Hacemos FFT a cada matriz M x DWEL y la ploteamos en un gráfico Velocidad
        % radial / Rango
                        
        % 2D map using view
        %figure(43);
        subplot(2,1,2),
        surf(vd, range/1000, fd,'EdgeColor', 'None', 'facecolor', 'interp'),
        xlabel('Velocidad Radial (m/s)'),
        ylabel('Range (Km)'),
        zlim([0,inf]),
        xlim([-6,6]),
        ylim([48,52]),
        view(2);
        colorbar;
        drawnow
        
        
        
        FramK = getframe(h10);
        writeVideo(mov,FramK);
end
close(mov);