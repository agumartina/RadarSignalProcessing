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
reZ=abs(real(Z));
imZ=abs(imag(Z));

%graph
figure(1);
subplot(1,2,1);
title('|Re(procNov11stare0.mat)|');
ph=pcolor(Y, X ,reZ);
set(ph,'edgecolor','none');
xlabel('Rango [m]');
ylabel('Nº PRF');
colorbar;

subplot(1,2,2);
title('|Im(procNov11stare0.mat)|');
ph=pcolor(Y, X ,imZ);
set(ph,'edgecolor','none');
xlabel('Rango [m]');
ylabel('Nº PRF');
colorbar;



