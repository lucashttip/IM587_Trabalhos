clear all
clc

% Constantes geométricas são declaradas no código mecanismo
load functions.mat

%% Declaração dos parâmetros

global functions = {Ff, iJf, Kf, Lf, P1f, KP1f, LP1f, P2f, KP2f, LP2f, P3f, KP3f, LP3f}

global g M1 M2 M3 I1 I2 I3 Fq Faux Pe functions


% Análise 1: Sem força de atuação
Fq = 0;               % Newton
q0 = deg2rad(-6.5);
dq0 = deg2rad(60);



% Análise 2: Com força de atuação
#Fq = 1;
#q0 = deg2rad(-6.5);
#dq0 = deg2rad(0);



g = 9.81; 
% Barra 1
M1 = 0.175;        % kg
%I1 = 1.316*10^6   % g/mm^2
I1 = 1.316*10^-3;   % kg/m^2



Faux = Fq;

% Barra 2
M2 = 0.263;        % kg
%I2 = 4.44*10^6    % g/mm^2
I2 = 4.44*10^-3;    % kg/m^2


% Barra 3
M3 = 0.697;        % kg
%I3 = 5.259*10^6   % g/mm^2
I3 = 5.259*10^-3; % kg/m^2

% Condições para simulação
ti = 0;
tf = 1;
temp = [ti tf];



Y0 = [q0
      dq0];

%% Resolução da equação de movimento
options = odeset('RelTol', 1e-8,'AbsTol',1e-10);   
[t,y] = ode45(@MEC_4Barras, temp, Y0, options);

A = zeros(size(y(:,1)));
B = zeros(size(y(:,2)));

for i =1:size(A,1)
  [A(i), B(i)] = newtonR2(functions{1},functions{2},y(i,1),pi,0,1e-5,15);
end


plota(t,rad2deg(y(:,1)),'Tempo [s]', 'Posição de q [graus]', strcat("imagens/Posq",num2str(Faux),".png"))

plota(t,y(:,2),'Tempo [s]', 'Velocidade de q [rad/s]', strcat("imagens/Velq",num2str(Faux),".png"))

plota(t,rad2deg(A),'Tempo [s]', 'Posição de A [graus]', strcat("imagens/PosA",num2str(Faux),".png"))

plota(t,rad2deg(B),'Tempo [s]', 'Posição de B [graus]', strcat("imagens/PosB",num2str(Faux),".png"))

%{
%}
