function dy = MEC_4Barras(t,y)
  # Obter as variáveis aqui em cima
global g M1 M2 M3 I1 I2 I3 Fq Faux Pe functions

Ff = functions{1};
iJf = functions{2};
Kf = functions{3};
Lf = functions{4};
P1f = functions{5};
KP1f = functions{6};
LP1f = functions{7};
P2f = functions{8};
KP2f = functions{9};
LP2f = functions{10};
P3f = functions{11};
KP3f = functions{12};
LP3f = functions{13};

% Eq. de Eksergian é: j*qdd + 0.5*djdq*qd = Q

% y=[q
%    dq];

q=y(1);
dq=y(2);

%Cinematica das Secundárias A e X
[A, B] = newtonR2(Ff,iJf,q,pi,0,1e-5,15);

%if B <=0 || B >=deg2rad(90)
%   dq = -dq;
%   Fq = -Fq;
%   q = q+dq;
%   [A, B] = newtonR2(Ff,iJf,q,pi,0,1e-5,15);
%end


if B >=deg2rad(25)
   Fq = 0;
end

if Fq == 0 && Faux > 0 && B < deg2rad(20)
   Fq = Faux;
end



K=Kf(A,B,q);
L = Lf(A,B,q);

%Cinematica do Ponto Int. Barra 1
P1 = P1f(q) .* 10^-2;

Kp1=KP1f(q).* 10^-2;

Lp1 = LP1f(q) .* 10^-2;

%Cinematica do Ponto Int. Barra 2
P2 = P2f(A,q) .* 10^-2;

Kp2=KP2f(A,B,q) .* 10^-2;

Lp2 = LP2f(A,B,q) .* 10^-2;

%Cinematica do Ponto Int. Barra 3
P3 = P3f(B) .* 10^-2;

Kp3=KP3f(A,B,q) .* 10^-2;

Lp3 = LP3f(A,B,q) .* 10^-2;

%Termo Energia Potencial dV/dq
% V = M1*g*X1 + M2*g*X2 + M3*g*X3
% dV/dq = M1*g*dX1/dq + M2*g*dX2/dq +M3*g*dX3/dq
%=> X3= X + delta => dX3/dq=dX/dq=Kx
dVdq = M1*g*Kp1(2) + M2*g*Kp2(2) +M3*g*Kp3(2);

%Força Generalizada
Qnc=  Fq*22.5*cos(q);

%display(Qnc)

%Inercia Generalizada
Iq= M1*(Kp1(1)^2 + Kp1(2)^2) + I1 + M2*(Kp2(1)^2 + Kp2(2)^2) + I2*(K(1)^2) + M3*(Kp3(1)^2 + Kp3(2)^2) + I3*(K(2)^2);

%Derivada Inercia Generalizada dIq = 0.5 * dIq/dq
dIq=M1*(Kp1(1)*Lp1(1) + Kp1(2)*Lp1(2)) + M2*(Kp2(1)*Lp2(1) + Kp2(2)*Lp2(2)) + I2*(K(1)*L(1)) + M3*(Kp3(1)*Lp3(1) + Kp3(2)*Lp3(2)) + I3*(K(2)*L(2));

d2q= (Qnc -dVdq -dIq*(dq^2))/Iq;

dy=[dq
    d2q];
  
end