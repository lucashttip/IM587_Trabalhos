clc
clear

### ====== Constantes =======
% Barra 1
l1v = 12.5;
U1v = -2.5;
V1v = 0;

% Barra 2
C2v = 45;
U2v = 22;
V2v = 0;

% Barra 3
C3v = 35;
U3v = -17.5;
V3v = 0;
l2v = 1.5;

% Distâncias dos apoios
L1v = 12.5*cosd(6.5) + 1.5;
L2v = 45 - l1v*sind(6.5);

### ====== Resolvendo o problema com simbolico ======

### ========================
### ========================
### ======= Entrada ========
### ========================
### ========================

# Pacote e variáveis
pkg load symbolic
syms l1 l2 L1 L2 C2 C3 A B q 
syms Aq(q) Bq(q) Ka Kaq(q) Kb Kbq(q) La Laq(q) Lb Lbq(q) qp qpp
syms U1 U2 U3 V1 V2 V3

# Equações de loop
f1 = l1*cos(q) + C2*sin(A) + l2*cos(B) - L1;
f2 = l1*sin(q) - C2*cos(A) - l2*sin(B) - L2;

F = [f1; f2];

# Equações dos pontos de interesse
X1 = U1*cos(q) - V1*sin(q);
Y1 = U1*sin(q) + V1*cos(q);
P1 = [X1; Y1];

X2 = l1*cos(q) + U2*sin(Aq) + V2*cos(Aq);
Y2 = l1*sin(q) - U2*cos(Aq) + V2*sin(Aq);
P2 = [subs(X2,Aq,A); subs(Y2,Aq,A)];

X3 = L1 + U3*cos(Bq) + V3*sin(Bq);
Y3 = L2 - U3*sin(Bq) + V3*cos(Bq);
P3 = [subs(X3,Bq,B); subs(Y3,Bq,B)];

### ================================
### ================================
### ======= Solução do Loop ========
### ================================
### ================================

# jacobiano
J = jacobian(F,[A B]);

# Inversa do jacobiano
iJ = simplify(inv(J));

# Vetor b
b = -jacobian(F,q);

# Coeficientes de velocidade das secundárias (K = [Ka, Kb])
K = simplify(iJ*b);

# Calculando os coeficientes de aceleração das secundárias
Kq = subs(K,{A,B},{Aq,Bq});
JL = jacobian(Kq,q);

# Coeficientes de aceleração das secundárias
L = subs(subs(JL,{diff(Aq,q),diff(Bq,q)},{Ka,Kb}),{Aq,Bq},{A,B});
L = simplify(L);


### =================================================
### =================================================
### ======= Soluções dos Pontos de Interesse ========
### =================================================
### =================================================

## Equações do ponto de interesse:
%xp = C3*cos(Bq);
%yp = C3*sin(Bq);
%P = [subs(xp,Bq,B); subs(yp,Bq,B)];

# coeficientes de velocidade do ponto de interese
# Ponto 1
KX1q = subs(diff(X1,q),{diff(Aq,q),diff(Bq,q)},{Kaq,Kbq});
KY1q = subs(diff(Y1,q),{diff(Aq,q),diff(Bq,q)},{Kaq,Kbq});

# Ponto 2
KX2q = subs(diff(X2,q),{diff(Aq,q),diff(Bq,q)},{Kaq,Kbq});
KY2q = subs(diff(Y2,q),{diff(Aq,q),diff(Bq,q)},{Kaq,Kbq});

# Ponto 3
KX3q = subs(diff(X3,q),{diff(Aq,q),diff(Bq,q)},{Kaq,Kbq});
KY3q = subs(diff(Y3,q),{diff(Aq,q),diff(Bq,q)},{Kaq,Kbq});



# Coeficiente de aceleração do ponto de interesse
# Ponto 1
LX1q = subs(subs(diff(KX1q,q),{diff(Kaq,q),diff(Kbq,q)},{Laq,Lbq}),{diff(Aq,q),diff(Bq,q)},{Kaq,Kbq});
LY1q = subs(subs(diff(KY1q,q),{diff(Kaq,q),diff(Kbq,q)},{Laq,Lbq}),{diff(Aq,q),diff(Bq,q)},{Kaq,Kbq});

# Ponto 2
LX2q = subs(subs(diff(KX2q,q),{diff(Kaq,q),diff(Kbq,q)},{Laq,Lbq}),{diff(Aq,q),diff(Bq,q)},{Kaq,Kbq});
LY2q = subs(subs(diff(KY2q,q),{diff(Kaq,q),diff(Kbq,q)},{Laq,Lbq}),{diff(Aq,q),diff(Bq,q)},{Kaq,Kbq});

# Ponto 3
LX3q = subs(subs(diff(KX3q,q),{diff(Kaq,q),diff(Kbq,q)},{Laq,Lbq}),{diff(Aq,q),diff(Bq,q)},{Kaq,Kbq});
LY3q = subs(subs(diff(KY3q,q),{diff(Kaq,q),diff(Kbq,q)},{Laq,Lbq}),{diff(Aq,q),diff(Bq,q)},{Kaq,Kbq});

# Simplificando
KX1 = simplify(subs(KX1q,{Aq,Bq,Kaq,Kbq},{A,B,Ka,Kb}));
KY1 = simplify(subs(KY1q,{Aq,Bq,Kaq,Kbq},{A,B,Ka,Kb}));

KX2 = simplify(subs(KX2q,{Aq,Bq,Kaq,Kbq},{A,B,Ka,Kb}));
KY2 = simplify(subs(KY2q,{Aq,Bq,Kaq,Kbq},{A,B,Ka,Kb}));

KX3 = simplify(subs(KX3q,{Aq,Bq,Kaq,Kbq},{A,B,Ka,Kb}));
KY3 = simplify(subs(KY3q,{Aq,Bq,Kaq,Kbq},{A,B,Ka,Kb}));

KP1 = [KX1; KY1];
KP2 = [KX2; KY2];
KP3 = [KX3; KY3];

LX1 = simplify(subs(LX1q,{Aq,Bq,Kaq,Kbq,Laq,Lbq},{A,B,Ka,Kb,La,Lb}));
LY1 = simplify(subs(LY1q,{Aq,Bq,Kaq,Kbq,Laq,Lbq},{A,B,Ka,Kb,La,Lb}));

LX2 = simplify(subs(LX2q,{Aq,Bq,Kaq,Kbq,Laq,Lbq},{A,B,Ka,Kb,La,Lb}));
LY2 = simplify(subs(LY2q,{Aq,Bq,Kaq,Kbq,Laq,Lbq},{A,B,Ka,Kb,La,Lb}));

LX3 = simplify(subs(LX3q,{Aq,Bq,Kaq,Kbq,Laq,Lbq},{A,B,Ka,Kb,La,Lb}));
LY3 = simplify(subs(LY3q,{Aq,Bq,Kaq,Kbq,Laq,Lbq},{A,B,Ka,Kb,La,Lb}));

LP1 = [LX1; LY1];
LP2 = [LX2; LY2];
LP3 = [LX3; LY3];


### ========================================
### ========================================
### ======= Substituindo os valores ========
### ========================================
### ========================================

## Transformando as expressões simbolicas em funções

# Expressões para calculo da posição:
iJf1 = function_handle(iJ,'vars', [A B l1 l2 L1 L2 C2 C3]);
iJf = @(A,B) iJf1(A, B, l1v, l2v, L1v, L2v, C2v, C3v);

Ff1 = function_handle(F,'vars', [A B q l1 l2 L1 L2 C2 C3]);
Ff = @(A,B,q) Ff1(A, B, q, l1v, l2v, L1v, L2v, C2v, C3v);

# Expressões para calculo da velocidade e aceleração das secundárias
Kf1 = function_handle(K,'vars', [A B q l1 l2 L1 L2 C2 C3]);
Kf = @(A, B, q) Kf1(A, B, q, l1v, l2v, L1v, L2v, C2v, C3v);

Lf1 = function_handle(L,'vars', [A B q Ka Kb l1 l2 L1 L2 C2 C3]);
Lf = @(A, B, q) Lf1(A, B, q, Kf(A,B,q)(1), Kf(A,B,q)(2), l1v, l2v, L1v, L2v, C2v, C3v);


# Expressões para calculo da posição, velocidade e aceleração do ponto de interesse
P1f1 = function_handle(P1,'vars', [q U1 V1 l1 l2 L1 L2 C2 C3]);
P2f1 = function_handle(P2,'vars', [A q U2 V2 l1 l2 L1 L2 C2 C3]);
P3f1 = function_handle(P3,'vars', [B U3 V3 l1 l2 L1 L2 C2 C3]);

P1f = @(q) P1f1(q, U1v, V1v,l1v, l2v, L1v, L2v, C2v, C3v);
P2f = @(A,q) P2f1(A,q, U2v, V2v,l1v, l2v, L1v, L2v, C2v, C3v);
P3f = @(B) P3f1(B, U3v, V3v,l1v, l2v, L1v, L2v, C2v, C3v);

KP1f1 = function_handle(KP1,'vars', [q U1 V1 l1 l2 L1 L2 C2 C3]);
KP2f1 = function_handle(KP2,'vars', [A q Ka U2 V2 l1 l2 L1 L2 C2 C3]);
KP3f1 = function_handle(KP3,'vars', [B Kb U3 V3 l1 l2 L1 L2 C2 C3]);

KP1f = @(q) KP1f1(q, U1v, V1v,l1v, l2v, L1v, L2v, C2v, C3v);
KP2f = @(A,B,q) KP2f1(A,q,Kf(A,B,q)(1), U2v, V2v,l1v, l2v, L1v, L2v, C2v, C3v);
KP3f = @(A,B,q) KP3f1(B, Kf(A,B,q)(2), U3v, V3v,l1v, l2v, L1v, L2v, C2v, C3v);

LP1f1 = function_handle(LP1,'vars', [q U1 V1 l1 l2 L1 L2 C2 C3]);
LP2f1 = function_handle(LP2,'vars', [A q Ka La U2 V2 l1 l2 L1 L2 C2 C3]);
LP3f1 = function_handle(LP3,'vars', [B Kb Lb U3 V3 l1 l2 L1 L2 C2 C3]);

LP1f = @(q) LP1f1(q, U1v, V1v,l1v, l2v, L1v, L2v, C2v, C3v);
LP2f = @(A,B,q) LP2f1(A,q,Kf(A,B,q)(1),Lf(A,B,q)(1), U2v, V2v,l1v, l2v, L1v, L2v, C2v, C3v);
LP3f = @(A,B,q) LP3f1(B, Kf(A,B,q)(2),Lf(A,B,q)(2), U3v, V3v,l1v, l2v, L1v, L2v, C2v, C3v);

%global functions = {Ff, iJf, Kf, Lf, P1f, KP1f, LP1f, P2f, KP2f, LP2f, P3f, KP3f, LP3f}

save functions.mat Ff iJf Kf Lf P1f KP1f LP1f P2f KP2f LP2f P3f KP3f LP3f Ff1 iJf1 Kf1 Lf1 P1f1 KP1f1 LP1f1 P2f1 KP2f1 LP2f1 P3f1 KP3f1 LP3f1