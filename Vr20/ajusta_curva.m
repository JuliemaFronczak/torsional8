% Teste da rotina fit do Matlab para uma fun��o de ajuste qualquer

clc, clear, close all

% Carrega dados:

load cmVr20.mat

% Tipo de fun��o

% N�o se esquecer de modificar a frequencia para cada sinal!!! Freq. em rad/s

ft = fittype('rho*cos(0.03927*x+theta)','independent','x');

% Ajusta com par�metros iniciais:

[model,gof,output]= fit(x,y,ft,'StartPoint',[1 1]);

disp(model)
disp(gof)

% Plota a curva original e a ajustada

plot(x,y), hold on
plot(model,'r')

syms a b

eqns = [a^2 + b^2 == model.rho^2, b == -a*tan(model.theta)];
S = solve(eqns,[a b]);
a = eval(S.a); 
a = a(1); disp(a)
b = eval(S.b); 
b = b(1); disp(b)