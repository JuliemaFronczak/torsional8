function sinal
clear all
close all
clc
%pkg load signal

desloc=load('coefficient.dat','-ascii'); 

[npontos, tam2] = size(desloc)
coef =6;  % (=2 (Cd) =4 (Cl) =6(Cm) )
Vr=2; %AQUI!!!!
f=1/(Vr*8); %Hz
omega=2*3.1416*f
T=1/f
  ttrans=14*T  %tempo inicial
 Tfinal=24*T %tempo final
  
k =1;


 desloc(:,coef)=desloc(:,coef)*k;
%fprintf ('Intervalo de tempo da resposta =');
%dint = input ('Intervalo de tempo da resposta =');
dint=.01;
icont    = 0;
deltat   = 0;
ds       = 0;

sinal_o=desloc(:,coef);
t_o=desloc(:,1);




plot(t_o(10:npontos),sinal_o(10:npontos));
title('Sinal original x sinal delta t constante');
hold on;                                    
%desloc(:,coef)=desloc(:,coef)*k;
for i=1:npontos;
    deltata = deltat;
    dxa = desloc(i,coef);
    deltat=desloc(i,1);
    while (icont*dint + dint) < deltat
        if deltat == (icont + 1)*dint
            icont = icont + 1;
            saida=[icont*dint desloc(i,coef)];
            bfil(icont,1)=saida(1,1);
            bfil(icont,2)=saida(1,2);
        else
            if deltat > (icont + 1)*dint
                icont = icont + 1;
                dix = (dint*icont - deltata)*...
                    (desloc(i,coef) - dxa)/(deltat - deltata)+ dxa;
                saida= [icont*dint dix];
                bfil(icont,1)=saida(1,1);
                bfil(icont,2)=saida(1,2);

            end
        end
    end

end
t=bfil(:,1);
b=bfil(:,2);

[npontos1, tam21] = size(b)
plot(t(10:npontos1),b(10:npontos1))


  ans=size(bfil);
  np=ans(1);[npontos, tam2] = size(desloc)
  icont=0;
   icont2=0;                          
  for i=1:np;
      if bfil(i,1) <= ttrans;
          icont=icont+1;
      end
           if bfil(i,1) <= Tfinal;
          icont2=icont2+1;
      end
  end
  np=icont2;
  bfil2(:,1)=bfil(icont:np,1)-ttrans;
  %save cl.dat bfil -ascii;
  media=mean(bfil(icont:np,2))
  bfil2(:,2)=bfil(icont:np,2)-media;
  mediacl=media;
    n=length(bfil2(:,2));
  rmscl=norm(bfil2(:,2))/sqrt(n);

 x=bfil2(:,1);
 y=bfil2(:,2);
 sy=size(y)

 d=2*2*3.1416/180*sin(omega*x);
 figure
 plot(x,d) 
 hold on
 
 
 plot(x,y)
 title('Deslocamento x sinal com dt constante')
 
 figure
 [a,b,yfit] = Fseries(x,y,100);
% Evaluate on finer grid
save cmVr2.dat bfil2 -ascii
% Visualize results

%plot(x,y,'x',x,yfit)
plot(x,yfit)
hold on
plot(x,y)
title('Sinal com delta t constante x sinal da serie de Fourier (yfit)')
figure
plot(a)
figure
plot(b)
ma=max(abs(a))
mb=max(abs(b))





Fs=1/dint;

t=bfil2(:,1);

xx=bfil2(:,2);


%Figure(2);
i=682*2;
window=i;
figure
[Pxx,f] =pwelch(xx,window,.5,[],Fs,'onesided');

plot(f,Pxx)
amp=max(Pxx)*1.2;  
 axis([0 2 0 amp]);
 espec=[f Pxx];

