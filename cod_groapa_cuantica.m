% Valori si functii proprii ale hamiltonianului in groapa cuantica
% dreptunghiulara (SQW) cu pereti finiti
clc; clear; close all;
h=6.626*1e-34; % J*s; constanta Planck
hbar=h/2/pi; % constanta Planck redusa
eV=1.602*1e-19; % J; electron-Volt (unitate atomica de energie)
m0=9.1*1e-31; % kg; masa electronului in vid
meff=0.067*m0; % masa efectiva a electronului de conductie in GaAs
U0=0.280*eV; % bariera de energie potentiala GaAs/AlGaAs
s=zeros(1000,20);
c1=zeros(20,20);
c2=c1;
c3=c1;
i=1;
j=1;
max=1;
Emax=U0; NE=50000;
for l=5:5:100
w=l*1e-9; % m; largimea gropii de potential
xs=-0.5*w; xd=0.5*w; Nx=1000; x=linspace(xs,xd,Nx); dx=x(2)-x(1); % discretizarea x
a=2*meff*dx^2/hbar^2; % vezi Curs 5 partea a II-a, rel. (***)
U=zeros(1,Nx); U(x<-w/2)=U0; U(x>w/2)=U0; % profilul de energie potentiala
phi=zeros(1,Nx); % prealocare functie de unda
phi(2)=1e-6; % pentru a evita o solutie identic nula

 
if l>40
    Emax=max;
end
E=Emax*logspace(-2,0,NE); % progresie exponentiala
phid=zeros(1,NE);

for iE=1:NE % stabileste energia de tir
    for ix=2:Nx-1 % implementeaza recurenta Schrodinger
        phi(ix+1)=2*phi(ix)-phi(ix-1)+a*(U(ix)-E(iE))*phi(ix); % Curs 5
    end
    phid(iE)=phi(Nx); % valoarea functiei de unda la frontiera dreapta
end
logphi=log(abs(phid)); % logaritmare pentru postselectie mai rapida
contor=1; Ep=zeros(1,1);
for iE=2:NE-1
    if (logphi(iE-1)>logphi(iE))&&(logphi(iE)<logphi(iE+1)) % test minim
        Ep(contor)=E(iE); contor=contor+1;
    end
end
NEp=length(Ep); % numarul starilor proprii
for iEp=1:NEp % ciclul energiilor proprii
    c1(i,iEp)=Ep(iEp);
    
end
    max=Ep(NEp);
i = i + 1;
end


i=9;
Emin=0; Emax=U0; NE=50000;
for l=45:5:100
w=l*1e-9; % m; largimea gropii de potential
xs=-0.5*w; xd=0.5*w; Nx=1000; x=linspace(xs,xd,Nx); dx=x(2)-x(1); % discretizarea x
a=2*meff*dx^2/hbar^2; % vezi Curs 5 partea a II-a, rel. (***)
U=zeros(1,Nx); U(x<-w/2)=U0; U(x>w/2)=U0; % profilul de energie potentiala
phi=zeros(1,Nx); % prealocare functie de unda
phi(2)=1e-6; % pentru a evita o solutie identic nula
if l> 85
Emax=max-0.01e-19;
end
E=Emax*logspace(-2,0,NE); % progresie exponentiala
phid=zeros(1,NE);

for iE=1:NE % stabileste energia de tir
    for ix=2:Nx-1 % implementeaza recurenta Schrodinger
        phi(ix+1)=2*phi(ix)-phi(ix-1)+a*(U(ix)-E(iE))*phi(ix); % Curs 5
    end
    phid(iE)=phi(Nx); % valoarea functiei de unda la frontiera dreapta
end
logphi=log(abs(phid)); % logaritmare pentru postselectie mai rapida
contor=1; Ep=zeros(1,1);
for iE=2:NE-1
    if (logphi(iE-1)>logphi(iE))&&(logphi(iE)<logphi(iE+1)) % test minim
        Ep(contor)=E(iE); contor=contor+1;
    end
end
NEp=length(Ep); % numarul starilor proprii
for iEp=1:NEp % ciclul energiilor proprii
    if(iEp+1<20)
    c2(i,iEp+1)=Ep(iEp);
    end
    
end
    max=Ep(NEp);
i = i + 1;

end

Emax=U0;i=19;
for l=95:5:100
w=l*1e-9; % m; largimea gropii de potential
xs=-0.5*w; xd=0.5*w; Nx=1000; x=linspace(xs,xd,Nx); dx=x(2)-x(1); % discretizarea x
a=2*meff*dx^2/hbar^2; % vezi Curs 5 partea a II-a, rel. (***)
U=zeros(1,Nx); U(x<-w/2)=U0; U(x>w/2)=U0; % profilul de energie potentiala
phi=zeros(1,Nx); % prealocare functie de unda
phi(2)=1e-6; % pentru a evita o solutie identic nula
E=Emax*logspace(-2,0,NE); % progresie exponentiala
phid=zeros(1,NE);

for iE=1:NE % stabileste energia de tir
    for ix=2:Nx-1 % implementeaza recurenta Schrodinger
        phi(ix+1)=2*phi(ix)-phi(ix-1)+a*(U(ix)-E(iE))*phi(ix); % Curs 5
    end
    phid(iE)=phi(Nx); % valoarea functiei de unda la frontiera dreapta
end
logphi=log(abs(phid)); % logaritmare pentru postselectie mai rapida
contor=1; Ep=zeros(1,1);
for iE=2:NE-1
    if (logphi(iE-1)>logphi(iE))&&(logphi(iE)<logphi(iE+1)) % test minim
        Ep(contor)=E(iE); contor=contor+1;
    end
end
NEp=length(Ep); % numarul starilor proprii
for iEp=1:NEp % ciclul energiilor proprii
    c3(i,iEp)=Ep(iEp);
end

    max=Ep(NEp);

i = i + 1;

end
disp(c1);
c2(:,1:8)=0;
c3(:,1:18)=0;
disp(c2);
disp(c3);
c1=c1+c2+c3;
disp(c1);
figure(1);
for i=1:9
 
    plot(i*5:5:100,c1(i:end,i),'-'); hold on;
    
end

for i=10:18
 
    plot(i*5:5:100,c1(i-1:end-1,i),'-'); hold on;
    
end

for i=19:20
 
    plot(i*5:5:100,c1(i:end,i),'-'); hold on;
    
end
xlabel('Largimea SQW ( nm )'); ylabel('Energie electron ( J )'); grid
%legend('Nivel1','Nivel2','Nivel3','Nivel4','Nivel5','Nivel6','Nivel7','Nivel8','Nivel9','Nivel10 '...
 %   ,'Nivel11','Nivel12','Nivel13','Nivel14','Nivel15','Nivel16','Nivel17','Nivel18','Nivel19','Nivel20');

% Se observa ca odata cu cresterea largimii creste si numarul de nivele
% ale energiei , ajungand la nivele mai mari de energie .
% La primul nivel de energie, pe care toate gropile de latimi diferite il au, se
% observa faptul ca este cel mai "ameloriat", datorita faptului ca la
% largimi mari valorile se stabilizeaza.
% Cu cat largimea creste , apar nivele noi cu energii mai mari ale caror
% oscilatii sunt "agresive".
% De asemenea cu cresterea largimii si a energiei se observa si cresterea
% frecventei energiei.

    








