clear all;
close all;

%% Dati Caratteristiche dei mezzi
epsilon0=8.854e-12;
mu0=4*pi*1.e-7;

epsilon = 1*epsilon0;
epsilon_p = 2*epsilon0;

mu = 1*mu0;
mu_p = 1*mu0;

f=1e9;
c0=1/sqrt(epsilon0*mu0);
w = 2*pi*f;

lambda=c0/f;
d=0.7*lambda;

Einc_p = 1;
n12 = sqrt((epsilon*mu)/(epsilon_p*mu_p));

%% Calcolo del numero d'onda

k=w.*sqrt(epsilon*mu);
k_p=w.*sqrt(epsilon_p*mu_p);

%% Punti di valutazione
nptx=101;
nptz=10001;

x=linspace(-2*lambda,2*lambda,nptx);
z=linspace(-2*lambda,2*lambda,nptz);

%% Calcolo dei punti di intersezione
%% Per i Dispari
%al^2+((mu^2/mu_p^2)*al^2*tan(al)^2)-(k_p*d)^2*(1-n12^2);
%% Per i pari
%al^2+((mu^2/mu_p^2)*al^2*atan(al)^2)-(k_p*d)^2*(1-n12^2);

N=ceil(k_p*d*sqrt(1-n12^2)/(pi/2))

for I = 1:N   
    if mod(I, 2) == 0  % Se I è pari
        alpha_v(I) = findzero('alp', pi*(I-1)/2, pi*I/2,1e-3,mu,mu_p,k_p,d,n12);
    else               % Se I è dispari
        alpha_v(I) = findzero('ald', pi*(I-1)/2, pi*I/2,1e-3,mu,mu_p,k_p,d,n12);
    end
end

disp('Alpa_v')
disp(alpha_v(:))

%% Valutazione delle coordinate di alpha e beta

raggio=k_p*d*sqrt(1-n12^2);
beta_v=sqrt(raggio^2 - alpha_v.^2);

disp('Beta_v')
disp(beta_v(:))

%% Selezione di alpha e beta per i vettori di propagazione
alpha=alpha_v(3);
beta=beta_v(3);

kzin_p=alpha/d;
kzrif_p=-kzin_p;
kz=-1*i*(beta/d);

kztra3=kz;
kztra1=-kz;

kx=sqrt(k_p^2+kzin_p^2)


zita=sqrt(mu/epsilon)
zita_p=sqrt(mu_p/epsilon_p)
z_e=(zita*k)/kz             % Z mezzo esterno
z_p=(zita_p*k_p)/kzin_p     %Z mezzo interno da controllare il valore di kzinp

%% Calcolo dei valori di campo

E1=(2*z_e/(z_e-z_p))*Einc_p*exp(1*i*(kz+kzin_p)*d);
E3=(2*z_e/(z_e+z_p))*Einc_p*exp(1*i*(kz-kzin_p)*d);
Erif_p=((z_e+z_p)/(z_e-z_p))*Einc_p*exp(2*1*i*kzin_p*d);

%% Calcolo del campo totale

for xn = 1:nptx
    for zn = 1:nptz
        Etot(xn,zn)=onde_piane_te2(x(xn),z(zn),kx,kzin_p,kzrif_p,kztra1,kztra3,Einc_p,E1,Erif_p,E3,d);
    end
end
figure(1)
imagesc(z,x,abs(Etot));
colorbar;

%% Plot del'arco di circonferenza
% Angoli tra 0 e pi/2
theta = linspace(0, pi/2, 10000); %Ho aumentato il numero di punti per poter trovare l'intersezione in maniera più agevole

% Equazioni parametriche della circonferenza
x = raggio * cos(theta);
y = raggio * sin(theta);

% Plot della circonferenza
figure;
plot(x, y, 'b');
hold on 

%% Plot con fplot delle funzioni ald e alp

ylim([0,raggio])
xlim([0,raggio])
betd=@(ald) (mu/mu_p)*ald.*tan(ald);
fplot(betd)
hold on
betp=@(alp) -(mu/mu_p)*alp.*cot(alp);
fplot(betp)
hold off

%% Plot della sezione trasversa 
figure(3);
hold on
plot(z,abs(Etot(1,:)));
plot(linspace(d,d,101),linspace(0,2,101),'.')
plot(linspace(-d,-d,101),linspace(0,2,101),'.')
hold off