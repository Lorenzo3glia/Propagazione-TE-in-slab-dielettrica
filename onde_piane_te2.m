function Etot = onde_piane_te2(x,z,kx,kzin_p,kzrif_p,kztra1,kztra3,Einc_p,E1,Erif_p,E3,d)
%Sto calcolando i campi elettrici
%% Costruisco Gli If

if z<-1*d
    Etot=E1*exp(-1*i*(kx*x+kztra1*z));
elseif (-1*d<z)&&(z<d)
    Etot=Einc_p*exp(-1*i*((kx*x)+(kzin_p*z)))+Erif_p*exp(-1*i*(kx*x + kzrif_p*z));
elseif z>=d
    Etot=E3*exp(-1*i*(kx*x+kztra3*z));
end