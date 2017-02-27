function [vec]=vecB(N)
% Fonction pour les deux dernière colonne 
% Paramètres d'entrée:
% N est la dimension de la matrice de températures
% T est la matrice des températures NxN au temps n 
% Paramètre de sortie
% a est une matrices remplies de zeros, mis-à-partles deux dernière
% colonnes
%  Ordre des colonnes:
% T(i+1,j)    []
% T(i+2,j)  []
% ...
% T(N,j+1)  []
% T(N+1,j+1)[]
% Auteurs: Mathieu Bergeron
% Fait le 25 février 2017
%% Section temporaire
[d11,d21]=h10000h;
[d12,d22]=h01000h;
[d13,d23]=h00100h;
[d14,d24]=h00010h;
[d15,d25]=h00001h;
% Vecteur dr
borne=4*2.1929*10^-7;
 inter=(0:(borne)^(1/4)/N :(borne)^(1/4)).^4;
 z=inter;
 r=inter;
 k=1;
 Cp=1;
 dt=1;
 rho=1;
% Température fake
T=ones(N)*25;
%% Section du calcul de B
vec=zeros(1,N^2);
%% Colonne de droite 
 % Première
    a=r(2)-r(1);
    b=r(3)-r(1);
    c=r(4)-r(1);
    d=r(5)-r(1);
    z1=z(2)-z(1);
    z2=z(3)-z(1);
    z3=z(4)-z(1);
    z4=z(5)-z(1);
    vec(1)=-k/2*(d11(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        r(1)*d21(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        d21(z1,z2,z3,z4,T(1),T(2),T(3),T(4),T(5)))-r(1)*Cp*rho/dt*T(1);
 %Deuxième
    z1=z(1)-z(2);
    z2=z(3)-z(2);
    z3=z(4)-z(2);
    z4=z(5)-z(2);
    vec(2)=-k/2*(d11(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        r(1)*d21(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        d22(z1,z2,z3,z4,T(1),T(2),T(3),T(4),T(5)))-r(1)*Cp*rho/dt*T(1);
 %Centre
    z1=z(1:N-4)-z(3:N-2);
    z2=z(2:N-3)-z(3:N-2);
    z3=z(4:N-1)-z(3:N-2);
    z4=z(5:N)-z(3:N-2);
    vec(3:N-2)=-k/2*(d11(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        r(1)*d21(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        d23(z1,z2,z3,z4,T(1:N-4),T(2:N-3),T(3:N-2),T(4:N-1),T(5:N)))-r(1)*Cp*rho/dt*T(1);
 % Avant dernière
    z1=z(N-4)-z(N-1);
    z2=z(N-3)-z(N-1);
    z3=z(N-2)-z(N-1);
    z4=z(N)-z(N-1);
    vec(N-1)=-k/2*(d11(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        r(1)*d21(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        d24(z1,z2,z3,z4,T(1),T(2),T(3),T(4),T(5)))-r(1)*Cp*rho/dt*T(1);
 % Dernière
    z1=z(N-4)-z(N);
    z2=z(N-3)-z(N);
    z3=z(N-2)-z(N);
    z4=z(N-1)-z(N);
    vec(N)=-k/2*(d11(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        r(1)*d21(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        d25(z1,z2,z3,z4,T(1),T(2),T(3),T(4),T(5)))-r(1)*Cp*rho/dt*T(1); 
%% Colonne de droite -1
 % Première
    a=r(1)-r(2);
    b=r(3)-r(2);
    c=r(4)-r(2);
    d=r(5)-r(2);
    z1=z(2)-z(1);
    z2=z(3)-z(1);
    z3=z(4)-z(1);
    z4=z(5)-z(1);
    vec(N+1)=-k/2*(d12(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        r(2)*d22(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        d21(z1,z2,z3,z4,T(N+1),T(N+2),T(N+3),T(N+4),T(5)))-r(2)*Cp*rho/dt*T(1);
 %Deuxième
    z1=z(1)-z(2);
    z2=z(3)-z(2);
    z3=z(4)-z(2);
    z4=z(5)-z(2);
    vec(N+2)=-k/2*(d12(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        r(2)*d22(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        d21(z1,z2,z3,z4,T(1),T(2),T(3),T(4),T(5)))-r(2)*Cp*rho/dt*T(1);
 %Centre
    z1=z(1:N-4)-z(3:N-2);
    z2=z(2:N-3)-z(3:N-2);
    z3=z(4:N-1)-z(3:N-2);
    z4=z(5:N)-z(3:N-2);
    vec(N+3:2*N-2)=-k/2*(d12(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        r(2)*d22(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        d21(z1,z2,z3,z4,T(1),T(2),T(3),T(4),T(5)))-r(2)*Cp*rho/dt*T(1);
 % Avant dernière
    z1=z(N-4)-z(N-1);
    z2=z(N-3)-z(N-1);
    z3=z(N-2)-z(N-1);
    z4=z(N)-z(N-1);
    vec(2*N-1)=-k/2*(d12(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        r(2)*d22(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        d21(z1,z2,z3,z4,T(1),T(2),T(3),T(4),T(5)))-r(2)*Cp*rho/dt*T(1);
 % Dernière
    z1=z(N-4)-z(N);
    z2=z(N-3)-z(N);
    z3=z(N-2)-z(N);
    z4=z(N-1)-z(N);
    vec(2*N)=-k/2*(d12(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        r(2)*d22(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        d21(z1,z2,z3,z4,T(1),T(2),T(3),T(4),T(5)))-r(2)*Cp*rho/dt*T(1);
%% Centre
for i=2N+1:N:N^2-3*N
    fac=(i-1)/N+1;
% Première
    a=r(fac-2)-r(fac);
    b=r(fac-1)-r(fac);
    c=r(fac+1)-r(fac);
    d=r(fac+2)-r(fac);
    z1=z(2)-z(1);
    z2=z(3)-z(1);
    z3=z(4)-z(1);
    z4=z(5)-z(1);
    vec(N+1)=-k/2*(d12(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        r(fac)*d22(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        d21(z1,z2,z3,z4,T(1),T(2),T(3),T(4),T(5)))-r(fac)*Cp*rho/dt*T(1);
 %Deuxième
    z1=z(1)-z(2);
    z2=z(3)-z(2);
    z3=z(4)-z(2);
    z4=z(5)-z(2);
    vec(N+2)=-k/2*(d12(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        r(2)*d22(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        d21(z1,z2,z3,z4,T(1),T(2),T(3),T(4),T(5)))-r(2)*Cp*rho/dt*T(1);
 %Centre
    z1=z(1:N-4)-z(3:N-2);
    z2=z(2:N-3)-z(3:N-2);
    z3=z(4:N-1)-z(3:N-2);
    z4=z(5:N)-z(3:N-2);
    vec(N+3:2*N-2)=-k/2*(d12(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        r(2)*d22(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        d21(z1,z2,z3,z4,T(1),T(2),T(3),T(4),T(5)))-r(2)*Cp*rho/dt*T(1);
 % Avant dernière
    z1=z(N-4)-z(N-1);
    z2=z(N-3)-z(N-1);
    z3=z(N-2)-z(N-1);
    z4=z(N)-z(N-1);
    vec(2*N-1)=-k/2*(d12(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        r(2)*d22(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        d21(z1,z2,z3,z4,T(1),T(2),T(3),T(4),T(5)))-r(2)*Cp*rho/dt*T(1);
 % Dernière
    z1=z(N-4)-z(N);
    z2=z(N-3)-z(N);
    z3=z(N-2)-z(N);
    z4=z(N-1)-z(N);
    vec(2*N)=-k/2*(d12(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        r(2)*d22(a,b,c,d,T(1),T(1+N),T(1+2*N),T(1+3*N),T(1+4*N))+...
        d21(z1,z2,z3,z4,T(1),T(2),T(3),T(4),T(5)))-r(2)*Cp*rho/dt*T(1);
end
%% Colonne de gauche -1

%% Colonne de gauche
   

