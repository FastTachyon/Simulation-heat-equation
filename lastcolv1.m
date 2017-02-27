function [a]=lastcol(N,d)
%Fonction pour les deux dernière colonne 
% Paramètres d'entrée:
% N est la dimension de la matrice de températures
% d est le vecteur représentant les différents dx
% Paramètre de sortie
% a est une matrices remplies de zeros, mis-à-partles deux dernière
% colonnes
%  Ordre des colonnes:
% T(i,j)    []
% T(i,j+1)  []
% ...
% T(i+1,j)  []
% T(i+1,j+1)[]
% Auteurs: Mathieu Bergeron
% Fait le 18 février 2017

%% Section Temporaire
[~,c21]=T10000T;
[~,c22]=T01000T;
[~,c23]=T00100T;
[c14,c24]=T00010T;
[c15,c25]=T00001T;
% Vecteur dr
 inter=(1*10^-9:(borne)^(1/4)/n :(borne)^(1/4)).^4;
 z=inter;
 r=inter;
%% Définition de la matrice
a=zeros(N^2);

%% Avant dernière colonne
 % Terme général de r
 vr=k/2*[1/r(end-1)*c14{1}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))+c24{1}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))...
    1/r(end-1)*c14{2}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))+c24{2}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))...
    1/r(end-1)*c14{3}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))+c24{3}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))...
    1/r(end-1)*c14{4}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))+c24{4}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))...
    1/r(end-1)*c14{5}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))+c24{5}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))];

% Première rangée
 % en z
a(N^2-2*N+1,(N^2-2*N+1):(N^2-2*N+5))=k/2*[c21{1}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))...
    ,c21{2}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)),c21{3}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))...
    ,c21{4}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)),c21{5}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))];
 % en r
a(N^2-2*N+1,N^2-2*N+1-3*N:N:(N^2-2*N+N+1))=vr;
% point milieu
a(N^2-2*N+1,N^2-2*N+1)=k/2*(vr(4)+c21{1}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)))-rho*Cp/dt;
% Deuxième rangée
 % en z
a(N^2-2*N+2,(N^2-2*N+1):(N^2-2*N+3))=k/2*[c22{1}(z(1)-z(2),z(3)-z(2),z(4)-z(1),z(5)-z(2))...
    ,c22{2}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)),c22{3}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2))...
    ,c22{4}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)),c22{5}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2))];
 % en r
a(N^2-2*N+2,(N^2-2*N-3*N+2):N:(N^2+2*N+N+2))=vr;
% point milieu
a(N^2-2*N+2,N^2-2*N+2)=k/2*(vr(4)+c22{2}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)))-rho*Cp/dt;
%Centre
for i=3:N^2-3
     % en z
  a(N^2-2*N+i,(N^2-2*N+i-2):(N^2-2*N+i+2))=k/2*[c23{1}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))...
    ,c23{2}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)),c23{3}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))...
    ,c23{4}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)),c23{5}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))];
     % en r
  a(N^2-2*N+i,(N^2-2*N-3*N+i):N:(N^2-2*N+N+i))=vr;
  % point milieu
  a(N^2-2*N+i,N^2-2*N+i)=k/2*(vr(4)+c23{3}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)))-rho*Cp/dt;
end
%Avant dernière rangée
 % en z
a(N^2-N+1,(N^2-N+4):(N^2-N))=k/2*[c24{1}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))...
    ,c24{2}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)),c24{3}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))...
    ,c24{4}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)),c24{5}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))];
 % en r
a(N^2-N+1,(N^2-3*N+1):N:(N^2-N+N+1))=vr;
% point milieu
a(N^2-N+1,N^2-N+1)=k/2*(vr(4)+c24{4}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)))-rho*Cp/dt;
%Dernière rangée
 % en z
a(N^2-N,(N^2-N+4):(N^2-N))=k/2*[c25{1}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))...
    ,c25{2}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)),c25{3}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))...
    ,c25{4}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)),c25{5}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))];
 % en r
a(N^2-N,(N^2-N-3*N):N:(N^2-N+N))=vr;
% point milieu
a(N^2-N,N^2-N)=k/2*(vr(4)+c25{5}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)))-rho*Cp/dt;
%% Dernière colonne
 % Terme général de r
 vr=k/2*[1/r(end-1)*c15{1}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))+c25{1}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))...
    1/r(end-1)*c15{2}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))+c25{2}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))...
    1/r(end-1)*c15{3}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))+c25{3}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))...
    1/r(end-1)*c15{4}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))+c25{4}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))...
    1/r(end-1)*c15{5}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))+c25{5}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))];

% Première rangée
 % en z
a(N^2-N+1,(N^2-N+1):(N^2-N+5))=k/2*[c21{1}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))...
    ,c21{2}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)),c21{3}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))...
    ,c21{4}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)),c21{5}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))];
 % en r
a(N^2-N+1,N^2-N+1-3*N:N:(N^2-N+N+1))=k/2*vr;
 % point milieu
a(N^2-N+1,N^2-N+1)=k/2*(vr(5)+c21{1}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)))-rho*Cp/dt;
% Deuxième rangée
 % en z
a(N^2-N+2,(N^2-2*N+1):(N^2-2*N+3))=[c22{1}(z(1)-z(2),z(3)-z(2),z(4)-z(1),z(5)-z(2))...
    ,c22{2}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)),c22{3}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2))...
    ,c22{4}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)),c22{5}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2))];
 % en r
a(N^2-N+2,(N^2-2*N-3*N+2):N:(N^2+2*N+N+2))=k/2*vr;
 % point milieu
a(N^2-N+2,N^2-N+2)=k/2*(vr(5)+c22{2}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)))-rho*Cp/dt;
%Centre
for i=3:N^2-3
     % en z
  a(N^2-N+i,(N^2-N+i-2):(N^2-N+i+2))=k/2*[c23{1}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))...
    ,c23{2}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)),c23{3}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))...
    ,c23{4}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)),c23{5}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))];
     % en r
  a(N^2-N+i,(N^2-N-3*N+i):N:(N^2-N+N+i))=k/2*vr;
     % point milieu
  a(N^2-N+i,N^2-N+i)=k/2*(vr(5)+c23{3}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)))-rho*Cp/dt;
end
%Avant dernière rangée
 % en z
a(N^2-1,(N^2-4):(N^2))=k/2*[c25{1}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))...
    ,c25{2}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)),c25{3}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))...
    ,c25{4}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)),c25{5}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))];
 % en r
a(N^2-1,(N^2-4*N-1):N:(N^2-1))=k/2*vr;
 % point milieu
a(N^2-1,N^2-1)=k/2*(vr(5)+c25{4}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1)))-rho*Cp/dt;
%Dernière rangée
 % en z
a(N^2,(N^2-4):(N^2))=k/2*[c25{1}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))...
    ,c25{2}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)),c25{3}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))...
    ,c25{4}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)),c25{5}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))];
 % en r
a(N^2,(N^2-5*N):N:(N^2))=k/2*vr;
 % point milieu
a(N^2,N^2)=kr/2*(vr(5)+c25{5}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)))-rho*Cp/dt;




















