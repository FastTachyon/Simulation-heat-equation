function [a]=lastcol(N)
%Fonction servant à la création de la matrice M'
% Paramètres d'entrée:
% N est la dimension de la matrice de températures
% d est le vecteur représentant les différents dx
% Paramètre de sortie
% a est la matrice servant à être inverser pour résoudre le problème.
%  Ordre des colonnes:
% T(i,j)    []
% T(i,j+1)  []
% ...
% T(i+1,j)  []
% T(i+1,j+1)[]
% Auteurs: Mathieu Bergeron
% Fait le 18 février 2017

%% Section Temporaire
[c11,c21]=T10000T;
[c12,c22]=T01000T;
[c13,c23]=T00100T;
[c14,c24]=T00010T;
[c15,c25]=T00001T;
% Vecteur dr
borne=4*2.1929*10^-7;
 inter=(0:(borne)^(1/4)/N :(borne)^(1/4)).^4;
 z=inter;
 r=inter;
 k=1;
 Cp=1;
 dt=1;
 rho=1;
 tic;
%% Définition de la matrice
a=zeros(N^2);
%% Première colonne
 % Terme général de r
 ap=r(2)-r(1);
 bp=r(3)-r(1);
 cp=r(4)-r(1);
 dp=r(5)-r(1);
 vr=k/2*[1*c11{1}(ap,bp,cp,dp)+r(2)*c21{1}(ap,bp,cp,dp)...
    1*c11{2}(ap,bp,cp,dp)+r(2)*c21{2}(ap,bp,cp,dp)...
    1*c11{3}(ap,bp,cp,dp)+r(2)*c21{3}(ap,bp,cp,dp)...
    1*c11{4}(ap,bp,cp,dp)+r(2)*c21{4}(ap,bp,cp,dp)...
    1*c11{5}(ap,bp,cp,dp)+r(2)*c21{5}(ap,bp,cp,dp)];
% Première rangée
 % en z
a(1,(1):(5))=r(1)*k/2*[c21{1}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))...
    ,c21{2}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)),c21{3}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))...
    ,c21{4}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)),c21{5}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))];
 % en r
a(1,1:N:(4*N+1))=vr;
% point milieu
a(1,1)=k/2*(vr(1)+r(1)*c21{1}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)))-r(1)*rho*Cp/dt;
% Deuxième rangée
 % en z
a(2,1:5)=r(1)*k/2*[c22{1}(z(1)-z(2),z(3)-z(2),z(4)-z(1),z(5)-z(2))...
    ,c22{2}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)),c22{3}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2))...
    ,c22{4}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)),c22{5}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2))];
 % en r
a(2,2:N:(4*N+2))=vr;
% point milieu
a(2,2)=k/2*(vr(1)+r(1)*c22{2}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)))-r(1)*rho*Cp/dt;
%Centre
for i=3:N-2
     % en z
  a(i,(i-2):(i+2))=r(1)*k/2*[c23{1}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))...
    ,c23{2}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)),c23{3}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))...
    ,c23{4}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)),c23{5}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))];
     % en r
  a(i,(i):N:(4*N+i))=k/2*vr;
  % point milieu
  a(i,i)=k/2*(vr(1)+r(1)*c23{3}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)))-r(1)*rho*Cp/dt;
end
%Avant dernière rangée
 % en z
a(N-1,(N-4):(N))=r(1)*k/2*[c24{1}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))...
    ,c24{2}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)),c24{3}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))...
    ,c24{4}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)),c24{5}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))];
 % en r
a(N-1,(N-1):N:(N-1+4*N+1))=k/2*vr;
% point milieu
a(N-1,N-1)=k/2*(vr(1)+r(1)*c24{4}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)))-r(1)*rho*Cp/dt;
%Dernière rangée
 % en z
a(N,(N-4):(N))=r(1)*k/2*[c25{1}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))...
    ,c25{2}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)),c25{3}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))...
    ,c25{4}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)),c25{5}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))];
 % en r
a(N,N:N:(5*N))=k/2*vr;
% point milieu
a(N,N)=k/2*(vr(1)+r(1)*c25{5}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)))-r(1)*rho*Cp/dt;

%% Deuxième colonne
 % Terme général de r
 ap=r(1)-r(2);
 bp=r(3)-r(2);
 cp=r(4)-r(2);
 dp=r(5)-r(2);
 vr=k/2*[1*c12{1}(ap,bp,cp,dp)+r(2)*c22{1}(ap,bp,cp,dp)...
    1*c12{2}(ap,bp,cp,dp)+r(2)*c22{2}(ap,bp,cp,dp)...
    1*c12{3}(ap,bp,cp,dp)+r(2)*c22{3}(ap,bp,cp,dp)...
    1*c12{4}(ap,bp,cp,dp)+r(2)*c22{4}(ap,bp,cp,dp)...
    1*c12{5}(ap,bp,cp,dp)+r(2)*c22{5}(ap,bp,cp,dp)];

% Première rangée
 % en z
a(N+1,(N+1):(N+5))=r(2)*k/2*[c21{1}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))...
    ,c21{2}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)),c21{3}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))...
    ,c21{4}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)),c21{5}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))];
 % en r
a(N+1,(1):N:(4*N+1))=k/2*vr;
 % point milieu
a(N+1,N+1)=k/2*(vr(2)+r(2)*c21{1}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)))-r(2)*rho*Cp/dt;
% Deuxième rangée
 % en z
a(N+2,(N+1):(N+5))=k/2*r(2)*[c22{1}(z(1)-z(2),z(3)-z(2),z(4)-z(1),z(5)-z(2))...
    ,c22{2}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)),c22{3}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2))...
    ,c22{4}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)),c22{5}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2))];
 % en r
a(N+2,(N-N+2):N:(4*N+2))=k/2*vr;
 % point milieu
a(N+2,N+2)=k/2*(vr(2)+r(2)*c22{2}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)))-r(2)*rho*Cp/dt;
%Centre
for i=3:N-2
     % en z
  a(N+i,(N+i-2):(N+i+2))=r(2)*k/2*[c23{1}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))...
    ,c23{2}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)),c23{3}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))...
    ,c23{4}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)),c23{5}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))];
     % en r
  a(N+i,(N-N+i):N:(4*N+i))=k/2*vr;
     % point milieu
  a(N+i,N+i)=k/2*(vr(2)+r(2)*c23{3}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)))-r(2)*rho*Cp/dt;
end
%Avant dernière rangée
 % en z
a(2*N-1,(2*N-4):(2*N))=r(2)*k/2*[c25{1}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))...
    ,c25{2}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)),c25{3}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))...
    ,c25{4}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)),c25{5}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))];
 % en r
a(2*N-1,(N-1):N:(5*N-1))=k/2*vr;
 % point milieu
a(N-1,N-1)=k/2*(vr(2)+r(2)*c25{4}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)))-r(2)*rho*Cp/dt;
%Dernière rangée
 % en z
a(2*N,(2*N-4):(2*N))=r(2)*k/2*[c25{1}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))...
    ,c25{2}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)),c25{3}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))...
    ,c25{4}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)),c25{5}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))];
 % en r
a(2*N,(N):N:(5*N))=k/2*vr;
 % point milieu
a(2*N,2*N)=k/2*(vr(2)+r(2)*c25{5}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)))-r(2)*rho*Cp/dt;
%% Milieu de la matrice
for j=2*N:N:(N^2-3*N)
% Terme général de r
 ap=r(j/N+1-2)-r(j/N+1);
 bp=r(j/N+1-1)-r(j/N+1);
 cp=r(j/N+1+1)-r(j/N+1);
 dp=r(j/N+1+2)-r(j/N+1);
 vr=k/2*[1*c13{1}(ap,bp,cp,dp)+r(j/N+1)*c23{1}(ap,bp,cp,dp)...
    1*c13{2}(ap,bp,cp,dp)+r(j/N+1)*c23{2}(ap,bp,cp,dp)...
    1*c13{3}(ap,bp,cp,dp)+r(j/N+1)*c23{3}(ap,bp,cp,dp)...
    1*c13{4}(ap,bp,cp,dp)+r(j/N+1)*c23{4}(ap,bp,cp,dp)...
    1*c13{5}(ap,bp,cp,dp)+r(j/N+1)*c23{5}(ap,bp,cp,dp)];
% Première rangée
 % en z
a(j+1,(j+1):(j+5))=r(j/N+1)*k/2*[c21{1}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))...
    ,c21{2}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)),c21{3}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))...
    ,c21{4}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)),c21{5}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))];
 % en r
a(j+1,(j-2*N+1):N:(j+2*N+1))=k/2*vr;
% point milieu
a(j+1,j+1)=k/2*(vr(3)+r(j/N+1)*c21{1}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)))-r(j/N+1)*rho*Cp/dt;
% Deuxième rangée
 % en z
a(j+2,(j+1):(j+5))=r(j/N+1)*k/2*[c22{1}(z(1)-z(2),z(3)-z(2),z(4)-z(1),z(5)-z(2))...
    ,c22{2}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)),c22{3}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2))...
    ,c22{4}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)),c22{5}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2))];
 % en r
a(j+2,(j-2*N+2):N:(j+2*N)+2)=k/2*vr;
% point milieu
a(j+2,j+2)=k/2*(vr(3)+r(j/N+1)*c22{2}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)))-r(j/N+1)*rho*Cp/dt;
%Centre
for i=3:N-2
     % en z
     
  a(j+i,(j+i-2):(j+i+2))=r(j/N+1)*k/2*[c23{1}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))...
    ,c23{2}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)),c23{3}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))...
    ,c23{4}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)),c23{5}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))];
     % en r
  a(j+i,(j-2*N+i):N:(j+2*N+i))=k/2*vr;
  % point milieu
  a(j+i,j+i)=k/2*(vr(3)+r(j/N+1)*c23{3}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)))-r(j/N+1)*rho*Cp/dt;
end
%Avant dernière rangée
 % en z
a(j+N-1,(j+N-4):(j+N))=r(j/N+1)*k/2*[c24{1}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))...
    ,c24{2}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)),c24{3}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))...
    ,c24{4}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)),c24{5}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))];
 % en r
a(j+N-1,(j+N-1-2*N):N:(j+N-1+2*N))=k/2*vr;
% point milieu
a(j+N-1,j+N-1)=k/2*(vr(3)+r(j/N+1)*c24{4}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)))-r(j/N+1)*rho*Cp/dt;
%Dernière rangée
 % en z
a(j+N,(j+N-4):(j+N))=r(j/N+1)*k/2*[c25{1}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))...
    ,c25{2}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)),c25{3}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))...
    ,c25{4}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)),c25{5}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))];
 % en r
a(j+N,(j+N-2*N):N:(j+N+2*N))=k/2*vr;
% point milieu
a(j+N,j+N)=k/2*(vr(3)+r(j/N+1)*c25{5}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)))-r(j/N+1)*rho*Cp/dt;
end

%% Avant dernière colonne
 % Terme général de r
 vr=k/2*[1*c14{1}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))+r(end-1)*c24{1}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))...
    1*c14{2}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))+r(end-1)*c24{2}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))...
    1*c14{3}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))+r(end-1)*c24{3}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))...
    1*c14{4}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))+r(end-1)*c24{4}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))...
    1*c14{5}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))+r(end-1)*c24{5}(r(end-4)-r(end-1),r(end-3)-r(end-1),r(end-2)-r(end-1),r(end)-r(end-1))];

% Première rangée
 % en z
a(N^2-2*N+1,(N^2-2*N+1):(N^2-2*N+5))=r(end-1)*k/2*[c21{1}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))...
    ,c21{2}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)),c21{3}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))...
    ,c21{4}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)),c21{5}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))];
 % en r
a(N^2-2*N+1,N^2-2*N+1-3*N:N:(N^2-2*N+N+1))=vr;
% point milieu
a(N^2-2*N+1,N^2-2*N+1)=k/2*(vr(4)+r(end-1)*c21{1}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)))-r(end-1)*rho*Cp/dt;
% Deuxième rangée
 % en z
a(N^2-2*N+2,(N^2-2*N+1):(N^2-2*N+5))=r(end-1)*k/2*[c22{1}(z(1)-z(2),z(3)-z(2),z(4)-z(1),z(5)-z(2))...
    ,c22{2}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)),c22{3}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2))...
    ,c22{4}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)),c22{5}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2))];
 % en r
a(N^2-2*N+2,(N^2-2*N-3*N+2):N:(N^2-N+2))=vr;
% point milieu
a(N^2-2*N+2,N^2-2*N+2)=k/2*(vr(4)+r(end-1)*c22{2}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)))-r(end-1)*rho*Cp/dt;
%Centre
for i=3:N-2
     % en z
  a(N^2-2*N+i,(N^2-2*N+i-2):(N^2-2*N+i+2))=r(end-1)*k/2*[c23{1}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))...
    ,c23{2}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)),c23{3}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))...
    ,c23{4}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)),c23{5}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))];
     % en r
  a(N^2-2*N+i,(N^2-2*N-3*N+i):N:(N^2-2*N+N+i))=k/2*vr;
  % point milieu
  a(N^2-2*N+i,N^2-2*N+i)=k/2*(vr(4)+r(end-1)*c23{3}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)))-r(end-1)*rho*Cp/dt;
end
%Avant dernière rangée
 % en z
a(N^2-N-1,(N^2-N-4):(N^2-N))=r(end-1)*k/2*[c24{1}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))...
    ,c24{2}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)),c24{3}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))...
    ,c24{4}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)),c24{5}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))];
 % en r
a(N^2-N-1,(N^2-4*N-1):N:(N^2-N+N-1))=k/2*vr;
% point milieu
a(N^2-N-1,N^2-N-1)=k/2*(vr(4)+r(end-1)*c24{4}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)))-r(end-1)*rho*Cp/dt;
%Dernière rangée
 % en z
a(N^2-N,(N^2-N-4):(N^2-N))=r(end-1)*k/2*[c25{1}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))...
    ,c25{2}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)),c25{3}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))...
    ,c25{4}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)),c25{5}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))];
 % en r
a(N^2-N,(N^2-N-3*N):N:(N^2-N+N))=k/2*vr;
% point milieu
a(N^2-N,N^2-N)=k/2*(vr(4)+r(end-1)*c25{5}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)))-r(end-1)*rho*Cp/dt;

%% Dernière colonne
 % Terme général de r
 ap=r(end-4)-r(end);
 bp=r(end-3)-r(end);
 cp=r(end-2)-r(end);
 dp=r(end-1)-r(end);
 vr=k/2*[1*c15{1}(ap,bp,cp,dp)+r(end)*c25{1}(ap,bp,cp,dp)...
    1*c15{2}(ap,bp,cp,dp)+r(end)*c25{2}(ap,bp,cp,dp)...
    1*c15{3}(ap,bp,cp,dp)+r(end)*c25{3}(ap,bp,cp,dp)...
    1*c15{4}(ap,bp,cp,dp)+r(end)*c25{4}(ap,bp,cp,dp)...
    1*c15{5}(ap,bp,cp,dp)+r(end)*c25{5}(ap,bp,cp,dp)];

% Première rangée
 % en z
a(N^2-N+1,(N^2-N+1):(N^2-N+5))=r(end)*k/2*[c21{1}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))...
    ,c21{2}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)),c21{3}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))...
    ,c21{4}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)),c21{5}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1))];
 % en r
a(N^2-N+1,(N^2-5*N+1):N:(N^2-N+1))=k/2*vr;
 % point milieu
a(N^2-N+1,N^2-N+1)=k/2*(vr(5)+r(end)*c21{1}(z(2)-z(1),z(3)-z(1),z(4)-z(1),z(5)-z(1)))-r(end)*rho*Cp/dt;
% Deuxième rangée
 % en z
a(N^2-N+2,(N^2-N+1):(N^2-N+5))=k/2*r(end)*[c22{1}(z(1)-z(2),z(3)-z(2),z(4)-z(1),z(5)-z(2))...
    ,c22{2}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)),c22{3}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2))...
    ,c22{4}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)),c22{5}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2))];
 % en r
a(N^2-N+2,(N^2-2*N-3*N+2):N:(N^2-N+2))=k/2*vr;
 % point milieu
a(N^2-N+2,N^2-N+2)=k/2*(vr(5)+r(end)*c22{2}(z(1)-z(2),z(3)-z(2),z(4)-z(2),z(5)-z(2)))-r(end)*rho*Cp/dt;
%Centre
for i=3:N-2
     % en z
  a(N^2-N+i,(N^2-N+i-2):(N^2-N+i+2))=r(end)*k/2*[c23{1}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))...
    ,c23{2}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)),c23{3}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))...
    ,c23{4}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)),c23{5}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i))];
     % en r
  a(N^2-N+i,(N^2-N-4*N+i):N:(N^2-2*N+N+i))=k/2*vr;
     % point milieu
  a(N^2-N+i,N^2-N+i)=k/2*(vr(5)+r(end)*c23{3}(z(i-2)-z(i),z(i-1)-z(i),z(i+1)-z(1),z(i+2)-z(i)))-r(end)*rho*Cp/dt;
end
%Avant dernière rangée
 % en z
a(N^2-1,(N^2-4):(N^2))=r(end)*k/2*[c25{1}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))...
    ,c25{2}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)),c25{3}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))...
    ,c25{4}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)),c25{5}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1))];
 % en r
a(N^2-1,(N^2-4*N-1):N:(N^2-1))=k/2*vr;
 % point milieu
a(N^2-1,N^2-1)=k/2*(vr(5)+r(end)*c25{4}(z(end)-z(end-1),z(end-2)-z(end-1),z(end-3)-z(end-1),z(end-4)-z(end-1)))-r(end)*rho*Cp/dt;
%Dernière rangée
 % en z
a(N^2,(N^2-4):(N^2))=r(end)*k/2*[c25{1}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))...
    ,c25{2}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)),c25{3}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))...
    ,c25{4}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)),c25{5}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end))];
 % en r
a(N^2,(N^2-4*N):N:(N^2))=k/2*vr;
 % point milieu
a(N^2,N^2)=k/2*(vr(5)+r(end)*c25{5}(z(end-1)-z(end),z(end-2)-z(end),z(end-3)-z(end),z(end-4)-z(end)))-r(end)*rho*Cp/dt;
toc
