function [a,b,c,d,e,f,g,h,i,j]=T01000T()
% fonction retrounant les formules de diff�rence finies comme function handle
% Ne prend aucune entr�e
% le 1 dans le titre repr�sente o� les formules de d�rivations sont
% centr�es.
% a � e sont les 5 param�tres de la d�riv�e premi�re
% f � j sont les 5 param�tres matriciels de la d�riv�e seconde
% Auteur: Mathieu Bergeron, Oscar Cespedes, Blandine Rippert
% Fait le 18 f�vrier 2017

%% Formule de d�rivation [0 1 0 0 0]
syms a b c d y1 y2 y3 y4 y5
l=[1 a a^2/2 a^3/6 a^4/24; 1 0 0 0 0; ...
    1 b b^2/2 b^3/6 b^4/24; 1 c c^2/2 c^3/6 c^4/24;...
    1 d d^2/2 d^3/6 d^4/24];
l=inv(l);
a=cell(1,5);
b=cell(1,5);
a{1}=matlabFunction(l(2,1));
a{2}=matlabFunction(l(2,2));
a{3}=matlabFunction(l(2,3));
a{4}=matlabFunction(l(2,4));
a{5}=matlabFunction(l(2,5));
b{1}=matlabFunction(l(3,1));
b{2}=matlabFunction(l(3,2));
b{3}=matlabFunction(l(3,3));
b{4}=matlabFunction(l(3,4));
b{5}=matlabFunction(l(3,5));