function [a,b]=T00100T()
% fonction retrounant les formules de différence finies comme function handle
% le 1 dans le titre représente où les formules de dérivations sont
% centrées.
% Ne prend aucune entrée
% a à e sont les 5 paramètres de la dérivée première
% f à j sont les 5 paramètres matriciels de la dérivée seconde
% grille non uniforme.
% Auteur: Mathieu Bergeron, Oscar Cespedes, Blandine Ripert
% Fait le 18 février 2017

%% Formule de dérivation centrée à 5 points [0 0 1 0 0]
syms a b c d 
l=[1 a a^2/2 a^3/6 a^4/24;1 b b^2/2 b^3/6 b^4/24; ...
    1 0 0 0 0; 1 c c^2/2 c^3/6 c^4/24;...
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