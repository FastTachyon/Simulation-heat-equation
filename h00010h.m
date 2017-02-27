function [d,d2]=h00010h()
% fonction retrounant les formules de différence finies comme function handle
% Ne prend aucune entrée
% le 1 dans le titre représente où les formules de dérivations sont
% centrées.
% d est la formule de différence finie pour la première dérivée dans une
% grille non-uniforme
% d2 est la formule de différence fini pour la dérivée seconde dans un
% grille non uniforme.
% Auteur: Mathieu Bergeron, Oscar Cespedes, Blandine Rippert
% Fait le 18 février 2017
%% Formule de dérivation [0 0 0 1 0]
syms a b c d y1 y2 y3 y4 y5
l=[1 a a^2/2 a^3/6 a^4/24;1 b b^2/2 b^3/6 b^4/24; ...
    1 c c^2/2 c^3/6 c^4/24;1 0 0 0 0 ;...
    1 d d^2/2 d^3/6 d^4/24];
l=inv(l);
d=matlabFunction(l(2,1)*y1+l(2,2)*y2+l(2,3)*y3+l(2,4)*y4+l(2,5)*y5);
d2=matlabFunction(l(3,1)*y1+l(3,2)*y2+l(3,3)*y3+l(3,4)*y4+l(3,5)*y5);