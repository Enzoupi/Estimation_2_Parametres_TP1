function [x0,Jx,GJx,nit] = GOPT(J,GJ,x0,epsil,nitmax,findic)
% Algorithme de descente du gradient à pas optimal avec intégration du
% calcul optimal par des formules analytiques pour les fonctions étudiées
%% Paramètres et déclarations utiles
error = 10000;
nit = 0;
i = 1:length(x0);
% solution exacte pour calcul de l'erreur
if findic == 1 || findic == 2
    solex = i;
end

%% On itère le processus jusqu'à se rapprocher de la solution
while (error > epsil) && (nit < nitmax)
    dk = -GJ(x0,findic);
    if findic == 1
        dk = -GJ(x0,findic);
        b = 2*sum( (x0-i) .* dk);
        c = sum(dk .* dk);
        pas = -b/(2*c);
    elseif findic == 2
        b = 2*sum( 1./i .* (x0-i) .* dk);
        c = sum(1./i .* dk .* dk);
        pas = -b/(2*c);
    else
        disp('---------------------------------------------------------------')
        disp('La fonction demandée nest pas encore implémentée !!')
        disp('---------------------------------------------------------------')
        return
    end
    x0 = x0 + pas .* dk;
    error = max(abs(x0-solex));
    nit = nit +1;
    Jx = J(x0,findic);
    GJx = GJ(x0,findic);
end

end

