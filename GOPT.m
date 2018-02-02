function [x,Jx,GJx,nit] = GOPT(J,GJ,x0,epsil,nitmax,findic)
% Algorithme de descente du gradient à pas optimal avec intégration du
% calcul optimal par des formules analytiques pour les fonctions étudiées
%% Paramètres et déclarations utiles
error = 10000;
nit = 0;
i = 1:length(x0);
%% Initialisation des premiers calculs de pas optimaux
if findic == 1
    b = 2 * sum( (x0 -i).*GJ(x0,findic) );
    c = sum( GJ(x0,findic).*GJ(x0,findic) );
    pas = -b/(2*c);
    solex = i;
elseif findic == 2
    
    b = 2 * sum( 1./i.*(x0-i).*GJ(x0,findic) );
    c = sum(1./i.*GJ(x0,findic).*GJ(x0,findic));
    pas = -b/(2*c);
    solex = i;
else
    disp('---------------------------------------------------------------')
    disp('La fonction demandée nest pas encore implémentée !!')
    disp('---------------------------------------------------------------')
    return 
end
%% On ré-itère le processus jusqu'à se rapprocher de la solution
while (error > epsil) && (nit < nitmax)
    x = x0 + pas .* GJ(x0,findic);
    error = max(abs(x-solex));
    nit = nit +1;
    x0=x;
    if findic == 1
        b = 2*sum((x0-(1:length(x0))).*GJ(x0,findic));
        c = sum(GJ(x0,findic).*GJ(x0,findic));
        pas = -b/(2*c);
    elseif findic == 2
        i = 1:length(x0);
        b = 2*sum(1./i.*(x0-(1:length(x0))).*GJ(x0,findic));
        c = sum(1./i.*GJ(x0,findic).*GJ(x0,findic));
        pas = -b/(2*c);
    end
    Jx = J(x,findic);
    GJx = GJ(x,findic);
end

end

