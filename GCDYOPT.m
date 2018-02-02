function [x,Jx,GJx,nit] = GCDYOPT(J,GJ,x0,epsil,nitmax,findic)
%% Paramètres divers
error = 10000;      % Initialisation d'une valeur élévée pour l'erreur
nit = 0;            % Nombre d'itérations
i = 1:length(x0);   % Commodité pour réaliser les calculs

%% Première étape : gradient à pas optimal normal
% Initialisation des premiers calculs de pas optimaux
if findic == 1
    b = 2 * sum( (x0 -i).*GJ(x0,findic) );
    c = sum( GJ(x0,findic).*GJ(x0,findic) );
    pas = -b/(2*c);
    solex = i;
elseif  findic == 2
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
% Calcul du premier pas par l'algo du gradient à pas optimal
x = x0 + pas .* GJ(x0,findic);
error = max(abs(x-solex));
nit = nit +1;
dkm1 = -GJ(x0,findic);
xm1 = x0;
x0=x;

%% Gradient conjugué de Dai Yuan
while (error > epsil) && (nit < nitmax)
    Bkm1 = sum(GJ(x,findic).*GJ(x,findic))/(sum(dkm1.*(GJ(x,findic)-GJ(xm1,findic))));
    dk = -GJ(x,findic) + Bkm1.*dkm1;
    xm1 = x;
    % Calcul du pas optimal
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
    x = x + pas .* dk;
    dkm1 = dk;
    error = abs(norm(GJ(x,findic)));
    nit = nit +1;
end
Jx = J(x,findic);
GJx = GJ(x,findic);

