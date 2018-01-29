clear all
close all
%% Zone de commentaire
% Programme principal pour le TP1 d'estimation de paramètres

% remarque : question 1 regarder aussi la convergence et non pas que le
% nombre d'iterations

%Remarque : fonction de Rosenvrock : le gradient à pas constant va
%fonctionner pour des pas assez petits

%% Cas test avec un vecteur 
pas = 0.05;
epsil = 0.01;
nitmax = 1000;
findic = 1;
x0 = [2,2,2,2,5];
[x,Jx,GJx,nit] = GCST(@J,@GJ,x0,pas,epsil,nitmax,findic);

%% Question 2
figure
hold on
for ndim = 1
    % Solution exacte
    sol = 1:length(ndim);
    findic=1;
    x0 = zeros(1,ndim);
    epsil = 1e-7;
    pas = 0.1:0.1 :5;
    nitmax = 5000;
    nitconv=zeros(1,length(pas));
    res=zeros(1,length(pas));
    for i = 1:length(pas)
        [res(i),~,~,nitconv(i)]=GCST(@J,@GJ,x0,pas(i),epsil,nitmax,findic);
    end
    plot(pas,nitconv)
end
title('Convergence en fonction du pas')
xlabel('Taille du pas')
ylabel('Nombre diterations avant convergence')
legend('n=1')
hold off

figure
loglog(pas,abs(sol-res)/sol)
title('Ecart relatif en fonction du pas')
xlabel('Taille du pas')
ylabel('Ecart relatif')

%% Cas Test avec Dai Yuan
pas = 0.05;
epsil = 1e-7;
nitmax = 1000;
findic = 1;
x0 = [4,4,4];
[x,Jx,GJx,nit] = GCDYCST(@J,@GJ,x0,pas,epsil,nitmax,findic);

%% Question 3 : Convergence avec Dai Yuan
findic = 1;
pas = 0.1:0.1 :5;
err = zeros(2,length(pas));
res = zeros(1,length(pas));
nitmax = 5000;
figure
hold on
for ndim = 1:2
    solex = 1:ndim;
    findic=1;
    x0 = zeros(1,ndim);
    epsil = 1e-7;
    nitconv=zeros(1,length(pas));
    for i = 1:length(pas)
        [res,~,~,nitconv(i)]=GCDYCST(@J,@GJ,x0,pas(i),epsil,nitmax,findic);
    end
    err(ndim,i) = max(abs(solex-res)/solex);
    plot(pas,nitconv)
end
title('Convergence Pour GC de Dai Yuan')
xlabel('Taille du pas')
ylabel('Nombre diterations avant convergence')
legend('n=1','n=2')
hold off

%% Question 4 :
findic = 2;
x0 = [0,0];
solex = 1:length(x0);
pas = 0.1:0.1 :5;
errgcdy = zeros(1,length(pas));
errgcst = zeros(1,length(pas));
for i=1:length(pas)
    [resgcdy,~,~,nitgcdy(i)]=GCDYCST(@J,@GJ,x0,pas(i),epsil,nitmax,findic);
    [resgcst,~,~,nitg(i)]=GCST(@J,@GJ,x0,pas(i),epsil,nitmax,findic);
    errgcdy(i)=max(abs(solex-resgcdy)/solex);
    errgcst(i)=max(abs(solex-resgcst)/solex); 
end
figure
hold on
plot(pas,nitgcdy)
plot(pas,nitg)
title('Comparaison des convergences pour deux versions de la descente du gradient')
xlabel('Pas')
ylabel('Nombre d iterations avant convergence')
legend('Gradient conjugué (DY)','Gradient à pas constant')
figure
hold on
plot(pas,errgcdy)
plot(pas,errgcst)
ylim([0 1e-7])
title('Ecart relatif maximum')
xlabel('Pas')
ylabel('Ecart relatif')
legend('Gradient conjugué (DY)','Gradient à pas constant')

%% Question 5
findic = 3; % On choisit la fonction oscillante
x0 = 0:0.01:4;
pas = 0.01;
epsil = 1e-7;
nitmax = 5000;
% Détermination des minimas
res = zeros(1,length(x0));
for i=1:length(x0)
    [res(i),~,~,nitg(i)]=GCST(@J,@GJ,x0(i),pas,epsil,nitmax,findic);
end

figure
hold on
plot(x0,J(x0,findic))
plot(x0,res)
title('Minimum obtenu en fonction du point initial')
xlabel('x')
ylabel('y')
legend('Fonction oscillante','Minima obtenu')