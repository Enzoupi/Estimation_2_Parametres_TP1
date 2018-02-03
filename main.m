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

%% Exercice 2 : Gradient à pas optimal
disp('---------------------------------------------------------------')
disp('Exercice 2')
disp('---------------------------------------------------------------')
% ==> Question 1 :

% Validation des descente de gradient à pas optimal : on sait que pour la
% fonction f on converge en une itération exactement et on retrouve ce
% comportement :
epsil = 1e-7;
nitmax = 1000;
findic = 1;
x0 = [5 3 64 15 -135];
disp('Avec le Gradient à pas optimal')
[x,Jx,GJx,nit,nbeval] = GOPT(@J,@GJ,x0,epsil,nitmax,findic);
disp(['Point de départ au hasard    : ' num2str(x0)])
disp(['Minimum obtenu               : ' num2str(x)])
disp(['En exactement ' num2str(nit) ' itération(s)'])
disp('Avec le Gradient Conjugué de Dai Yuan à pas optimal')
[x,Jx,GJx,nit,nbeval] = GCDYOPT(@J,@GJ,x0,epsil,nitmax,findic);
disp(['Point de départ au hasard    : ' num2str(x0)])
disp(['Minimum obtenu               : ' num2str(x)])
disp(['En exactement ' num2str(nit) ' itération(s)'])

% Comparaison avec GCST pour n=2 et pas = 0.5
epsil = 1e-7;
nitmax = 1000;
findic = 2;
x0 = [15 -15];
pas = 0.5;
[resgopt,Jx,GJx,nit_gopt,nbeval] = GOPT(@J,@GJ,x0,epsil,nitmax,findic);
[resgcst,~,~,nit_gcst]=GCST(@J,@GJ,x0,pas,epsil,nitmax,findic);
[resgcdyopt,~,~,nit_gcdyopt,nbeval] = GCDYOPT(@J,@GJ,x0,epsil,nitmax,findic);
disp(['Comparaison sur la fonction g pour n = 2 avec une tolerance de ' num2str(epsil)])
disp(['Nombre d itérations pour GCST : ' num2str(nit_gcst)])
disp(['Nombre d itérations pour GOPT : ' num2str(nit_gopt)])
disp(['Nombre d itérations pour GCDYOPT : ' num2str(nit_gcdyopt)])


% ==> Question 2 :
disp('Exercice 2 Question 2:')
% Pour f
nn = [1,2,5,10,100,1000,10000];
epsil = 1e-7;
findic =2;
err_gopt = zeros(size(nn));
err_gcdyopt = zeros(size(nn));
nit_gopt = zeros(size(nn));
nit_gcdyopt = zeros(size(nn));
nbeval_gopt = zeros(size(nn));
nbeval_gcdyopt = zeros(size(nn));
for i=1:length(nn)
    n=nn(i);
    x0 = zeros(1,n);
    solex = 1:n;
    [resgopt,~,~,nit_gopt(i),nbeval_gopt(i)]=GOPT(@J,@GJ,x0,epsil,nitmax,findic);
    [resgcdyopt,~,~,nit_gcdyopt(i),nbeval_gcdyopt(i)] = GCDYOPT(@J,@GJ,x0,epsil,nitmax,findic);
    err_gopt(i) = max(abs(resgopt-solex));
    err_gcdyopt(i) = max(abs(resgcdyopt-solex));
end

% 
% figure
% hold on
% subplot(2,1,1)
% plot(nn,log10(err_gcdyopt))
% plot(nn,log10(err_gopt))
% legend('Gradient Conjugué DY Optimal','Gradient Optimal')
% subplot(2,1,2)
% plot(nn,nit_gcdyopt)
% plot(nn,nit_gopt)
% legend('Gradient Conjugué DY Optimal','Gradient Optimal')

figure
hold on
plot(log10(nn),log10(nbeval_gopt))
plot(log10(nn),log10(nbeval_gcdyopt))
legend('Gradient Optimal','Gradient Conjugué DY Optimal')
xlabel('log10(n)')
ylabel('log10 nombre d évaluation de J et GJ')

% représentation des trajets
epsil = 1e-7;
nitmax = 1000;
findic = 2;
x0 = [5 -5];
pas = 0.5;
[resgopt,Jx,GJx,nit_gopt,nbeval] = GOPT_allx(@J,@GJ,x0,epsil,nitmax,findic);
[resgcdyopt,Jx,GJx,nit,nbeval] = GCDYOPT_allx(@J,@GJ,x0,epsil,nitmax,findic);

figure
subplot(2,1,1)
hold on
plot(resgopt(:,1),resgopt(:,2))
plot(resgcdyopt(:,1),resgcdyopt(:,2))
axis([-1 6 -6 6])
plot(1,2,'ro')
text(1.1,2.1,'Minimum de la fonction en [1 2]')
xlabel('x')
ylabel('y')
legend('Gradient optimal','Gradient Conjugué DY Optimal')

subplot(2,1,2)
hold on
plot(resgopt(:,1),resgopt(:,2))
plot(resgcdyopt(:,1),resgcdyopt(:,2))
axis([0.99 1.01 1.99 2.01])
plot(1,2,'ro')
text(1.1,2.1,'Minimum de la fonction en [1 2]')
xlabel('x')
ylabel('y')
legend('Gradient optimal','Gradient Conjugué DY Optimal')
