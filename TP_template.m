mode(-1);
clear all;
close ();

%! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%! Methode GC: Gradient a pas fixe projete
%! Methode UZAWA: Methode d'UZAWA
%! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
global PRINT

% EXERCICE='Exo1';
% EXERCICE='Exo2';
% EXERCICE='Exo2b';
 EXERCICE='Exo3';

figure (1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%! Definition des matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=500;          % taille du probleme
fprintf('N=%i\n',N);

Mtot=40.0;             %masse totale de la chaîne
mi=Mtot/(N+2);         %masse de chaque point
omega=70*(N+1);        % coefficient de ressort
g0=9.81;

% coordonnées des extrémités
x0=-2;
y0=1.0;
xf=2.0;
yf=1.0;


J=diag(ones(N-1,1),1)*omega;

A=(2*eye(N,N)*omega-J-J');
ON=zeros(N,N);

%matrice H 
H=[A ON  ; ON A];

% le vecteur g  : compléter

 g=zeros(2*N,1);
 g(1)=omega*x0;
 g(N)=omega*xf;
 for i=N+1:2*N
     g(i)=-g0*mi;
 end
 g(N+1)=g(N+1)+omega*y0;
 g(2*N)=g(2*N)+omega*yf;
 

% valeurs propres min et max de la matrice H : 2N x 2N

L1= 4*omega*sin(pi/(2*(2*N+1)))^2;
L2N= 4*omega*sin(N*pi/(2*N+1))^2;

eta=1.e-4;      % residu desire'
Imax=1000;	% nombre d'iterations maximal



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%! Reponses par exercice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch EXERCICE
    
    case 'Exo1'
        
        fprintf('RESOLUTION PAR Gradient conjuge:\n');
        
        
        zInit=ones(2*N,1);  % Configuation initiale pour GPF
        %timer();
        tic;
        [u,n,Residu_u]=GC(H,g, zInit,eta, Imax);
        t=toc;
        err=Residu_u(end);
        fprintf('temps=%5.2f, Iterations:%5i, ||U^k-U^{k-1}||=%10.2e\n',t,n,Residu_u(end));
        
        xx=[x0; u(1:N); xf];
        yy=[y0; u(N+1:2*N); yf];
        
        hold on;
        plot(xx,yy,'r', 'Linewidth', 2);
        hold on
        
        
    case 'Exo2'
        fprintf('RESOLUTION PAR Gradient projete:\n');
        rho=2/(L1+L2N); % pas pour la methode GPF
        zini=ones(2*N,1);  % Configuation initiale pour GPF
        % Calcul du rayon r=|| I - \rho A || = (L_N-L_1)/(L_N+L_1)
         
        r=(L2N-L1)/(L2N+L1); fprintf('r  :=%f\n',r);
        
        % vecteur f : compléter 
        f= 0.5*ones(N,1);
        
        eta2=(1-r)/r*eta;
        fprintf('eta=%5.2e; eta2=%5.2e\n',eta,eta2);
        timer();
        tic;
        [u,n,Residu_u]=GP(H,g,f,zini,eta2,rho,10*Imax);
        t=toc;
        fprintf('temps=%5.2f, Iterations:%5i, ||U^k-U^{k-1}||=%10.2e\n',t,n,Residu_u(end));
        
        err=Residu_u(end);
        
        
        xx=[x0; u(1:N); xf];
        yy=[y0; u(N+1:2*N); yf];
        
       
        plot(xx,yy,'r', 'Linewidth', 2);
        hold on;
        legend('Solution');
        hold on
        plot(xx,  0.5*ones(N+2,1), 'k', 'Linewidth', 2);
        legend('Solution','Obstacle');
       
    case 'Exo2b'
        t_N=10:10:1000;
        log_N=log(t_N);
        L1= 4*omega*sin(pi./(2*(2*t_N+1))).^2;
        L2N= 4*omega*sin(t_N*pi./(2*t_N+1)).^2;
        
        
        gamma=(L2N-L1)./(L2N+L1);
        eta=(1-gamma)./gamma * 10^(-4);
        
%        semilogx(t_N,gamma,'k','Linewidth',2);
%        xlabel('N');
%        ylabel('\gamma');
        
       semilogx(t_N,eta,'k','Linewidth',2);
       xlabel('N');
       ylabel('\eta');  
        
    case 'Exo3'
        
        % Définir la matrice C, vecteur f
%         f= 0.75:0.05/(N-1):0.8;
%         f=-f';
        f= -0.5*ones(N,1);
       
        C=[ON, -eye(N)];
        
        rho=L1;
        gamma=(L2N-L1)./(L2N+L1);
        eta=(1-gamma)./gamma * eta;
  
        Eps=1e-5;
        Lambda0=zeros(size(f));
        fprintf('RESOLUTION PAR UZAWA:\n');
        tic;
        [u,Lambda,n,Residu_Lambda,Residu_u] = UZAWA(H,C,-g,f,Lambda0,rho,eta,Eps,Imax);
        t=toc;

        %  compléter pour les graphiques 
        xx=[x0; u(1:N); xf];
        yy=[y0; u(N+1:2*N); yf];
        plot(xx,yy,'r', 'Linewidth', 2);
        hold on;
        legend('Solution');
        hold on
        fprintf('temps=%5.2f, Iterations:%5i, ||U^k-U^{k-1}||=%10.2e\n',t,n,Residu_u(end));
        
        plot(xx,  [0.5;-f;0.5], 'k', 'Linewidth', 2);
        legend('Solution','Obstacle');
        
end



