function [u,lambda,i,resi_l,resi]=UZAWA(A,C,g,f,Lambda0,rho,eta,Eps,Imax)

lambda1=Lambda0;
i=0;
resi=[]; % vecteur des ||u_k-u_{k-1}||
resi_l=[]; % vecteur des ||lambda_k-lambda_{k-1}||

%initialisation
resi_u=1;
N=size(g,1)/2;
x0=zeros(2*N,1);

b=-g-C'*lambda1;
[u1,~,~]=GC(A,b,x0,eta,Imax);



while (i<Imax && norm(max(0,C*u1-f))>eta || resi_u>Eps)
  i=i+1;
  u=u1; % stock u_{k-1}
  lambda=lambda1; % stock lambda_{k-1}
  
  % calcul de u_k et lambda_k
  lambda1=max(lambda+rho*(C*u-f),0);
  b=-g-C'*lambda1;
  [u1,~,~]=GC(A,b,x0,eta,Imax);
  nr=norm(A*u1-b);
  resi_u = norm(u1-u);
  
  %affichage des résidus toutes les 10 itérations
    if(mod(i,10) ==0)
      fprintf('||z_k - z_(k-1)||=%e, grad(J).z_k=%d\n',norm(u1-u),nr);
    end
    
  % On insère le nouveau résidu dans notre tableau pour le tracking ; 
  n1=norm(lambda1-lambda);
  resi=[resi;resi_u]; 
  resi_l=[resi_l;n1];

  % Si le résidu devient trop grand, sans doute il est temps d'en rester
  % là ...
  if nr>1e10; fprintf('  explosion !\n'); break;
  end
end

end

