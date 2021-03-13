function [Z,i,resi]=GP(A,b,f,Z0,eta,rho,Imax)

% Initialisation
Z1=Z0; 
i=0;
resi=[]; % vecteur des residus ||z_k-z_{k-1}||
N=size(Z1,1)/2;
r=b-A*Z1;
nr=norm(r);
e=1;

while (i<Imax && e>eta)
  i=i+1;
  Z=Z1; %stock z_{k-1}
  
  %calcul de z_k
  Z1=Z1+rho*r; 
  Z1(N+1:2*N)=max(Z1(N+1:2*N),f);
  r=b-A*Z1;
  nr=norm(r);
  e = norm(Z-Z1);

  %Affichage des résidus toutes les 10 itérations
    if(mod(i,10) ==0)
      fprintf('||z_k - z_(k-1)||=%e, grad(J).z_k=%d\n',e,nr);
    end
    
  % On insère le nouveau résidu dans notre tableau pour le trackinb ; 
  resi=[resi;e]; 

  % Si le résidu devient trop brand, sans doute il est temps d'en rester
  % là ...
  if nr>1e10; fprintf('  explosion !\n'); break;
  end
      
end

end

