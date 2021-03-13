function [Z,i,resi]=GP(A,b,f,Z0,eta,rho,Imax)

Z1=Z0;
i=0;
resi=[];
N=size(Z1,1)/2;

% On insère le nouveau résidu dans notre tableau pour le trackinb ; 
r=b-A*Z1;
nr=norm(r);
e=1;

while (i<Imax && e>eta)
  i=i+1;
  % Algorithme et mise à jour : à remplir
  Z=Z1;
  Z1=Z1+rho*r;
  Z1(N+1:2*N)=max(Z1(N+1:2*N),f);
  r=b-A*Z1;
  nr=norm(r);
  e = norm(Z-Z1);

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

