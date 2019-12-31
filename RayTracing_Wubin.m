clear all;

%% raytracing trough arbitrary lenses program by Wubin Pang
% parameters setting up
alpha=10;
kappa=1/2500;
delta=4;
sigma=1;
N=1.3;   %% INDEX of glass
D=20;    %% aperture size of each surface
sig=1e-3;    
num_surf=6;  %% number of total surfaces
nu_ray=20;   %% number of rays 

% Build ray matrix
ray_x=linspace(-D/2+1,D/2-1,nu_ray);
ray_z=zeros(1,nu_ray);
ray_vecx=zeros(1,nu_ray);
ray_vecz=ones(1,nu_ray);
ray_mat=[ray_x',ray_z',ray_vecx',ray_vecz']; 


for j=1:nu_ray  %% iterate through each ray from various positions
  
      x0=ray_mat(j,1);  %% initial x coordinate of ray
      z0=ray_mat(j,2);  %% initial z coordinate of ray
      vecx0=ray_mat(j,3);  %% optical direction cosine of x
      vecz0=ray_mat(j,4);  %% optical direction consine of z 
     
    
  for i=2:num_surf+1  %% iterate through every sequential surface
      
     if rem(i,2)==0  %% even surface number
        
      if abs(vecx0)<sig;
          x=x0;
          z=i*delta+kappa*(x-alpha)*(x+alpha)*x^2;          
          
      else
          syms x_s;
          f=i*delta+kappa*(x_s-alpha)*(x_s+alpha)*x_s^2-(x_s-x0)*vecz0/vecx0-z0;   %% solve for intersection point between ray and surface
          x_roots=solve(f);
          roots=double(x_roots);
          if abs(roots(1))>D/2&&abs(roots(2))>D/2
              continue;
          elseif abs(roots(1))<=D/2;
              x=roots(1);        
          elseif abs(roots(2))<=D/2;
              x=roots(2);
          end
          
           z=(x-x0)*vecz0/vecx0+z0;
          
      end
      n_x=-4*kappa*x^3+2*kappa*alpha^2*x;        %% solve for emergent ray direction vector
      n_z=1;
      nor_x=n_x/sqrt(n_x^2+n_z^2);
      nor_z=n_z/sqrt(n_x^2+n_z^2);
      P=sqrt(N^2-1+(nor_x*vecx0+nor_z*vecz0)^2)-nor_x*vecx0+nor_z*vecz0;
      vecx=vecx0+P*nor_x;
      vecz=vecz0+P*nor_z;      
      
      plot([z0,z],[x0,x],'-r');
      hold on;
      
      vecx0=vecx;
      vecz0=vecz;
      x0=x;
      z0=z;
      
     else                              %% for odd surface number
          if abs(vecx0)<sig;          
          x=x0;
          z=(i-1)*delta+sigma-kappa*(x-alpha)*(x+alpha)*x^2;     
          
          else
          syms x_s;
          f=(i-1)*delta+sigma-kappa*(x_s-alpha)*(x_s+alpha)*x_s^2-(x_s-x0)*vecz0/vecx0-z0;
          x_roots=solve(f);
          roots=double(x_roots);
          if abs(roots(1))>D/2&&abs(roots(2))>D/2
              continue;
          elseif abs(roots(1))<=D/2;
              x=roots(1);          
         
          elseif abs(roots(2))<=D/2;
              x=roots(2);
             
          end  
         z=(x-x0)*vecz0/vecx0+z0;
          
          end
      n_x=4*kappa*x^3-2*kappa*alpha^2*x;
      n_z=1;
      nor_x=n_x/sqrt(n_x^2+n_z^2);
      nor_z=n_z/sqrt(n_x^2+n_z^2);
      P=sqrt(1-N^2+N^2*(nor_x*vecx0+nor_z*vecz0)^2)-N*(nor_x*vecx0+nor_z*vecz0);
      vecx=vecx0+P*nor_x;
      vecz=vecz0+P*nor_z;
          
     plot([z0,z],[x0,x],'-r');
      hold on;
      
      vecx0=vecx;
      vecz0=vecz;
      x0=x;
      z0=z;
         
     end      
       
  end 
  
  z=50;
  x=(z-z0)*vecx0/vecz0+x0;
  plot([z0,z],[x0,x],'-r');
  hold on;
   
end

xx=linspace(-D/2,D/2,50);                 %% draw the surfaces
for n=2:num_surf+1
    if rem(n,2)==0
      zz=n*delta+kappa*(xx-alpha).*(xx+alpha).*xx.^2;
      plot(zz,xx,'b');
      hold on;  
    else
        zz=(n-1)*delta+sigma-kappa*(xx-alpha).*(xx+alpha).*xx.^2;
        plot(zz,xx,'k');
        hold on;  
    end
  
end




