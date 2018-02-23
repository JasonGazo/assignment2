function [Jsum] = getVmap( L,width, sigvar )
%Function performs just like Assignment2P2 script except it takes in
%input parameters, outputs total current, and doesn't plot anything
clc
%Fixed Voltage
vx=1;

%Set Boundaries (L is input argument)
nx = 2*L/3;
ny = L;

%Make Sure bottleneck is less than or equal to frame length

if(width>ny)
    error('Width is too large')
end

%Create G Matrix
G = sparse(nx*ny);
%Boundaries
v = zeros(1,nx*ny);

%Conductivities in Frame (Sigvar is input argument for function)
sig1 = 1;
sig2 = sigvar;

%Establish dimensions of boxes
%[Left side, right side, top, bottom]
box1 = [nx*2/5 nx*3/5 ny (ny/2 + width/2)];
box2 = [nx*2/5 nx*3/5 (ny/2-width/2) 0];

%Conductivity Matrix
sigma=ones(nx,ny);
for i=1:nx
    for j=1:ny
        
        if(i > box1(1) && i < box1(2) && (j < box2(3)||j > box1(4)))
            sigma(i,j)=sigvar;
        end
        
    end
end

for i=1:nx
    
   for j = 1:ny
       n=j+(i-1)*ny;  
       
       if(i==1)
           G(n,:)=0;
           G(n,n)=1;
           v(n)=vx;
           
       elseif(i==nx)
           G(n,:)=0;
           G(n,n)=1;
           v(n)=0;
           
       elseif (j == 1)
           
            if (i > box1(1) && i < box1(2))
                G(n, n) = -3;
                G(n, n+1) = sig2;
                G(n, n+ny) = sig2;
                G(n, n-ny) = sig2;
                
            else
                G(n, n) = -3;
                G(n, n+1) = sig1;
                G(n, n+ny) = sig1;
                G(n, n-ny) = sig1;
                
            end
            
       elseif (j == ny)
           
           if (i > box1(1) && i < box1(2))
                G(n, n) = -3;
                G(n, n-1) = sig2;
                G(n, n+ny) = sig2;
                G(n, n-ny) = sig2;
                
            else
                G(n, n) = -3;
                G(n, n-1) = sig1;
                G(n, n+ny) = sig1;
                G(n, n-ny) = sig1;
                
           end
            
       else
            
            if (i > box1(1) && i < box1(2) && (j < box2(3)|| j > box1(4)))
                G(n, n) = -4;
                G(n, n+1) = sig2;
                G(n, n-1) = sig2;
                G(n, n+ny) = sig2;
                G(n, n-ny) = sig2;
                
            else
                G(n, n) = -4;
                G(n, n+1) = sig1;
                G(n, n-1) = sig1;
                G(n, n+ny) = sig1;
                G(n, n-ny) = sig1;
                
            end
            
       end
       
   end
   
end

temp=G\v';
vmap=zeros(nx,ny);


%backmap 'temp' vector into a physical matrix 
for i=1:nx
    for j=1:ny
        
         n=j+(i-1)*ny;
         vmap(i,j)=temp(n);  
         
    end
        
end


%Electric Field of Voltage Map
[Ex,Ey]=gradient(vmap);


%Current Density of Electric Field
Jx=Ex.*sigma;
Jy=Ey.*sigma;

%Solve for total Current Density
Jtotal=(Jx.^2+Jy.^2).^0.5;

%Solve for Current (Sum Current Density across 1D length)
Jsum=0;

for j=1:ny
    Jsum=Jsum+Jtotal(nx/2,j);
end

end