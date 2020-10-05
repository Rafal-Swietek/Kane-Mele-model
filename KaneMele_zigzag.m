clear all
clc

%initializing parameters & Bravais vectors
a1 = 0.5*[-sqrt(3),3];
a2 = 0.5*[sqrt(3),3];
a = norm(a1)/sqrt(3); %hexagon cell side
%%
%Generating Hamiltonian
d1 = a*[0,1]; %up
d2 = a*[sqrt(3)/2,1/2]; %right-down
d3 = a*[-sqrt(3)/2,1/2]; %left-down
d4 = a*[-1,0]; %left
        sig_x = [0 1;1 0];
        sig_y = [0 -i;i 0];
        sig_z = [1 0;0 -1];
        I = [1 0;0 1];

%Generating Hamiltonian
        %Parameters:
            t = 1;
            LSO = 0.15;
            V = 0.05;
            LR = 0.0;
N = 100; %number of atoms in unit cell
a = a*sqrt(3); %now its unit cell width
figure('Renderer', 'painters', 'Position', [300 100 1000 800]);

grid = 1000;
K = 0:pi/grid/a:2*pi/a;
E = zeros(length(K),2*N);
cl = grid - 220;
cr = grid + 220;
crr = grid - 114;
cll = grid + 114;
for g = 1:length(K)
        k = K(g);
        u = exp(i*k*a); %u is phase going right
            H_NN = zeros(2*N,2*N);
            H_NNN = zeros(2*N,2*N);
            H_V = zeros(2*N,2*N);
            H_R = zeros(2*N,2*N);

            if(mod(N,4) == 0)
                H_NN(2*N-2,2*N) = t*(1+u); %spin up
                H_NN(2*N-3,2*N-1) = t*(1+u); %spin down
            else
                H_NN(2*N-2,2*N) = t*(1+conj(u)); %spin up
                H_NN(2*N-3,2*N-1) = t*(1+conj(u)); %spin down
            end
            H_NNN(2*N-3,2*N-3) = -2*LSO*sin(k*a); %spin down
            H_NNN(2*N-2,2*N-2) = 2*LSO*sin(k*a); %spin up
            H_NNN(2*N-1,2*N-1) = 2*LSO*sin(k*a); %spin down
            H_NNN(2*N,2*N) = -2*LSO*sin(k*a); %spin up

            H_V(2*N-1,2*N-1) = -V; %because width is always even
            H_V(2*N,2*N) = -V;
            H_R(2*N-2,2*N-1) = i*LR*(d3(2)-i*d3(1)) + i*LR*(d2(2)-i*d2(1)); %spin up->down
            H_R(2*N-3,2*N) = i*LR*(d3(2)+i*d3(1)) + i*LR*(d2(2)+i*d2(1)); %spin down->up

  ID = 1; %iterator over sites  
  for ii=1:2:2*N-2
        %Next Nearest Neighbours-----------------
          H_V(ii,ii) = -(-1)^ID*V; 
          H_V(ii+1,ii+1) = -(-1)^ID*V;
        %----------------------------------------
        if(ii<=2*N-5)   
          %Nearest Neighbours & Rashba & Next Nearest Neughbour------------ 
            if(mod(ID,4)==1)
                H_NN(ii,ii+2) = t*(1+conj(u)); %spin down
                H_NN(ii+1,ii+3) = t*(1+conj(u)); %spin up
                H_R(ii,ii+3) = i*LR*( (d2(2)+i*d2(1) + (d3(2)+i*d3(1))*conj(u))); %spin down->up
                H_R(ii+1,ii+2) = i*LR*( (d2(2)-i*d2(1) + (d3(2)-i*d3(1))*conj(u))); %spin up->down

                H_NNN(ii,ii) = -2*LSO*sin(k*a); %spin down
                H_NNN(ii+1,ii+1) = 2*LSO*sin(k*a); %spin up
                H_NNN(ii,ii+4) = i*LSO; %spin down
                H_NNN(ii+1,ii+5) = -i*LSO; %spin up
            end
            if(mod(ID,4)==2)
                H_NN(ii,ii+2) = t; %spin down
                H_NN(ii+1,ii+3) = t; %spin up
                H_R(ii,ii+3) = -i*LR*(d1(2)+i*d1(1)); %spin down->up
                H_R(ii+1,ii+2) = -i*LR*(d1(2)-i*d1(1)); %spin up->down

                H_NNN(ii,ii) = 2*LSO*sin(k*a); %spin down
                H_NNN(ii+1,ii+1) = -2*LSO*sin(k*a); %spin up
                H_NNN(ii,ii+4) = i*LSO; %spin down
                H_NNN(ii+1,ii+5) = -i*LSO; %spin up
            end
            if(mod(ID,4)==3)
                H_NN(ii,ii+2) = t*(1+u); %spin down
                H_NN(ii+1,ii+3) = t*(1+u); %spin up
                H_R(ii,ii+3) = i*LR*( (d2(2)+i*d2(1) + (d3(2)+i*d3(1))*u)); %spin down->up
                H_R(ii+1,ii+2) = i*LR*( (d2(2)-i*d2(1) + (d3(2)-i*d3(1))*u)); %spin up->down

                H_NNN(ii,ii) = -2*LSO*sin(k*a); %spin down
                H_NNN(ii+1,ii+1) = 2*LSO*sin(k*a); %spin up
                H_NNN(ii,ii+4) = -i*LSO; %spin down
                H_NNN(ii+1,ii+5) = i*LSO; %spin up

            end
            if(mod(ID,4)==0)
                H_NN(ii,ii+2) = t; %spin down
                H_NN(ii+1,ii+3) = t; %spin up
                H_R(ii,ii+3) = -i*LR*(d1(2)+i*d1(1)); %spin down->up
                H_R(ii+1,ii+2) = -i*LR*(d1(2)-i*d1(1)); %spin up->down

                H_NNN(ii,ii) = 2*LSO*sin(k*a); %spin down
                H_NNN(ii+1,ii+1) = -2*LSO*sin(k*a); %spin up
                H_NNN(ii,ii+4) = -i*LSO; %spin down
                H_NNN(ii+1,ii+5) = i*LSO; %spin up
            end
      %----------------------------------------------------------------
    end
    ID = ID+1;
  end  
    H = H_NN + H_V + H_NNN + H_R;
    H = H + H';
    E(g,:) = eig(H);
end
    %%
    %Plotting energy
    for ii=1:2*N
        if(ii==N-1)
            plot(K,E(:,ii),'k-');
            plot(K(cl:cr),E(cl:cr,ii),'g-','Linewidth',2);
            plot(K(crr:cll),E(crr:cll,ii),'g-','Linewidth',2);
        elseif(ii==N)
            plot(K,E(:,ii),'k-');
            plot(K(cl:cr),E(cl:cr,ii),'r-','Linewidth',2);%
            plot(K(crr:cll),E(crr:cll,ii),'g-','Linewidth',2);
        elseif(ii==N+1)
            plot(K,E(:,ii),'k-');
            plot(K(cl:cr),E(cl:cr,ii),'g-','Linewidth',2);%
            plot(K(crr:cll),E(crr:cll,ii),'r-','Linewidth',2);
        elseif(ii==N+2)
            plot(K,E(:,ii),'k-');
            plot(K(cl:cr),E(cl:cr,ii),'r-','Linewidth',2);
            plot(K(crr:cll),E(crr:cll,ii),'r-','Linewidth',2);
        else
            plot(K,E(:,ii),'k-');
        end
%         if(ii==N-1)
%             plot(K,E(:,ii),'k-');
%             plot(K(cl:cr),E(cl:cr,ii),'g-','Linewidth',2);
%         elseif(ii==N)
%             plot(K,E(:,ii),'k-');
%             plot(K(cl:cr),E(cl:cr,ii),'g-','Linewidth',2);
%         elseif(ii==N+1)
%             plot(K,E(:,ii),'k-');
%             plot(K(cl:cr),E(cl:cr,ii),'r-','Linewidth',2);
%         elseif(ii==N+2)
%             plot(K,E(:,ii),'k-');
%             plot(K(cl:cr),E(cl:cr,ii),'r-','Linewidth',2);
%         else
%             plot(K,E(:,ii),'k-');
%         end
        %plot(K,E(:,ii),'k-');
        hold on
    end
    drawnow
    %title(sprintf('Energy spectrum for graphene ribbon with zizzag edge for %0.0f atoms width\n using parameters: t=%0.2f, \\lambda=%0.2f and V=%0.2f \n Topological phase transition',N,t,L,V)); 
    title(sprintf('Energy spectrum for graphene ribbon with zizzag edge for %0.0f atoms width\n using parameters: t=%0.2f, \\lambda_{SO}=%0.2f, V=%0.2f and \\lambda_R=%0.2f',N,t,LSO,V,LR)); drawnow
    xlabel(sprintf('k [\\pi/a]')); drawnow
    ylabel(sprintf('E [eV]')); drawnow
    xticks([-3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi]/a); drawnow
    xticklabels({'-3/2','-1','-1/2','0','1/2','1','3/2'}); drawnow
    %axis([min(K) max(K) min(E(:,1)) max(E(:,2*N))]); drawnow
    axis([min(K) max(K) -1 1]); drawnow
    hold off








