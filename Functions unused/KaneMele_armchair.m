%initializing parameters & Bravais vectors
a1 = 0.5*[-sqrt(3),3];
a2 = 0.5*[sqrt(3),3];
a = sqrt(3)*norm(a1); %unit cell width is 3a0
a0 = a/3;
d1 = a0*[0,1]; %up
d2 = a0*[sqrt(3)/2,1/2]; %right-down
d3 = a0*[-sqrt(3)/2,1/2]; %left-down
d4 = a0*[-1,0]; %left
        sig_x = [0 1;1 0];
        sig_y = [0 -i;i 0];
        sig_z = [1 0;0 -1];
        I = [1 0;0 1];

%Generating Hamiltonian
        %Parameters:
            t = 1;      %NN hopping
            LSO = 0.0; %NNN, spin orbit coupling
            V = 0.0;    %stagered potential
            LR = 0.0;   %rashba term (external field)
        %
figure('Renderer', 'painters', 'Position', [200 100 800 700]);
%filename = 'ChangingL_V=0,2.gif';
grid = 200;

K = -pi/a:pi/grid/a:pi/a;
n = 17; %width in number of atoms
%for n=10:1:30
N = 2*n; %number of atoms in unit cell

E = zeros(length(K),2*N);
    for g = 1:length(K)
        k = K(g);
        u = exp(i*k*a); %u is phase going right
        
            H_NN  = zeros(2*N,2*N);
            H_NNN = zeros(2*N,2*N);
            H_V = zeros(2*N,2*N);
            H_R = zeros(2*N,2*N);
     
           ID = 1;
           for ii=1:2:N-3
                H_NN(ii,ii+2) = t;
                H_NN(ii+1,ii+3) = t;
                H_NN(2*N-ii+1,2*N-ii-1) = t; %spin up
                H_NN(2*N-ii,2*N-ii-2) = t; %spin down
              
                if(mod(ID,2)==0)%even site
                  %Nearest neighbours---------
                    H_NN(ii,2*N-ii) = t;
                    H_NN(ii+1,2*N-ii+1) = t;
                  %Stagerred potential--------  
                    H_V(ii,ii) = V; %spin down
                    H_V(ii+1,ii+1) = V; %spin up
                  %Next Nearest Neighbours----
                   if(ii<=N-5)
                        H_NNN(ii,ii+4) = -i*LSO; %2->4; spin down
                        H_NNN(ii+1,ii+5) = i*LSO; %spin up
                        H_NNN(ii,2*N-ii-2) = i*LSO*(1+conj(u)); %2->12; spin down
                        H_NNN(ii,2*N-ii-2) = -i*LSO*(1+conj(u)); %spin up
                        H_NNN(2*N-ii,2*N-ii-4) = -i*LSO;%14->12; spin down
                        H_NNN(2*N-ii+1,2*N-ii-3) = i*LSO;%spin up
                        H_NNN(2*N-ii,ii+2) = i*LSO*(1+u); %14->2; spin down
                        H_NNN(2*N-ii+1,ii+3) = -i*LSO*(1+u); %spin up
                   end
                  %Rashba term----------------  
                    H_R(ii,ii+3) = i*LR*(d3(2)+i*d3(1)); %spin down->up
                    H_R(ii+1,ii+2) = i*LR*(d3(2)-i*d3(1)); %spin up->down
                    
                    H_R(2*N-ii+1,2*N-ii-2) = i*LR*(d2(2)+i*d2(1)); %spin down->up                       
                    H_R(2*N-ii,2*N-ii-1) = i*LR*(d2(2)-i*d2(1)); %spin up->down
                    
                    H_R(ii,2*N-ii-2) = i*LR*(d2(2)+i*d2(1)); %spin down->up
                    H_R(ii+1,2*N-ii-1) = i*LR*(d2(2)-i*d2(1)); %spin up->down
                else %odd site
                  %Nearest neighbours---------
                    H_NN(ii,2*N-ii) = t*conj(u);
                    H_NN(ii+1,2*N-ii+1) = t*conj(u);
                  %Stagerred potential--------  
                    H_V(ii,ii) = -V; %spin down
                    H_V(ii+1,ii+1) = -V; %spin up
                  %Next Nearest Neighbours---- 
                   if(ii<=N-5)
                        H_NNN(ii,ii+4) = i*LSO; %1->3; spin down
                        H_NNN(ii+1,ii+5) = -i*LSO; %spin up
                        H_NNN(ii,2*N-ii-2) = -i*LSO*(1+conj(u)); %1->13; spin down
                        H_NNN(ii,2*N-ii-2) = i*LSO*(1+conj(u)); %spin up
                        H_NNN(2*N-ii,2*N-ii-4) = i*LSO;%13->11; spin down
                        H_NNN(2*N-ii+1,2*N-ii-3) = -i*LSO;%spin up
                        H_NNN(2*N-ii,ii+2) = -i*LSO*(1+u); %13->3; spin down
                        H_NNN(2*N-ii+1,ii+3) = i*LSO*(1+u); %spin up
                   end
                  %Rashba term---------------- 
                    H_R(ii,ii+3) = i*LR*(d2(2)+i*d2(1)); %spin down->up
                    H_R(ii+1,ii+2) = i*LR*(d2(2)-i*d2(1)); %spin up->down
                    
                    H_R(2*N-ii+1,2*N-ii-2) = i*LR*(d3(2)+i*d3(1)); %spin down->up                       
                    H_R(2*N-ii,2*N-ii-1) = i*LR*(d3(2)-i*d3(1)); %spin up->down
                    
                    H_R(ii,2*N-ii-2) = i*LR*(d2(2)+i*d2(1))*conj(u); %spin down->up
                    H_R(ii+1,2*N-ii-1) = i*LR*(d2(2)-i*d2(1))*conj(u); %spin up->down
                end
                ID = ID+1; %determines odd or even site
           end
       %-----------------------------
       %Edge points
        if(mod(n,2)==1) %odd width
            H_NN(2*n,2*n+2) = t*conj(u);
            H_NN(2*n-1,2*n+1) = t*conj(u);
            H_NNN(2*n-3,2*n+1) = i*LSO*(1+conj(u));
            H_NNN(2*n-2,2*n+2) = i*LSO*(1+conj(u));
            
            H_R(2*n-1,2*n+2) = i*LR*(d4(2)+i*d4(1))*conj(u); %spin down->up
            H_R(2*n,2*n+1) = i*LR*(d4(2)-i*d4(1))*conj(u); %spin up->down
        else %even width
            H_NN(2*n,2*n+2) = t;
            H_NN(2*n-1,2*n+1) = t;
            H_NNN(2*n-3,2*n+1) = -i*LSO*(1+conj(u));
            H_NNN(2*n-2,2*n+2) = -i*LSO*(1+conj(u));
            
            H_R(2*n-1,2*n+2) = -i*LR*(d4(2)+i*d4(1)); %spin down->up
            H_R(2*n,2*n+1) = -i*LR*(d4(2)-i*d4(1)); %spin up->down
        end
        
        H = H_NN + H_V + H_NNN + H_R;
        H = H + H';
        E(g,:) = eig(H);
    end
    %%
    %Plotting energy
    for ii=1:2*N
        if((ii==N-1) || (ii==N+1))
            plot(K,E(:,ii),'r-');
        elseif((ii==N) || (ii==N+2))
            plot(K,E(:,ii),'r-');
        else
            plot(K,E(:,ii),'k-','Linewidth',1);
        end
        hold on
    end
    drawnow
    title(sprintf('Energy spectrum for graphene ribbon with armchair edge for %0.0f atoms width\n using parameters: t=%0.2f, \\lambda_{SO}=%0.2f, V=%0.2f and \\lambda_R=%0.2f',N/2,t,LSO,V,LR)); drawnow
    xlabel(sprintf('k [\\pi/a]')); 
    ylabel(sprintf('E [eV]')); 
    xticks([-pi -pi/2 0 pi/2 pi]/a); 
    xticklabels({'-1','-1/2','0','1/2','1'});
    axis([min(K) max(K) -1 1]); drawnow 
    %axis([min(K) max(K) min(E(:,1)) max(E(:,2*N))]); drawnow
    hold off

 %end









