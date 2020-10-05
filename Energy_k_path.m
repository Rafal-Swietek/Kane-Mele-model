function Energy_k_path(k_path,parameters,a1,a2)
%Function solving the eigenvalue problem along the Gamma-K-M-K'-Gamma path
%->for any advice using this program mail me at: 77swietek77@gmail.com
%
%Copyright Rafal Schwientek
%--------------------------------------------------------------
%Initial parameters
    t = parameters(1); %hopping integral -> Nearest Neighbours
    V = parameters(2); %staggered potential -> breaks inversion symmetry
    LSO = parameters(3); %spin-orbit coupling -> Next Nearest Neigbours
    LR = parameters(4); %Rashba term -> perpendicular electric field or interaction with substrate
    a = norm(a1)/sqrt(3); %lattice parameter
    
%   % Vectors d_ij in rashba term
%     d1 = a*[0,1]; %up
%     d2 = a*[sqrt(3)/4,-1/2]; %right-down
%     d3 = a*[-sqrt(3)/4,-1/2]; %left-down
%--------------------------------------------------------------    
%Matrices 
    %Pauli matrices   
        sig_x = [0 1;1 0];
        sig_y = [0 -i;i 0];
        sig_z = [1 0;0 -1];
        I = [1 0;0 1];
    %Time-reversal operator    
        T = i*kron(sig_y,I);
    %Gamma matrices (only a bunch, the rest is at the end of the program)
        G1 = kron(sig_x,I);
        G2 = kron(sig_z,I);
        G3 = kron(sig_y,sig_x);
        G4 = kron(sig_y,sig_y);
        G5 = kron(sig_y,sig_z);
        G15 = 1/(2*i)*( G1*G5 - G5*G1 );
        G23 = 1/(2*i)*( G2*G3 - G3*G2 );
        G24 = 1/(2*i)*( G2*G4 - G4*G2 );
        G34 = 1/(2*i)*( G3*G4 - G4*G3 );
        G35 = 1/(2*i)*( G3*G5 - G5*G3 );
        G45 = 1/(2*i)*( G4*G5 - G5*G4 );
%--------------------------------------------------------------
  %Initializing arrays used in the program  
    E1 = zeros(1,length(k_path));
    E2 = zeros(1,length(k_path));
    E3 = zeros(1,length(k_path));
    E4 = zeros(1,length(k_path));
    TRIM  = zeros(1,length(k_path));
    TRIM2  = zeros(1,length(k_path));
    valence = 0;
    cond = 0;
    %Kramers_pair = zeros(1,length(k_path));
    
  %Generating Hamiltonian and diagonilizing him   
    for ii=1:length(k_path)
       kx = k_path(1,ii);
       ky = k_path(2,ii);
       k = [kx, ky];
       k1 = dot(k,a1);  k2 = dot(k,a2);
       %Double Haldane model for different spin
           d45 = t*(1 + cos(k1) + cos(k2)); % Re( NN )
           d35 = -t*( sin(k1) + sin(k2) ); % Im( NN )
           d15 =  2*LSO*( sin(k1) - sin(k2) - sin(k1-k2)); % spin-orbit coupling
           d34 = V; %stagerred potential
        
           H = d45*G45 + d35*G35 + d15*G15 + d34*G34; %d2*G2 if V is the potential between differen spin
       
       %Rashba term
           x = (k1 + k2)/2; %temporary variables
           y = (k2 - k1)/2;
           d3 = sqrt(3)*LR*sin(y)*cos(x);
           d4 = -sqrt(3)*LR*sin(x)*sin(y);
           d23 = -LR*sin(x)*cos(y);
           d24 = LR*( 1 - cos(x)*cos(y) );
           
       H_R = d3*G3 + d4*G4 + d23*G23 + d24*G24;
       H = H + H_R + H_R';
       [~,E] = eig(H);
       %Properly reordering the bands
           if(( abs(E(1,1) - E(2,2)) <= 1e-14) && (kx > 0) )
               valence = valence+1;
           end
           if((abs(E(3,3) - E(4,4)) <= 1e-14) && (kx > 0))
               cond = cond+1;
           end
           if(mod(valence,2) == 1)
               temp = E(2,2);
               E(2,2) = E(1,1);
               E(1,1) = temp;
           end
           if(mod(cond,2) == 1)
               temp = E(4,4);
               E(4,4) = E(3,3);
               E(3,3) = temp;
           end
       %-------------------
       E1(ii) = E(1,1);
       E2(ii) = E(2,2);
       E3(ii) = E(3,3);
       E4(ii) = E(4,4);
       %------------------------------------
       k = [-kx, -ky];
       k1 = dot(k,a1);  k2 = dot(k,a2);
       %Double Haldane model for different spin
           d45 = t*(1 + cos(k1) + cos(k2)); % Re( NN )
           d35 = -t*( sin(k1) + sin(k2) ); % Im( NN )
           d15 =  2*LSO*( sin(k1) - sin(k2) - sin(k1-k2)); % spin-orbit coupling
           d34 = V; %stagerred potential
        
           H2 = d45*G45 + d35*G35 + d15*G15 + d34*G34; %d2*G2 if V is the potential between differen spin
       
       %Rashba term
           x = (k1 + k2)/2; %temporary variables
           y = (k2 - k1)/2;
           d3 = sqrt(3)*LR*sin(y)*cos(x);
           d4 = -sqrt(3)*LR*sin(x)*sin(y);
           d23 = -LR*sin(x)*cos(y);
           d24 = LR*( 1 - cos(x)*cos(y) );
           
       H_R = d3*G3 + d4*G4 + d23*G23 + d24*G24;
       H2 = H2 + H_R + H_R';
       [~,E] = eig(H2);
       %Properly reordering the bands
           if(( abs(E(1,1) - E(2,2)) <= 1e-14) && (kx > 0) )
               valence = valence+1;
           end
           if((abs(E(3,3) - E(4,4)) <= 1e-14) && (kx > 0))
               cond = cond+1;
           end
           if(mod(valence,2) == 1)
               temp = E(2,2);
               E(2,2) = E(1,1);
               E(1,1) = temp;
           end
           if(mod(cond,2) == 1)
               temp = E(4,4);
               E(4,4) = E(3,3);
               E(3,3) = temp;
           end
       %------------------------------------
       
       %Hamiltonian and Time-reversal operator commutation - finding TRIM
            TRIM(ii) = norm( H - T*conj(H)*T^(-1) ) ;
            TRIM2(ii) = norm( H2 - T*conj(H)*T^(-1) ) ;
       %-------------------------------------------------------------- 

    end
%--------------------------------------------------------------    
   %Plotting energy solution
    figure(135);
    y = min(E1):1:max(E4); 
    x1 = length(k_path)/3*ones(1,length(y))+1; 
    x2 = length(k_path)/3*1.5*ones(1,length(y))+1; 
    x3 = length(k_path)/3*2*ones(1,length(y))+1; %dashed vertical lines
    
    plot(x1,y,'k--',x2,y,'k--',x3,y,'k--',1:length(k_path),E1,'k-', 1:length(k_path),E2,'r-'); hold on
    plot(x1,y,'k--',x2,y,'k--',x3,y,'k--',1:length(k_path),E3,'k-', 1:length(k_path),E4,'r-'); drawnow
    xticks([0 length(k_path)/3+1 1.5*length(k_path)/3+1 2*length(k_path)/3+1 length(k_path)]); drawnow
    xticklabels({'\Gamma','K','M','K"','\Gamma'}); drawnow
    axis([1 length(k_path) min(E1) max(E4)]); drawnow
    ylabel('E [eV]'); drawnow
    title(sprintf('Energy spectrum for given k-path:\\Gamma-K-M-K"-\\Gamma \n using parameters: t =%1.1f , V=%1.2f, \\lambda_{SO} = %1.2f and \\lambda_R = %1.2f ',t,V,LSO,LR)); drawnow
    hold off
%--------------------------------------------------------------   
   %Plotting commutation
        figure(11);
        plot(1:length(k_path),TRIM,'ro-',1:length(k_path),TRIM2,'ko-');
        title(sprintf('Commutation of Hamiltonian \nand Time-reversal operator \nalong \\Gamma-K-M-K"-\\Gamma path\n'));
        xticks([0 length(k_path)/3+1 1.5*length(k_path)/3+1 2*length(k_path)/3+1 length(k_path)]);
        xticklabels({'\Gamma','K','M','K"','\Gamma'});
        ylabel(sprintf(' Re([ H,T ]) '));
        %axis([0 length(k_path) 0 10]); %draw one with axis on and with axis off
        legend('TRIM','Comutation');
%--------------------------------------------------------------   
%    %Plotting Kramers pairs
%         figure(22);
%         plot(1:length(k_path),Kramers_pair,'ro-');
%         title(sprintf('Kramers pairs \nalong \\Gamma-K-M-K"-\\Gamma path\n'));
%         xticks([0 length(k_path)/3+1 1.5*length(k_path)/3+1 2*length(k_path)/3+1 length(k_path)]);
%         xticklabels({'\Gamma','K','M','K"','\Gamma'});
%         ylabel(sprintf(' H(-k) - TH(k)T^-1'));

end

%--------------------------------------------------------------
% %Other Gamma matrices
%         G12 = -0.5*i*( G1*G2 - G2*G1 ); Unecessary for now
%         G13 = -0.5*i*( G1*G3 - G3*G1 );
%         G14 = 1/(2*i)*( G1*G4 - G4*G1 );
%         G15 = -0.5*i*( G1*G5 - G5*G1 );
%         G25 = -0.5*i*( G2*G5 - G5*G2 );
%         G34 = 1/(2*i)*( G3*G4 - G4*G3 );
%         G35 = -0.5*i*( G3*G5 - G5*G3 );
%--------------------------------------------------------------
% %Hamiltonian as double Haldane model
%        d_x = t*(1 + cos(k1) + cos(k2));
%        d_y = t*( sin(k1) + sin(k2) );
%        d_z = V - 2*LSO*( sin(k1) - sin(k2) - sin(k1-k2)) ;
%        H_Haldane = [d_z, d_x-i*d_y; d_x+i*d_y, -d_z];
%        A = [1 0; 0 0];  B = [0 0;0 1];
%        H = kron(A,H_Haldane) + kron(B,conj(H_Haldane));
%-------------------------------------------------------------- 
% For finding Kramers pairs one has to write the hamiltonian down for -k!!
%
%        %Hamiltonian and Time-reversal operator commutation - finding Kramers pairs
%             Kramers_pair(ii) = det( H2*T - T*conj(H) ) ;
%-------------------------------------------------------------- 
%--------------------------------------------------------------     
% % Creating .gif file for animation
%     frame = getframe(1);
%     im = frame2im(frame);
%     [imind, cm] = rgb2ind(im,256);
%     if L == -L_0;
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%     else
%           imwrite(imind,cm,filename,'gif','WriteMode','append');
%     end
%--------------------------------------------------------------
%--------------------------------------------------------------
% %Double Haldane model for different spin - WRONG MAY IT BE!!!
%            dx = t*(1+cos(k1)+cos(k2));
%            dy = -t*( sin(k1) + sin(k2) );
%            dz =  V - 2*LSO*( sin(k1) - sin(k2) - sin(k1-k2));
%         
%            H = dx*G45 + dy*G35 + dz*G34;
%--------------------------------------------------------------
%--------------------------------------------------------------
% %Anihilation operators on each site with each spin
%         a_Adown = [0 0 0 0;1 0 0 0];
%         a_Aup = [0 0 0 0;0 0 1 0];
%         a_Bdown = [0 0 0 0;0 1 0 0];
%         a_Bup = [0 0 0 0;0 0 0 1];
%--------------------------------------------------------------
% %Rashba term using anihilation/ creation operators
%            H_R = i*LR*(a_Adown'*(d1(2)+i*d1(1))*a_Bup + a_Aup'*(d1(2)-i*d1(1))*a_Bdown);
%            H_R = H_R + i*LR*(a_Adown'*(d2(2)+i*d2(1))*a_Bup + a_Aup'*(d2(2)-i*d2(1))*a_Bdown)*exp(-i*k1);
%            H_R = H_R + i*LR*(a_Adown'*(d3(2)+i*d3(1))*a_Bup + a_Aup'*(d3(2)-i*d3(1))*a_Bdown)*exp(-i*k2);
