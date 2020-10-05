function Energy3D(parameters,a1,a2, grid,KK)
%--------------------------------------------------------------
%Initial parameters
    t = parameters(1);
    V = parameters(2);
    LSO = parameters(3);
    LR = parameters(4);
    a = norm(a1)/sqrt(3); %lattice parameter
%   % Vectors d_ij in rashba term
%     d1 = a*[0,0.5]; %up
%     d2 = 0.5*a*[sqrt(3)/4,1/2]; %right-down
%     d3 = 0.5*a*[-sqrt(3)/4,1/2]; %left-down
%     d4 = a*[-0.5,0]; %left
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
    kx = -2*pi/norm(a1):pi/norm(a1)/grid:2*pi/norm(a1); %kx grid
    ky = -2*pi/norm(a2):pi/norm(a2)/grid:2*pi/norm(a2); %ky grid
        E1 = zeros(length(kx),length(ky));
        E2 = zeros(length(kx),length(ky));
        E3 = zeros(length(kx),length(ky));
        E4 = zeros(length(kx),length(ky));
        out_re = zeros(length(kx),length(ky));
  %Simulation
    for jj = 1:length(ky)
        k_y = ky(jj);
        for ii = 1:length(kx)
            k_x = kx(ii);
            k1 = k_x*a1(1)+k_y*a1(2);
            k2 = k_x*a2(1)+k_y*a2(2);
       %Double Haldane model for different spin
           d45 = t*(1 + cos(k1) + cos(k2)); % Re( NN )
           d35 = -t*( sin(k1) + sin(k2) ); % Im( NN )
           d15 =  4*LSO*( sin(k1) - sin(k2) - sin(k1-k2)); % spin-orbit coupling
           d34 = 2*V; %stagerred potential
        
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
       [~,energy] = eig(H);
       E1(ii,jj) = energy(1,1);
       E2(ii,jj) = energy(2,2);
       E3(ii,jj) = energy(3,3);
       E4(ii,jj) = energy(4,4);
       %Hamiltonian and Time-reversal operator commutation
            out_re(ii,jj) = norm( (H - T*conj(H)*T^(-1)));
        end
    end
%--------------------------------------------------------------       
% %Plotting 3D plot E(kx,ky) and 2D plots (E(kx,0) & E(0,ky)
%     figure(2); % 3D plot
%     Colour1 = gradient(E1);
%     Colour2 = gradient(E4);
%     mesh(kx,ky,E1,Colour1); hold on
%     mesh(kx,ky,E2,Colour2); hold on
%     mesh(kx,ky,E3, Colour1); hold on 
%     mesh(kx,ky,E4, Colour2); hold on
%     plot(KK(:,1)/sqrt(3),KK(:,2)/sqrt(3),'k-','Linewidth',2);hold off
%     xticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a1));   
%     xticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
%     yticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a2));   
%     yticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
%     xlabel(sprintf('k_x [\\pi/a]')); ylabel(sprintf('k_y [\\pi/a]')); zlabel('E [eV]');
%     title(sprintf('Energy spectrum \n using parameters: t =%1.1f, V = %1.2f and \\lambda = %1.2f ',t,V,LSO));
%    
    %Plotting commutation
    figure(1111);
    mesh(kx,ky,out_re);hold on
    plot(KK(:,1)*1.15,KK(:,2)*1.15,'k-','Linewidth',2);hold off
    
    title(sprintf('Commutation of Hamiltonian \nand Time-reversal operator\n for parameters: t =%1.1f , V=%1.2f, \\lambda_{SO} = %1.2f and \\lambda_R = %1.2f\n\n\n\n ',t,V,LSO,LR));
    xticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a1));   
    xticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
    yticks([-2*pi -3/2*pi -pi -pi/2 0 pi/2 pi 3/2*pi 2*pi]/norm(a2));   
    yticklabels({'-2','-3/2','-1','-1/2','0','1/2','1','3/2','2'});
    xlabel(sprintf('k_x [\\pi/a]')); ylabel(sprintf('k_y [\\pi/a]'));
    zlabel(sprintf(' Re([ H,T ]) '));
    axis([-2*pi/norm(a1) 2*pi/norm(a1) -2*pi/norm(a1) 2*pi/norm(a1) 0 max(max(out_re))]); %draw one with axis on and with axis off

end

%        %Rashba term
%            d3 = LR*(1-cos((k1+k2)/2)*cos((k2-k1)/2));
%            d4 = -sqrt(3)*LR*sin((k2-k1)/2)*sin((k1+k2)/2);
%            d23 = -LR*cos((k2-k1)/2)*sin((k1+k2)/2);
%            d24 = sqrt(3)*LR*sin((k2-k1)/2)*sin((k1+k2)/2);
%         
%            H_R = d3*G3 + d4*G4 + d23*G23 + d24*G24; 

