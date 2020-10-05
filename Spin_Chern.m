function Spin_Chern(parameters,a1,a2,grid,M)
%This function allows to calculate the spin chern numbers for both spin-up
%and spin-down bands
%   The formula used in this calulation is as follows:
%   c_up,down = sum_{k in EBZ}(berry_up,down(kx,ky) )
%-------------------------------------------------------------
%Initial parameters
    t = parameters(1); %Nearest Neighbours
    V = parameters(2); %staggered potential 
    LSO = parameters(3); %Next Nearest Neigbours
    LR = parameters(4); %Rashba term
%--------------------------------------------------------------    
%Matrices 
    %Pauli matrices   
        sig_x = [0 1;1 0];
        sig_y = [0 -i;i 0];
        sig_z = [1 0;0 -1];
        I = [1 0;0 1];
        
        T = i*kron(sig_y,I); %Time-reversal operator
    %Gamma matrices
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
%Calculations
    dkx = pi/norm(a1)/grid;
    dky = pi/norm(a2)/grid;
    k_x = -pi/norm(a1):dkx:2*pi/norm(a1); %kx grid
    k_y = -pi/norm(a2):dky:2*pi/norm(a2); %ky grid
    
    Gamma1 = 2*M(1,:);
    Gamma2 = 2*M(2,:);
    
    c_Vdown = 0; %Chern numbers for each band
    c_Vup = 0;
    c_Cdown = 0;
    c_Cup = 0;
for jj = 1:length(k_y)
    ky = k_y(jj);
    for ii = 1:length(k_x)
        kx = k_x(ii);
        if((kx>=0) && (kx<=Gamma1(1))) %Conditions for being inside of EBZ
           if( ky >= Gamma1(2)/Gamma1(1)*kx )
              if(ky <=( (Gamma2(2)-Gamma1(2))/(Gamma2(1)-Gamma1(1))*kx + Gamma2(2)) )
                  %Eigenvectors first----------------
                    k = [kx, ky];
                    k1 = dot(k,a1);  k2 = dot(k,a2);
                   %Double Haldane model for different spin
                       d45 = t*(1 + cos(k1) + cos(k2));    % Re( NN )
                       d35 = t*( sin(k1) + sin(k2) ); % Im( NN )
                       d15 =  2*LSO*( sin(k1) - sin(k2) - sin(k1-k2)); % spin-orbit coupling
                       d34 = V; %stagerred potential
                       H = d45*G45 + d35*G35 + d15*G15 + d34*G34;
                   %Rashba term
                       x = (k1 + k2)/2;
                       y = (k2 - k1)/2;
                       d3 = sqrt(3)*LR*sin(y)*cos(x);
                       d4 = -sqrt(3)*LR*sin(x)*sin(y);
                       d23 = -LR*sin(x)*cos(y);
                       d24 = LR*( 1 - cos(x)*cos(y) );
                   H_R = d3*G3 + d4*G4 + d23*G23 + d24*G24;
                   H = H + H_R;
                   H = H + H';
                   [Pk1,~] = eig(H);
                   Pk1(2,:) = Pk1(2,:)*exp(i*(k1+k2)/2);
                   Pk1(4,:) = Pk1(4,:)*exp(i*(k1+k2)/2);
                    %Properly reordering the bands
                        if(abs( Pk1(1,1) )==0)
                           temp = Pk1(:,1);
                           Pk1(:,1) = Pk1(:,2);
                           Pk1(:,2) = temp;
                        end
                        if(abs( Pk1(1,4) )==0)
                           temp = Pk1(:,3);
                           Pk1(:,3) = Pk1(:,4);
                           Pk1(:,4) = temp;
                        end
                    %-------------------
             %-----------------------------------------
             %Eigenvectors second----------------------
                    k = [kx+dkx, ky];
                    k1 = dot(k,a1);  k2 = dot(k,a2);
                   %Double Haldane model for different spin
                       d45 = t*(1 + cos(k1) + cos(k2));    % Re( NN )
                       d35 = t*( sin(k1) + sin(k2) ); % Im( NN )
                       d15 =  2*LSO*( sin(k1) - sin(k2) - sin(k1-k2)); % spin-orbit coupling
                       d34 = V; %stagerred potential
                       H = d45*G45 + d35*G35 + d15*G15 + d34*G34;
                   %Rashba term
                       x = (k1 + k2)/2;
                       y = (k2 - k1)/2;
                       d3 = sqrt(3)*LR*sin(y)*cos(x);
                       d4 = -sqrt(3)*LR*sin(x)*sin(y);
                       d23 = -LR*sin(x)*cos(y);
                       d24 = LR*( 1 - cos(x)*cos(y) );
                   H_R = d3*G3 + d4*G4 + d23*G23 + d24*G24;
                   H = H + H_R;
                   H = H + H';
                   [Pk2,~] = eig(H);
                   Pk2(2,:) = Pk2(2,:)*exp(i*(k1+k2)/2);
                   Pk2(4,:) = Pk2(4,:)*exp(i*(k1+k2)/2);
                   %Properly reordering the bands

                        if(abs( Pk2(1,1) )==0)
                           temp = Pk2(:,1);
                           Pk2(:,1) = Pk2(:,2);
                           Pk2(:,2) = temp;
                        end
                        if(abs( Pk2(1,4) )==0)
                           temp = Pk2(:,3);
                           Pk2(:,3) = Pk2(:,4);
                           Pk2(:,4) = temp;
                        end
                    %-------------------
             %Eigenvectors third--------------------
                    k = [kx+dkx, ky+dky];
                    k1 = dot(k,a1);  k2 = dot(k,a2);
                   %Double Haldane model for different spin
                       d45 = t*(1 + cos(k1) + cos(k2));    % Re( NN )
                       d35 = t*( sin(k1) + sin(k2) ); % Im( NN )
                       d15 =  2*LSO*( sin(k1) - sin(k2) - sin(k1-k2)); % spin-orbit coupling
                       d34 = V; %stagerred potential
                       H = d45*G45 + d35*G35 + d15*G15 + d34*G34;
                   %Rashba term
                       x = (k1 + k2)/2;
                       y = (k2 - k1)/2;
                       d3 = sqrt(3)*LR*sin(y)*cos(x);
                       d4 = -sqrt(3)*LR*sin(x)*sin(y);
                       d23 = -LR*sin(x)*cos(y);
                       d24 = LR*( 1 - cos(x)*cos(y) );
                   H_R = d3*G3 + d4*G4 + d23*G23 + d24*G24;
                   H = H + H_R;
                   H = H + H';
                   [Pk3,~] = eig(H);
                   Pk3(2,:) = Pk3(2,:)*exp(i*(k1+k2)/2);
                   Pk3(4,:) = Pk3(4,:)*exp(i*(k1+k2)/2);
                   %Properly reordering the bands
                        if(abs( Pk3(1,1) )==0)
                           temp = Pk3(:,1);
                           Pk3(:,1) = Pk3(:,2);
                           Pk3(:,2) = temp;
                        end
                        if(abs( Pk3(1,4) )==0)
                           temp = Pk3(:,3);
                           Pk3(:,3) = Pk3(:,4);
                           Pk3(:,4) = temp;
                        end
                   %-------------------
              %Eigenvectors fourth--------------------
                    k = [kx, ky+dky];
                    k1 = dot(k,a1);  k2 = dot(k,a2);
                   %Double Haldane model for different spin
                       d45 = t*(1 + cos(k1) + cos(k2));    % Re( NN )
                       d35 = t*( sin(k1) + sin(k2) ); % Im( NN )
                       d15 =  2*LSO*( sin(k1) - sin(k2) - sin(k1-k2)); % spin-orbit coupling
                       d34 = V; %stagerred potential
                       H = d45*G45 + d35*G35 + d15*G15 + d34*G34;
                   %Rashba term
                       x = (k1 + k2)/2;
                       y = (k2 - k1)/2;
                       d3 = sqrt(3)*LR*sin(y)*cos(x);
                       d4 = -sqrt(3)*LR*sin(x)*sin(y);
                       d23 = -LR*sin(x)*cos(y);
                       d24 = LR*( 1 - cos(x)*cos(y) );
                   H_R = d3*G3 + d4*G4 + d23*G23 + d24*G24;
                   H = H + H_R;
                   H = H + H';
                   [Pk4,~] = eig(H);
                   Pk4(2,:) = Pk4(2,:)*exp(i*(k1+k2)/2);
                   Pk4(4,:) = Pk4(4,:)*exp(i*(k1+k2)/2);
                   %Properly reordering the bands

                        if(abs( Pk4(1,1) )==0)
                           temp = Pk4(:,1);
                           Pk4(:,1) = Pk4(:,2);
                           Pk4(:,2) = temp;
                        end
                        if(abs( Pk4(1,4) )==0)
                           temp = Pk4(:,3);
                           Pk4(:,3) = Pk4(:,4);
                           Pk4(:,4) = temp;
                        end
                    %------------------- 
                    %---------------------------------------
                    %Calculating chern numbers for each band
                        x12 = dot(Pk1(:,1),Pk2(:,1));
                        x23 = dot(Pk2(:,1),Pk3(:,1));
                        x34 = dot(Pk3(:,1),Pk4(:,1));
                        x41 = dot(Pk4(:,1),Pk1(:,1));
                        c_Vdown = c_Vdown + angle(x12/abs(x12)*x23/abs(x23)*x34/abs(x34)*x41/abs(x41));
                        
                        x12 = dot(Pk1(:,2),Pk2(:,2));
                        x23 = dot(Pk2(:,2),Pk3(:,2));
                        x34 = dot(Pk3(:,2),Pk4(:,2));
                        x41 = dot(Pk4(:,2),Pk1(:,2));
                        c_Vup = c_Vup + angle(x12/abs(x12)*x23/abs(x23)*x34/abs(x34)*x41/abs(x41));
                        
                        x12 = dot(Pk1(:,3),Pk2(:,3));
                        x23 = dot(Pk2(:,3),Pk3(:,3));
                        x34 = dot(Pk3(:,3),Pk4(:,3));
                        x41 = dot(Pk4(:,3),Pk1(:,3));
                        c_Cdown = c_Cdown + angle(x12/abs(x12)*x23/abs(x23)*x34/abs(x34)*x41/abs(x41));
                        
                        x12 = dot(Pk1(:,4),Pk2(:,4));
                        x23 = dot(Pk2(:,4),Pk3(:,4));
                        x34 = dot(Pk3(:,4),Pk4(:,4));
                        x41 = dot(Pk4(:,4),Pk1(:,4));
                        c_Cup = c_Cup + angle(x12/abs(x12)*x23/abs(x23)*x34/abs(x34)*x41/abs(x41));
                    %-----------------------------
              end %endif
           end %endif
        end %endif
    end
end
% fprintf('[c_Vdown, c_Vup, c_Cdown, c_Cup]');
% X = [c_Vdown, c_Vup, c_Cdown, c_Cup]/(2*pi)
C_up = (c_Vup + c_Cup)/(2*pi)
C_down = (c_Vdown + c_Cdown)/(2*pi)
Cs = (C_up - C_down)/2
% Spin_Chern_valence = 0.5*( c_Vup - c_Vdown)/(2*pi)
% Spin_Chern_conducting = 0.5*( c_Cup - c_Cdown)/(2*pi)
end

