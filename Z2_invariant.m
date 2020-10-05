function Z2 = Z2_invariant(EBZ,parameters,a1,a2,M, grid)
% This function allows to determine the Z2 index for the Kane-Mele model
% using the formula:
%       v = 1/2pi*( int_EBZ(F*d^2k) - int_dEBZ( A*dk ) ) mod 2
%
%--------------------------------------------------------------
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
%Calculating second term: int_dEBZ( A*dk )  
L = zeros(4,4); 
     %Eigenvectors for k=Gamma
        k = [0,0];
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
       [Pk,~] = eig(H);
       Pk(2,:) = Pk(2,:)*exp(i*(k1+k2)/2);
       Pk(4,:) = Pk(4,:)*exp(i*(k1+k2)/2);
A = 1;
for ii = 2:length(EBZ)
       %Eigenvectors for k
        kx = EBZ(1,ii);
        ky = EBZ(2,ii);
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
       [Pk_next,~] = eig(H);
       Pk_next(2,:) = Pk_next(2,:)*exp(i*(k1+k2)/2);
       Pk_next(4,:) = Pk_next(4,:)*exp(i*(k1+k2)/2);
        %Properly reordering the bands
            if(abs( Pk_next(1,1) )==0)
               temp = Pk_next(:,1);
               Pk_next(:,1) = Pk_next(:,2);
               Pk_next(:,2) = temp;
            end
            if(abs( Pk_next(1,4) )==0)
               temp = Pk_next(:,3);
               Pk_next(:,3) = Pk_next(:,4);
               Pk_next(:,4) = temp;
            end
        %-------------------
%       %Calculating matrix
%        for m=1:4
%            for n=1:4
%                %w_mn = < um(-k)|T|un(k) >
%                if( abs( conj(Pk(:,n)')*G2*Pk(:,n)+conj(Pk(:,m)')*G2*Pk(:,m))==0 )
%                    L(m,n) = dot(Pk(:,m),Pk_next(:,n));
%                end
%            end
%        end
       %L = Pk'*Pk_next;
       L(1,2) = dot(Pk(:,1),Pk_next(:,2));
       L(3,4) = dot(Pk(:,3),Pk_next(:,4));
       L(2,1) = dot(Pk(:,2),Pk_next(:,1));
       L(4,3) = dot(Pk(:,4),Pk_next(:,3));
       A = A*det(L)/abs(det(L));
      %Assigning values for next iteration 
       Pk = Pk_next;
end
A = angle(A);
%--------------------------------------------------------------     
%Calculating second term: int_EBZ(F*d^2k)
    dkx = pi/norm(a1)/grid;
    dky = pi/norm(a2)/grid;
    k_x = -pi/norm(a1):dkx:2*pi/norm(a1); %kx grid
    k_y = -pi/norm(a2):dky:2*pi/norm(a2); %ky grid
    
    Gamma1 = 2*M(1,:);
    Gamma2 = 2*M(2,:);
    F = 0; %berry curvature;
for jj = 1:length(k_y)
    ky = k_y(jj);
    for ii = 1:length(k_x)
        kx = k_x(ii);
        F_temp = 0;
      if((kx>=0)) %Conditions for being inside of EBZ
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
             L(1,2) = dot(Pk1(:,1),Pk2(:,2));
             L(3,4) = dot(Pk1(:,3),Pk2(:,4));
             L(2,1) = dot(Pk1(:,2),Pk2(:,1));
             L(4,3) = dot(Pk1(:,4),Pk2(:,3));
             F_temp = det(L) / abs(det(L));
             
             L(1,2) = dot(Pk2(:,1),Pk3(:,2));
             L(3,4) = dot(Pk2(:,3),Pk3(:,4));
             L(2,1) = dot(Pk2(:,2),Pk3(:,1));
             L(4,3) = dot(Pk2(:,4),Pk3(:,3));
             F_temp = F_temp*det(L) / abs(det(L));
             
             L(1,2) = dot(Pk3(:,1),Pk4(:,2));
             L(3,4) = dot(Pk3(:,3),Pk4(:,4));
             L(2,1) = dot(Pk3(:,2),Pk4(:,1));
             L(4,3) = dot(Pk3(:,4),Pk4(:,3));
             F_temp = F_temp*det(L) / abs(det(L));
             
             L(1,2) = dot(Pk4(:,1),Pk1(:,2));
             L(3,4) = dot(Pk4(:,3),Pk1(:,4));
             L(2,1) = dot(Pk4(:,2),Pk1(:,1));
             L(4,3) = dot(Pk4(:,4),Pk1(:,3));
             F_temp = F_temp*det(L) / abs(det(L));
%                    %Calculating matrix
%                        L = Pk1'*Pk2;
%                        F_temp = det(L) / abs(det(L));
%                        L = Pk2'*Pk3;
%                        F_temp = F_temp*det(L) / abs(det(L));
%                        L = Pk3'*Pk4;
%                        F_temp = F_temp*( det(L) / abs(det(L)) );
%                        L = Pk4'*Pk1;
%                        F_temp = F_temp*( det(L) / abs(det(L)) );
%                    %-----------------------------
              end %endif
          end %endif
      end %endif
      F = F + angle(F_temp)/(2*pi);
    end
end
F;
A;
Z2 = mod(F-A,2)

    if( abs(Z2-1) <= 1e-1 )
        disp('QSH phase');    
    end

    if( abs(Z2) <= 1E-1)
        disp('Normal phase');
    end
    
end
%Now you have calculated the Z2 index
%--------------------------------------------------------------------------
%            for m=1:4
%                for n=1:4
%                    %w_mn = < um(-k)|T|un(k) >
%                    w(m,n) = dot(psi_k(:,m),T*conj( psik(:,n) ));
%                end
%            end
%            Pf_w = w(1,2)*w(3,4) - w(1,3)*w(2,4) + w(2,3)*w(1,4);
%--------------------------------------------------------------------------
% Calculating berry curvature F for given k
%                        for m=1:4
%                            for n=1:4
%                                %L_mn = < um(k_i)|un(k_i+1) >
%                                if( abs( conj(Pk1(:,n)')*G2*Pk1(:,n)+conj(Pk1(:,m)')*G2*Pk1(:,m))==0 )
%                                    L(m,n) = dot(Pk1(:,n),Pk2(:,m));
%                                    F_temp = det(L) / abs(det(L));
% 
%                                    L(m,n) = dot(Pk2(:,n),Pk3(:,m));
%                                    F_temp = F_temp*det(L) / abs(det(L));
% 
%                                    L(m,n) = dot(Pk3(:,n),Pk4(:,m));
%                                    F_temp = F_temp*det(L) / abs(det(L));
% 
%                                    L(m,n) = dot(Pk4(:,n),Pk1(:,m));
%                                    F_temp = F_temp*det(L) / abs(det(L));
%                                end
%                            end
%                        end
%--------------------------------------------------------------------------


