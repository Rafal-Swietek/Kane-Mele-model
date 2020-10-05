function Z2 = Z2_TRIM(parameters,M,a1,a2)
% This function allows to determine the Z2 index for the Kane-Mele model
% using the sewing matrix at TRIM points
%
% This method is only applicable if inversion symmetry is preserved, thus
% only for:
%                               V=0!!!
%--------------------------------------------------------------
%Initial parameters
    t = parameters(1); %Nearest Neighbours
    V = parameters(2); %staggered potential 
    LSO = parameters(3); %Next Nearest Neigbours
    LR = parameters(4)*0; %Rashba term
%--------------------------------------------------------------    
%Matrices 
    %Pauli matrices   
        sig_x = [0 1;1 0];
        sig_y = [0 -i;i 0];
        sig_z = [1 0;0 -1];
        I = [1 0;0 1];
        
        T = i*kron(I,sig_y); %Time-reversal operator
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
 T = i*kron(sig_y,I); %Time-reversal operator  
 P = kron(I,sig_x); %Parity operator
    
%Calculating Z2 using TRIM points: all M points and Gamma   
w = zeros(4,4);
X = 1;
    for ii = 1:7
        if(ii<=6)
            kx = M(ii,1);
            ky = M(ii,2);
        else
            kx = 0;
            ky = 0;
        end
         %Eigenvectors for k
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
           H = H + H_R + H_R';
           [psik,~] = eig(H);
           psik(2,:) = psik(2,:)*exp(i*x);
           psik(4,:) = psik(4,:)*exp(i*x);
            %Properly reordering the bands
                if(abs( psik(1,1) )==0)
                   temp = psik(:,1);
                   psik(:,1) = psik(:,2);
                   psik(:,2) = temp;
                end
                if(abs( psik(1,4) )==0)
                   temp = psik(:,3);
                   psik(:,3) = psik(:,4);
                   psik(:,4) = temp;
                end
          %Calculating matrix
%               w(1,2) = dot(psik(:,1),T*conj( psik(:,2) ));
%               w(2,1) = dot(psik(:,2),T*conj( psik(:,1) ));
%               Pf_w = 0.5*( w(1,2) - w(2,1));
           for m=1:4
               for n=1:4
                   %w_mn = < um(-k)|T|un(k) >
                   if( abs( conj(psik(:,n)')*G2*psik(:,n)+conj(psik(:,m)')*G2*psik(:,m))==0 )
                        w(m,n) = dot(psik(:,m),T*conj( psik(:,n) ));
                   end
               end
           end
           Pf_w = w(1,2)*w(3,4) - w(1,3)*w(2,4) + w(2,3)*w(1,4);
           X = X*sqrt( det(w) )/Pf_w;
%            %Using the parity operator---------
%                Xi = Xi*dot( (psik(:,1)),P*psik(:,1) );
%                Xi = Xi*dot( (psik(:,2)),P*psik(:,2) );
%                Xi = Xi*dot( (psik(:,3)),P*psik(:,3) );
%                Xi = Xi*dot( (psik(:,4)),P*psik(:,4) );
%            %----------------------------------
    end
    % Now we have (-1)^Z2 = X, thus:
    % calculating Z2 invariant
    %Z2 = real(X);
    if( abs(X+1) <= 1e-5 )
        Z2 = 1
        disp('QSH phase');    
    elseif( abs(X-1) <= 1e-5 )
        Z2 = 0
        disp('Normal phase');
    else
        Z2 = X
    end
end

%Now you have calculated the Z2 index
%---------------------------------------------
%            for m=1:4
%                for n=1:4
%                    %w_mn = < um(-k)|T|un(k) >
%                    w(m,n) = dot(psi_k(:,m),T*conj( psik(:,n) ));
%                end
%            end
%            Pf_w = w(1,2)*w(3,4) - w(1,3)*w(2,4) + w(2,3)*w(1,4);



