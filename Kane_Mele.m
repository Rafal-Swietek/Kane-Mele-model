clear all
clc

%This program was made to solve the Kane - Mele problem on the Honeycomb
%lattice of graphene. Varoius parameters are included.

%K-space base vectors
    a1 = 0.5*[-3^(1/2),3]; %Bravais lattice vectors
    a2 = 0.5*[3^(1/2),3];
%     %Solving for K-space vectors    
%         syms b_1x b_1y b_2x b_2y; %variables
%         equations1 = [b_1x*a1(1) + b_1y*a1(2)==2*pi, b_1x*a2(1) + b_1y*a2(2)==0   ];
%         equations2 = [b_2x*a1(1) + b_2y*a1(2)== 0  , b_2x*a2(1) + b_2y*a2(2)==2*pi];
%         sol =   solve(equations1, [b_1x b_1y]); %solving for b1
%         solut = solve(equations2, [b_2x b_2y]); %solving for b2
%         
%         b1 = [sol.b_1x sol.b_1y]; %assigning solution to vector
%         b2 = [solut.b_2x solut.b_2y];%-||-
b1 = 2*pi*[-1/sqrt(3),1/3];
b2 = 2*pi*[1/sqrt(3),1/3];
%%
% Finding K & M points
    KK = [b1; b1+b2; b2; -b1; -b1-b2; -b2; b1]; %Unit cell
    r = norm(b1)/3^0.5; %radius of circle on Wigner-Seitz cell
    alfa = pi/3;  
    Rot = [cos(alfa) -sin(alfa); sin(alfa) cos(alfa)]; %matrix rotating around 60 degrees (pi/3)
    K1 = [r 0]'; %oncircle of hexagon
    K = [ K1'; (Rot*K1)'; (Rot*(Rot*K1))'; -K1'; -(Rot*K1)'; -(Rot*(Rot*K1))'; K1']; %high symmetry K-points
    M1 = b2'/2; %half lattice vector
    M = [M1'; (Rot*M1)'; (Rot*(Rot*M1))'; -M1'; -(Rot*M1)'; -(Rot*(Rot*M1))'; M1']; %all M points
%%
%Energy eigenvalue calculation

        t = 1.; %hopping integral -> Nearest Neighbours
        V = 0.0; %staggered potential -> breaks inversion symmetry
        LSO = 0.0; %spin-orbit coupling -> Next Nearest Neigbours
        LR = 0.0; %Rashba term -> perpendicular electric field or interaction with substrate
 
%%       
%Defining k-path & Energy along k-path
grid = 100; %grid of k-space (k_path)
parameters = [t, V, LSO, LR]; %array of parameters as input for the functions

%   k_path = K_gridGMKG(b1,b2,K,M,grid); %creates a G-K-M-K'-G path in k-space alog high symmetry points
%   Energy_k_path(k_path,parameters,a1,a2); %calculates the energy along the given path in k-space


% EBZ = Effective_Brillouin_Zone(b1,b2,M,K,grid); %creates the Effective Brillouin Zone
% if(V == 0)  %determines the Z2 index in topological insulators
%     Z2 = Z2_TRIM(parameters,M,a1,a2); %non-broken inversion symmetry
% else
%     Z2 = Z2_invariant(EBZ,parameters,a1,a2,M,grid); %broken inversion symmetry
% end
% hold off
% grid = 100;
% Spin_Chern(parameters,a1,a2,grid,M);

%Energy3D(parameters,a1,a2,grid,M); %energy output in 3D: solving on kx, ky plane

%--------------------------------------------------------------------------

EBZ = Effective_Brillouin_Zone(b1,b2,M,K,grid); %creates the Effective Brillouin Zone
parameters = [-1, 0.0, 0.06, 0]; %V = LR =0
Z2 = Z2_TRIM(parameters,M,a1,a2); %Normal phase

parameters = [-1, 0.0, 0.15, 0.05]; %V = 0
Z2 = Z2_TRIM(parameters,M,a1,a2); %?

parameters = [-1, 0.1, 0.0, 0.0]; % LSO = LR = 0
Z2 = Z2_invariant(EBZ,parameters,a1,a2,M,grid); %Normal phase

parameters = [-1, 0.1, 0.0, 0.05]; %LSO = 0
Z2 = Z2_invariant(EBZ,parameters,a1,a2,M,grid); %Normal phase

%--------------------------------------------------------------------------



    
    
    
    