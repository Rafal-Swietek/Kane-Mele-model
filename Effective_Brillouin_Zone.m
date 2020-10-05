function EBZ = Effective_Brillouin_Zone(b1,b2,MM,KK,discret)
    M1 = MM(1,:); %M1 - point
    M2 = MM(2,:); %M2 - point
    K = KK(1,:); %K-point
    G1 = 2*M1; %Gamma'
    G2 = 2*M2; %Gamma"
    EBZ = zeros(2,3*discret); %Effective Brillouin Zone points
    for ii = 0:discret-1
        EBZ(:,ii+1) = [G1(1)/discret*ii; G1(2)/discret*ii]; %Gamma->M1->Gamma' along x-axis
        EBZ(:,ii+1+discret) = [G1(1)-(G1(1)-G2(1))/discret*ii; G1(2)+(G2(2)-G1(2))/discret*ii]; %Gamma'->Gamma"
        EBZ(:,ii+1+2*discret) = [G2(1)-G2(1)/discret*ii; G2(2)-G2(2)/discret*ii];%Gamma"->M2->Gamma
    end

end

%   %-------------------------------------------------- 
%   %Plotting k-space & reduced BZ---------------------
%     theta = 0:60:360; fi = 0:60:360;
%     x1 = norm(K)*cosd(theta);  y1 = norm(K)*sind(theta);
%     x2 = KK(1,1)*cosd(fi);   y2 = KK(1,1)*sind(fi);
%     
%     figure(99);
%     %Extended Brillouin Zone----
%         p1 = plot(x1, y1, 'bo-','MarkerSize',5); hold on                      
%         plot(x1+G1(1), y1+G1(2), 'bo-','MarkerSize',5); hold on
%         plot(x1+G2(1), y1+G2(2), 'bo-','MarkerSize',5); hold on
%     %---------------------------
%     p2 = plot(EBZ(1,:),EBZ(2,:),'k-o','MarkerSize',1); hold on  
%         xticks([-pi -pi/2 0 pi/2 pi 1.5*pi 2*pi]);   
%         xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi','3/2\pi','2\pi'});
%         yticks([-pi -pi/2 0 pi/2 pi 1.5*pi 2*pi]);   
%         yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi','3/2\pi','2\pi'});
%     text(-0.4,0.0,'\Gamma','fontsize',14); %labeling points with Gamma
%     text(G1(1)+0.4,G1(2),'\Gamma"','fontsize',14); 
%     text(G2(1),G2(2)+0.4,'\Gamma""','fontsize',14); 
%     axis([-3.5 6.5 -3.5 6.5]);
%     legend([p1 p2],'k-space lattice','Extended Brillouin Zone');
%   %-------------------------------------------------- 