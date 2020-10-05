function k_track = K_gridGMKG(b1,b2,KK,MM,discret)
    K = KK(1,:); %K - point
    K2 = KK(2,:);
    k_track = zeros(2,3*discret);
    for ii = 0:discret-1
        k_track(:,ii+1) = [K(1)/discret*ii; K(2)/discret*ii]; %Gamma->K along x-axis
        k_track(:,ii+1+discret) = [K(1)-(K(1)-K2(1))/discret*ii; K(2)+(K2(2)-K(2))/discret*ii]; %K->K'
        k_track(:,ii+1+2*discret) = [K2(1)-K2(1)/discret*ii; K2(2)-K2(2)/discret*ii];%K'->Gamma
    end
    
%     %Plotting k-space & reduced BZ
%     theta = 30:60:390; fi = 0:60:360;
%     x1 = norm(b1)*cosd(theta);  y1 = norm(b1)*sind(theta);
%     x2 = KK(1,1)*cosd(fi);   y2 = KK(1,1)*sind(fi);
%     
%     figure(99);
%     p1 = plot(x1, y1, 'bo-','MarkerSize',5); hold on                      %Unit cell hexagon
%     p2 = plot(x2, y2, 'gs-','MarkerSize',4); hold on                      %Wigner-Seitz cell hexagon
%     plot(MM(:,1),MM(:,2),'ko','MarkerSize',4); hold on                    %M-point
%     p4 = plot(k_track(1,:),k_track(2,:),'k-o','MarkerSize',1); hold on
%     p5 = quiver([0 0],[0 0],1.1*[b1(1) b2(1)],1.1*[b1(2) b2(2)],'r-');    %base vector   
%         xticks([-1.5*pi -pi -pi/2 0 pi/2 pi 1.5*pi]);   
%         xticklabels({'-3/2\pi','-\pi','-\pi/2','0','\pi/2','\pi','3/2\pi'});
%         yticks([-1.5*pi -pi -pi/2 0 pi/2 pi 1.5*pi]);   
%         yticklabels({'-3/2\pi','-\pi','-\pi/2','0','\pi/2','\pi','3/2\pi'});
%     
%     A = [KK(1,1) KK(3,1) KK(5,1)];
%     B = [KK(1,2)-0.5 KK(3,2) KK(5,2)];
%     C = [KK(2,1) KK(4,1)-0.3 KK(6,1)];
%     D = [KK(2,2) KK(4,2) KK(6,2)];
%     text(A + 0.05*A/norm(A),B + 0.5*B/norm(B),'K','fontsize',14); %labeling points with K
%     text(C + 0.05*C/norm(C),D + 0.5*D/norm(D),'K"','fontsize',14); %labeling points with K'
%         AA = [MM(1,1)-0.15,MM(2,1),MM(3,1)-0.15,MM(4,1)+0.1,MM(5,1),MM(6,1)+0.1];
%         BB = [MM(1,2)+0.4,MM(2,2)+0.4,MM(3,2)-0.4,MM(4,2)+0.1,MM(5,2)-0.3,MM(6,2)-0.1];
%     text(AA,BB,'M','fontsize',14); %labeling points with K'
%     %text(M(1)+0.2,M(2)-0.15,'M','fontsize',14); %labeling points with M
%     text(-0.4,0.0,'\Gamma','fontsize',14); %labeling points with M
%     axis([-6 6 -6 6]);
%     legend([p1 p2 p4 p5],'Unit Cell','Wigner-Seitz cell','Reduced Brillouin Zone','Base vectors');
end

