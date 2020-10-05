a1=0.5*[-sqrt(3), 3];
a2=0.5*[sqrt(3), 3];
Lso=0.2;
Lr=0.0;
t=1.0;
V=0.0;
Lambda=0;
%plik=fopen("coœ.txt", "w");
step=pi/10;
Cn1=0; Cn2=0; Cn3=0; Cn4=0;
Kxx=[0,2*pi*sqrt(3)/3,0];
Kyy=[0,2*pi/3,4*pi/3];
kx=-pi;
for kx=-pi:step:pi
  ky=-pi;
  k=[kx,ky];
  k1=dot(k,a1);
  k2=dot(k,a2);
  H=zeros(4,4);
  H(1,1)=-2*Lso*(sin(k1)-sin(k2)+sin(k2-k1))+V;
  H(2,2)=-H(1,1)+2*V;
  H(3,3)=-H(1,1);
  H(4,4)=H(1,1)-2*V;
  H(3,1)=t*(1+exp(i*k1)+exp(i*k2));
  H(4,2)=H(3,1);
  H(3,2)=i*Lr*((0.5+i*sqrt(3)/2)*exp(i*k2)+(0.5-i*sqrt(3)/2)*exp(i*k1)-1);
  H(4,1)=i*Lr*((0.5-i*sqrt(3)/2)*exp(i*k2)+(0.5+i*sqrt(3)/2)*exp(i*k1)-1);
  H=H+H';
  [V1,e]=eig(H);
  
  if(V1(1,1)==0)
    temp = V1(:,1);
    V1(:,1) = V1(:,2);
    V1(:,2) = temp;
  end
  if(V1(1,4)==0)
    temp = V1(:,4);
    V1(:,4) = V1(:,3);
    V1(:,3) = temp;
  end
  
  k=[kx+step,ky];
  k1=dot(k,a1);
  k2=dot(k,a2);
  H=zeros(4,4);
  H(1,1)=-2*Lso*(sin(k1)-sin(k2)+sin(k2-k1))+V;
  H(2,2)=-H(1,1)+2*V;
  H(3,3)=-H(1,1);
  H(4,4)=H(1,1)-2*V;
  H(3,1)=t*(1+exp(i*k1)+exp(i*k2));
  H(4,2)=H(3,1);
  H(3,2)=i*Lr*((0.5+i*sqrt(3)/2)*exp(i*k2)+(0.5-i*sqrt(3)/2)*exp(i*k1)-1);
  H(4,1)=i*Lr*((0.5-i*sqrt(3)/2)*exp(i*k2)+(0.5+i*sqrt(3)/2)*exp(i*k1)-1);
  H=H+H';
  [V2,e]=eig(H);
  if(V2(1,1)==0)
    temp = V2(:,1);
    V2(:,1) = V2(:,2);
    V2(:,2) = temp;
  end
  if(V2(1,4)==0)
    temp = V2(:,4);
    V2(:,4) = V2(:,3);
    V2(:,3) = temp;
  end
  
  for ky=-pi:step:pi
    k=[kx+step,ky+step];
    k1=dot(k,a1);
    k2=dot(k,a2);
    H=zeros(4,4);
    H(1 ,1)=-2*Lso*(sin(k1)-sin(k2)+sin(k2-k1))+V;
    H(2,2)=-H(1,1)+2*V;
    H(3,3)=-H(1,1);
    H(4,4)=H(1,1)-2*V;
    H(3,1)=t*(1+exp(i*k1)+exp(i*k2));
    H(4,2)=H(3,1);
    H(3,2)=i*Lr*((0.5+i*sqrt(3)/2)*exp(i*k2)+(0.5-i*sqrt(3)/2)*exp(i*k1)-1);
    H(4,1)=i*Lr*((0.5-i*sqrt(3)/2)*exp(i*k2)+(0.5+i*sqrt(3)/2)*exp(i*k1)-1);
    H=H+H';
    [V3,e]=eig(H);
    if(V3(1,1)==0)
    temp = V3(:,1);
    V3(:,1) = V3(:,2);
    V3(:,2) = temp;
    end
    if(V3(1,4)==0)
      temp = V3(:,4);
      V3(:,4) = V3(:,3);
      V3(:,3) = temp;
    end
    
    k=[kx,ky+step];
    k1=dot(k,a1);
    k2=dot(k,a2);
    H=zeros(4,4);
    H(1,1)=-2*Lso*(sin(k1)-sin(k2)+sin(k2-k1))+V;
    H(2,2)=-H(1,1)+2*V;
    H(3,3)=-H(1,1);
    H(4,4)=H(1,1)-2*V;
    H(3,1)=t*(1+exp(i*k1)+exp(i*k2));
    H(4,2)=H(3,1);
    H(3,2)=i*Lr*((0.5+i*sqrt(3)/2)*exp(i*k2)+(0.5-i*sqrt(3)/2)*exp(i*k1)-1);
    H(4,1)=i*Lr*((0.5-i*sqrt(3)/2)*exp(i*k2)+(0.5+i*sqrt(3)/2)*exp(i*k1)-1);
    H=H+H';
    [V4,e]=eig(H);
    if(V4(1,1)==0)
    temp = V4(:,1);
    V4(:,1) = V4(:,2);
    V4(:,2) = temp;
    end
    if(V4(1,4)==0)
      temp = V4(:,4);
      V4(:,4) = V4(:,3);
      V4(:,3) = temp;
    end
    
    x=dot(V1(:,1),V2(:,1))/abs(dot(V1(:,1),V2(:,1)));
    x=x*dot(V2(:,1),V3(:,1))/abs(dot(V2(:,1),V3(:,1)));
    x=x*dot(V3(:,1),V4(:,1))/abs(dot(V3(:,1),V4(:,1)));
    x=x*dot(V4(:,1),V1(:,1))/abs(dot(V4(:,1),V1(:,1)));
    dS=step^2;
    krzywizna1 = angle(x)/dS; 
      
    x=dot(V1(:,2),V2(:,2))/abs(dot(V1(:,2),V2(:,2)));
    x=x*dot(V2(:,2),V3(:,2))/abs(dot(V2(:,2),V3(:,2)));
    x=x*dot(V3(:,2),V4(:,2))/abs(dot(V3(:,2),V4(:,2)));
    x=x*dot(V4(:,2),V1(:,2))/abs(dot(V4(:,2),V1(:,2)));
    dS=step^2;
    krzywizna2 = angle(x)/dS; 
    
    x=dot(V1(:,3),V2(:,3))/abs(dot(V1(:,3),V2(:,3)));
    x=x*dot(V2(:,3),V3(:,3))/abs(dot(V2(:,3),V3(:,3)));
    x=x*dot(V3(:,3),V4(:,3))/abs(dot(V3(:,3),V4(:,3)));
    x=x*dot(V4(:,3),V1(:,3))/abs(dot(V4(:,3),V1(:,3)));
    dS=step^2;
    krzywizna3 = angle(x)/dS; 
    
    x=dot(V1(:,4),V2(:,4))/abs(dot(V1(:,4),V2(:,4)));
    x=x*dot(V2(:,4),V3(:,4))/abs(dot(V2(:,4),V3(:,4)));
    x=x*dot(V3(:,4),V4(:,4))/abs(dot(V3(:,4),V4(:,4)));
    x=x*dot(V4(:,4),V1(:,4))/abs(dot(V4(:,4),V1(:,4)));
    dS=step^2;
    krzywizna4 = angle(x)/dS; 
    
    if(inpolygon(kx,ky,Kxx,Kyy)==true)
      Cn1=Cn1+krzywizna1*dS/(2*pi);
      Cn2=Cn2+krzywizna2*dS/(2*pi);
      Cn3=Cn3+krzywizna3*dS/(2*pi);
      Cn4=Cn4+krzywizna4*dS/(2*pi);
    end
    %fprintf(plik, "%f\t %f\t %f\n", kx, ky,krzywizna);
    V1 = V4; V2 = V3;
    end
end
Cherndown=Cn1+Cn3
Chernup=Cn2+Cn4
Chernspin=(Chernup-Cherndown)/2
%fclose(plik)