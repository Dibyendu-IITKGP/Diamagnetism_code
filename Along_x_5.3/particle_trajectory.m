clear;
close all;
clc;

load('field.mat');

%Time step
dt=0.4;
%Number of steps
step=50000;

%Initial position
x_i=5.3/1000;
y_i=10^-6;
z_i=zb+0.0002;

%Plot variables
xp=zeros(step+1,1);
yp=zeros(step+1,1);
zp=zeros(step+1,1);

xp(1,1)=x_i;
yp(1,1)=y_i;
zp(1,1)=z_i;


for it=1:step
               
        x=x_i; y=y_i; z=z_i;
        
        B1=double(subs(Bxx));
        B2=double(subs(Byy));
        B3=double(subs(Bzz));
        
        
        %Updating the position of the particle
        
        x_f=x_i+(ki_p-ki_f)*m_p*B1*(1/(Ro_p*6*pi*r*Ro_f*nu_f*mu_0))...
            *(dt+(m_p+0.5*m_f)*(1/(6*pi*r*Ro_f*nu_f))*exp(-dt*6*pi*r*Ro_f*nu_f*(1/(m_p+0.5*m_f))));

        y_f=y_i+(ki_p-ki_f)*m_p*B2*(1/(Ro_p*6*pi*r*Ro_f*nu_f*mu_0))...
            *(dt+(m_p+0.5*m_f)*(1/(6*pi*r*Ro_f*nu_f))*exp(-dt*6*pi*r*Ro_f*nu_f*(1/(m_p+0.5*m_f))));

        z_f=z_i+(((ki_p-ki_f)*m_p*B3*(1/(Ro_p*6*pi*r*Ro_f*nu_f*mu_0)))-((m_p-Ro_f*V_p)*g*(1/(6*pi*r*Ro_f*nu_f))))...
            *(dt+(m_p+0.5*m_f)*(1/(6*pi*r*Ro_f*nu_f))*exp(-dt*6*pi*r*Ro_f*nu_f*(1/(m_p+0.5*m_f))));

        x_i=x_f;
        y_i=y_f;
        if z_f > zb+0.0002 
            z_i=z_f;
        end
        
        if abs(x_i)<=0.0001 && abs(y_i)<=0.0001 && abs(z_i)<=(zb+0.0003)
            break
        end
            
        
    
        xp(it+1,1)=x_i;
        yp(it+1,1)=y_i;
        zp(it+1,1)=z_i;
        
        

end
    
save('trajectory.mat');

plot3(xp,yp,zp);
        
