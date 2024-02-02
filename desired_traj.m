function [djd,jd,ad,vd,xd,b1,db1,ddb1]=desired_traj(t,flag)
    % need forth derivative of position x
    if flag ==1
        xd=[0.4*t;0.4*sin(pi*t);0.6*cos(pi*t)];
        b1=[cos(pi*t);sin(pi*t);0];
        vd=[0.4;0.4*pi*cos(pi*t);0.6*pi*(-sin(pi*t))];
        ad=[0;0.4*pi^2*(-sin(pi*t));-0.6*pi^2*cos(pi*t)];
        jd=[0;-0.4*pi^3*cos(pi*t);0.6*pi^3*sin(pi*t)];
        djd=[0;0.4*pi^4*sin(pi*t);0.6*pi^4*cos(pi*t)];
        db1=[-pi*sin(pi*t);pi*cos(pi*t);0];
        ddb1=[-pi^2*cos(pi*t);-pi^2*sin(pi*t);0];
    
    elseif flag==2
    % flip
        xd=[0;0;0];
        b1=[1;0;0];
        vd=zeros(3,1);
        ad=zeros(3,1);
        jd=zeros(3,1);
        djd=zeros(3,1);
        db1=zeros(3,1);
        ddb1=zeros(3,1);
    end
    


end