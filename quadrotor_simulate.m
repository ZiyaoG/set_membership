clear all;
close all;
%% This is a set membership test on quadrotor system with a geometric controller. 
dt=0.002; % 500 Hz
T=10; % Simulation interval
flag.traj=1; % 1 for sin traj, 2 for flip
outputFolder='datas\quadrotors\tg';

% Quadrotor stats
J=diag([0.082 0.0845 0.1377]);
m=4.34;
d=0.315;
c_tauf=8.004*10^(-4);
g=9.8;
e3=[0;0;1];
% controller parameters
kx=16;
kv=5.6;
kR=8.81;
kw=2.54;

% States Innitiation
x=zeros([3,T/dt]);
v=zeros([3,T/dt]);
R=zeros([3,3,T/dt]);
w=zeros([3,T/dt]);
xd=zeros([3,T/dt]);

f=zeros([1,T/dt]);
M=zeros([3,T/dt]);
M_expl=0.5*normrnd(0,1,[3,T/dt]);
phi=zeros([6,T/dt]);


% Innitial condition
x(:,1)=[0;0;0];
v(:,1)=[0;0;0];
if flag.traj ==1 
    R(:,:,1)=eye(3);
elseif flag.traj == 2
    R(:,:,1)=[1 0 0;0 -0.995 -0.0314;0 0.0314 -0.9995];
else
    error("please choose a task")
end
w(:,1)=[0;0;0];

% Gaussian noise
pd = makedist('Normal');
tru_gp = truncate(pd,-1,1);

% uniform noise
% f_dist = -0.2+0.4*rand(T/dt,3);
% M_dist = -0.2+0.4*rand(T/dt,3);

for j=1:99
    fprintf('Running simulation "%2.2d".\n', j);
    f_dist = random(tru_gp,T/dt,3);
    M_dist = random(tru_gp,T/dt,3);
    for i=1:(T/dt-1)
        % at time t
        t=i*dt;
        xt=x(:,i);
        vt=v(:,i);
        Rt=R(:,:,i);
        wt=w(:,i);
        b3=Rt*e3;
    
        % error systems
        [djd,jd,ad,vd,xdt,b1,db1,ddb1]=desired_traj(t,flag.traj);
        xd(:,i)=xdt;
        ex=xt-xdt;
        ev=vt-vd;
    
        A=-kx*ex-kv*ev-m*g*e3+m*ad;
        ft=-dot(A,b3);
        ea=g*e3-ft/m*b3-ad;
    
        dA=-kx*ev-kv*ea+m*jd;
        db3=Rt*wedge(wt)*e3;
        dft=-dot(dA,b3)-dot(A,db3);
        ej=-dft/m*b3-ft/m*db3-jd;
        ddA=-kx*ea-kv*ej+m*djd;
        [b3d,db3d,ddb3d]=deriv_unit_vector(-A,-dA,-ddA);
    
        
        A2 = -wedge(b1) * b3d;
        dA2 = -wedge(db1) * b3d - wedge(b1) * db3d;
        ddA2 = - wedge(ddb1) * b3d - 2 * wedge(db1) * db3d - wedge(b1) * ddb3d;
        
        [b2d, db2d, ddb2d] = deriv_unit_vector(A2, dA2, ddA2);
        
        b1d = wedge(b2d) * b3d;
        db1d = wedge(db2d) * b3d + wedge(b2d) * db3d;
        ddb1d = wedge(ddb2d) * b3d + 2 * wedge(db2d) * db3d + wedge(b2d) * ddb3d;
        
        Rd = [b1d, b2d, b3d];
        dRd = [db1d, db2d, db3d];
        ddRd = [ddb1d, ddb2d, ddb3d];
    
        wd = vee(Rd' * dRd);
        dwd = vee(Rd' * ddRd - wedge(wd)^2);
        ew=wt-Rt'*Rd*wd;
        eR=0.5*vee(Rd'*Rt-Rt'*Rd);
    
        % control input total thrust and moment
        
        Mt=-kR*eR-kw*ew+cross(wt,J*wt)-J*(wedge(wt)*Rt'*Rd*wd-Rt'*Rd*dwd);
        % Mt=-kR*eR-kw*ew+wedge(Rt'*Rd*wd)*J*Rt'*Rd*wd+J*Rt'*Rd*dwd;
    
    
    
        % dynamics update
        % f_dist_t=0;
        % M_dist_t=[0;0;0];
        % f_dist_t=0.5*[sin(5*t);cos(3*t);cos(7*t)];
        f_dist_t=0.5*f_dist(i,:)';
        % M_dist_t=0.5*[1;1;1];
        % M_dist_t=0.5*[sin(t);cos(2*t);cos(5*t)];
        M_dist_t=0.5*M_dist(i,:)';
        M_expl_t=M_expl(:,i);
        [x_next,v_next,R_next,w_next]=quadrotor_dyn(xt,vt,Rt,wt,ft,Mt+M_expl_t,f_dist_t,M_dist_t,dt,m,g,J,e3);
        
    
        % logging
        f(:,i)=ft;
        M(:,i)=Mt;
        x(:,i+1)=x_next;
        v(:,i+1)=v_next;
        R(:,:,i+1)=R_next;
        w(:,i+1)=w_next;
        phi(1,i)=w_next(1,1)-wt(1,1);
        phi(2,i)=w_next(2,1)-wt(2,1);
        phi(3,i)=w_next(3,1)-wt(3,1);
        phi(4,i)=dt*wt(2,1)*wt(3,1);
        phi(5,i)=dt*wt(1,1)*wt(3,1);
        phi(6,i)=dt*wt(1,1)*wt(2,1);

        % for debug
        tmp=abs(m/dt*(v(:,i+1)-v(:,i))-m*g*e3+ft*Rt*e3);
    
    end
    baseFileName = sprintf('Trajecory_%2.2d_Task1_tg_f05_M05.mat', j);
    fullFileName = fullfile(outputFolder, baseFileName);
    save(fullFileName,'M','M_expl','M_dist','w','dt','phi');
end

%% plot
% figure('Name','position')
% plot3(x(1,:),x(2,:),x(3,:))
% grid on
% hold on
% plot3(xd(1,:),xd(2,:),xd(3,:))
% legend('x','xd');
% axis equal
% 
% figure('Name','angles')
% subplot(3,1,1)
% plot(M(1,:))
% subplot(3,1,2)
% plot(M(2,:))
% subplot(3,1,3)
% plot(M(3,:))


% save('datas\quadrotors\tg\Trajecory_Task1_tg_f05_M05.mat');
% save('datas\quadrotors\Trajecory_Task1_test.mat'); % sin cos noise
% save('datas\quadrotors\Trajecory_augdyn_Task1_uniform_f02_M02.mat');
% save('datas\quadrotors\Trajecory_augdyn_Task1_nonoise.mat');


% 
% fname = sprintf('figs/quadrotors/Trajecory_Task%0.3f.png', flag.traj);
% exportgraphics(gcf,fname)





