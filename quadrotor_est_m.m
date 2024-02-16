clear all;
close all;
load("datas\quadrotors\Trajecory_Task1_constMdist05.mat");
% % Analysis the datas from geometric controlled quadrototrs and use 
% % set membership to give estimation of the m and J. 
% clear all;
% close all;
% load("datas\quadrotors\Trajecory_Task1_constMdist05.mat");
% m_est.value=4:0.002:5;
% msize=size(m_est.value);
% m_est.flag=ones(msize);
% Bound=2.5;
% 
% for i=1:msize(2)
%     for j=1:(T/dt-1)
%         m_tmp=abs(m_est.value(i)/dt*(v(:,j+1)-v(:,j))-m_est.value(i)*g*e3+f(j)*R(:,:,j)*e3);
%         if sum(double(m_tmp>Bound))>0
%             m_est.flag(i)=0;
%             break;
%         end
%     end
% end
% 
% %% plot
% plot(m_est.value(1,:),m_est.flag(1,:),'o')
% hold on
% plot(m,1,'*')
% legend('estimation','Ground Truth')
% 
% % fname = sprintf('figs/quadrotors/Task%0.3f_Bound%0.3f.png', flag,Bound);
% % exportgraphics(gcf,fname)
%% estimate the m_max_bar and m_max_est
m_est=sdpvar(1,1);
Constraints= -m_est<=0;
data_range=10:2:100;
i=1000;
step=7;
tmp=[];
nx=3;
nz=3;
O1=3; % parameter in tilde O notation
% O2=7; % parameter in tilde O notation
m_max_bar=[];
m_max_est=[];
for data_number=data_range
    for n=1:data_number
        j=step*n;
        tmp=[tmp,max(abs(m_est/dt * (v(:,i+j+1) - v(:,i+j)) - m_est*g*e3+f(1,i+j)*R(:,:,i+j)*e3))];
    end
    Objective= max(tmp);
    options = sdpsettings('verbose', 0, 'solver', 'mosek'); % Choose a solver that is available and appropriate
    sol = optimize(Constraints, Objective, options);
    if sol.problem == 0
            m_est_value = double(m_est);
            % disp(m_est_value)
            obj_value=double(Objective);
            m_max_bar=[m_max_bar,obj_value];
            % disp(obj_value);
            g_T=nx^1.5*nz^2/data_number;
            delta_T=max(0,O1*g_T*log(g_T));
            % disp(delta_T);
            m_max_est=[m_max_est,obj_value+delta_T];
            % disp(m_max_est)
    else
            disp('No solution found');
            return;
    end
    
end
% % check the convergence of delta_T
% delta_T_series=[];
% for i=10:100
%     g_T=nx^1.5*nz^2/i;
%     delta_T_series=[delta_T_series,O1*g_T*log(g_T)]
% end
% plot(delta_T_series)
%% plot
plot(data_range,m_max_bar)
hold on
plot(data_range,m_max_est)
legend('$\bar{m}_{max}$','$\hat{m}_{max}$','Interpreter','latex')
xlabel('number of datapoints')
ylabel('$m_{max}$','Interpreter','latex')
% ylim([0,1])

fname = sprintf('figs/quadrotors/convergenceofm_max.png');
exportgraphics(gcf,fname)
return;
%% Another way to do this

loose=2;   % 1 for tight bound, >1 for conservative bound
B=loose*[0.5;0.5;0.5];
i=1000;  % first datapoint index to start with
data_number=20;   % number of datapoints used to generate the sets
step=7;  % step between two adjacent datapoints 
Diameter_m=zeros([data_number,1]);
m_minmax=zeros([data_number,2]);
m1=sdpvar(1,1);
m2=sdpvar(1,1);
Constraints_m1=-m1<=0;
Constraints_m2=-m2<=0;
Constraints=[Constraints_m1,Constraints_m2];
for n=1:data_number
    for j=step*(1:n)
        LMI_m1=abs(m1/dt * (v(:,i+j+1) - v(:,i+j)) - m1*g*e3+f(1,i+j)*R(:,:,i+j)*e3)-B<=0;
        LMI_m2=abs(m2/dt * (v(:,i+j+1) - v(:,i+j)) - m2*g*e3+f(1,i+j)*R(:,:,i+j)*e3)-B<=0;
        % m_tmp=abs(m_est.value(i)/dt*(v(:,j+1)-v(:,j))-m_est.value(i)*g*e3+f(j)*R(:,:,j)*e3);
        Constraints=[Constraints,LMI_m1,LMI_m2];
    end
    Distance=m1-m2;
    Objective = -abs(Distance);
    options = sdpsettings('verbose', 0, 'solver', 'mosek'); % Choose a solver that is available and appropriate
    sol = optimize(Constraints, Objective, options);
    
    if sol.problem == 0
        m1_value = double(m1);
        m2_value = double(m2);
        % disp('J1:');
        % disp(J1_value);
        % disp('J2:');
        % disp(J2_value);
        disp('--------------------------------------------')
        fprintf('Diameter at %.1f: max|m1-m2|=%f\n',n,abs(m1_value-m2_value))
        % disp(sum(abs(J1_value-J2_value),'all'));
    else
        disp('No solution found');
        return;
    end
    Diameter_m(n,1)=abs(m1_value-m2_value);
    m_minmax(n,1)=m1_value;
    m_minmax(n,2)=m2_value;
end

%%
% plot(dt*(i+(1:data_number)),Diameter_J)
plot(i+(1:data_number),m_minmax(:,1))
hold on
plot(i+(1:data_number),m_minmax(:,2))
plot(i+(1:data_number),m*ones(size(1:data_number)))
xlabel('Datapoints');
ylabel('Diameter of J: max|m1-m2|')
title('Estimation Bound of m')
% 
% fname = sprintf('figs/quadrotors/mbound_Task%0.1f_bound%0.1f_step%0.1f.png', flag.traj,loose,step);
% exportgraphics(gcf,fname)
