clear all;
close all;
%% Graph: compare point estimation for LSE and gaussian
% load("datas\quadrotors\Trajecory_Task1_uniform_f05_M05.mat");
load("datas\quadrotors\tg\Trajecory20_Task1_tg_f05_M05.mat");

% load("datas\quadrotors\Trajecory_augdyn_Task1_nonoise.mat");
J_vec=[0.082;0.0845;0.1377];

%% LSE
m_est_list=[];
J_est_list=[];
err_m_list=[];
err_J_list=[];
start_point=5;
end_point=200;
step=5;

for T=start_point:step:end_point
    m_est=sdpvar(1,1);
    J_est=sdpvar(3,1);
    Constraints=[];
    Objective_m=0;
    Objective_J=0;
    
    for i=1:T
        Objective_m=Objective_m+( v(:,i+1)-g*e3*dt-v(:,i)+f(1,i)/m_est*R(:,:,i)*e3*dt )'*( v(:,i+1)-g*e3*dt-v(:,i)+f(1,i)/m_est*R(:,:,i)*e3*dt );
        Objective_J=Objective_J+( w(:,i+1)-w(:,i)-diag(1./J_est)*(-cross(w(:,i),diag(J_est)*w(:,i))+M(:,i))*dt )'*( w(:,i+1)-w(:,i)-diag(1./J_est)*(-cross(w(:,i),diag(J_est)*w(:,i))+M(:,i))*dt );
    end
    
    
    options_m = sdpsettings('verbose', 0); 
    options_J = sdpsettings('verbose', 0); 
    sol_m = optimize([], Objective_m, options_m);
    sol_J = optimize([], Objective_J, options_J);
    err_m = abs(m-double(m_est));
    err_J = norm(J_vec-double(J_est),1);
    err_m_list=[err_m_list,err_m];
    err_J_list=[err_J_list,err_J];
    m_est_list=[m_est_list,double(m_est)];
    J_est_list=[J_est_list,double(J_est)];

end

%% Gaussian Process




%% plot
% figure(1)
% plot(start_point:step:end_point,err_m_list)
% xlabel('Datapoints');
% ylabel('$||m-m^*||$','interpreter','latex')
% title('SM Estimation of m')
% legend('LSE','GP')
% ylim([0 0.5])
% J
% 
% figure(2)
% plot(start_point:step:end_point,err_J_list)
% xlabel('Datapoints');
% ylabel('$||J-J^*||_2^2$','interpreter','latex')
% title('SM Estimation of J')
% legend('LSE','GP')
% ylim([0 0.01])

save('datas\quadrotors\tg\LSE\Trajecory20_Task1.mat','m_est_list','J_est_list','err_J_list','err_m_list');
% Theta_hat=((Z'*Z+lambda*eye(size(Z,2)))^(-1))*Z'*X(2:T+1,:);