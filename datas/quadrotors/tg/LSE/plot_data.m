clear all;
close all;

%%
m=4.34;
J_vec=[0.082;0.0845;0.1377];
J_vec_list=repmat(J_vec,1,40);
start_point=5;
end_point=200;
step=5;

Traj1=load('Trajecory1_Task1.mat');
Traj2=load('Trajecory2_Task1.mat');
Traj3=load('Trajecory3_Task1.mat');
Traj4=load('Trajecory4_Task1.mat');
Traj5=load('Trajecory5_Task1.mat');
Traj6=load('Trajecory6_Task1.mat');
Traj7=load('Trajecory7_Task1.mat');
Traj8=load('Trajecory8_Task1.mat');
Traj9=load('Trajecory9_Task1.mat');
Traj10=load('Trajecory10_Task1.mat');
Traj11=load('Trajecory11_Task1.mat');
Traj12=load('Trajecory12_Task1.mat');
Traj13=load('Trajecory13_Task1.mat');
Traj14=load('Trajecory14_Task1.mat');
Traj15=load('Trajecory15_Task1.mat');
Traj16=load('Trajecory16_Task1.mat');
Traj17=load('Trajecory17_Task1.mat');
Traj18=load('Trajecory18_Task1.mat');
Traj19=load('Trajecory19_Task1.mat');
Traj20=load('Trajecory20_Task1.mat');

m_est_mean=0.05*(Traj1.m_est_list+Traj2.m_est_list+Traj3.m_est_list+Traj4.m_est_list+ ...
    Traj5.m_est_list+Traj6.m_est_list+Traj7.m_est_list+Traj8.m_est_list+Traj9.m_est_list+ ...
    Traj10.m_est_list+Traj11.m_est_list+Traj12.m_est_list+Traj13.m_est_list+Traj14.m_est_list+ ...
    Traj15.m_est_list+Traj16.m_est_list+Traj17.m_est_list+Traj18.m_est_list+Traj19.m_est_list+ ...
    Traj20.m_est_list);
J_est_mean=0.05*(Traj1.J_est_list+Traj2.J_est_list+Traj3.J_est_list+Traj4.J_est_list+ ...
    Traj5.J_est_list+Traj6.J_est_list+Traj7.J_est_list+Traj8.J_est_list+Traj9.J_est_list+ ...
    Traj10.J_est_list+Traj11.J_est_list+Traj12.J_est_list+Traj13.J_est_list+Traj14.J_est_list+ ...
    Traj15.J_est_list+Traj16.J_est_list+Traj17.J_est_list+Traj18.J_est_list+Traj19.J_est_list+ ...
    Traj20.J_est_list);

err_m = abs(m-m_est_mean);
tmp = abs(J_vec_list-J_est_mean);
err_J = tmp(1,:)+tmp(2,:)+tmp(3,:);

%%
figure(1)
plot(start_point:step:end_point,err_m)
xlabel('Datapoints');
ylabel('$||m-m^*||$','interpreter','latex')
title('LSE Estimation of m')
legend('LSE','GP')
ylim([0 0.05])

figure(2)
plot(start_point:step:end_point,err_J)
xlabel('Datapoints');
ylabel('$||J-J^*||_2^2$','interpreter','latex')
title('LSE Estimation of J')
legend('LSE','GP')
ylim([0 0.001])