clear all;
close all;

%% load data
T_list=[50 100 200 300 500 1000 2000 5000];
diameter_list=zeros([6,size(T_list,2)]);
numfiles=99;
dt=0.002;
Traj=cell(1,numfiles);
for j=1:numfiles
    inputFolder='datas\quadrotors\tg';
    baseFileName = sprintf('Trajecory_%2.2d_Task1_tg_f05_M05.mat', j);
    fullFileName = fullfile(inputFolder, baseFileName);
    Traj{j}=load(fullFileName);
end

Jxx_truth=0.082;
Jyy_truth=0.0845;
Jzz_truth=0.1377;
Jzzyy_truth=Jzz_truth-Jyy_truth;
Jxxzz_truth=Jxx_truth-Jzz_truth;
Jyyxx_truth=Jyy_truth-Jxx_truth;
J_truth=[Jxx_truth 0 0 Jzzyy_truth 0 0;
         0 Jyy_truth 0 0 Jxxzz_truth 0;
         0 0 Jzz_truth 0 0 Jyyxx_truth];
T_index=1;  % for logging T_list
for T=T_list
    %% create X and Y
    for j=1:numfiles
        m_tilde{j}=Traj{j}.M(:,1:T)+Traj{j}.M_expl(:,1:T);
        % tmp=Traj{j}.M_dist(1:T,:);
        % w_tilde{j}=tmp';
        Y{j}=(dt*(m_tilde{j}))';
        tmp=Traj{j}.phi(:,1:T);
        X{j}=tmp';
        X_inv{j}=(X{j}'*X{j})^(-1)*X{j}';
        theta_hat{j}=Y{j}'*X_inv{j}';
    end
    
    %% mean 
    sum_theta=zeros([3,6]);
    for j=1:numfiles
        sum_theta=sum_theta+theta_hat{j};
    end
    mean_theta=sum_theta/numfiles;
    
    %% diameters
    m_hat=zeros([3,6]);
    s_hat=zeros([3,6]);
    tmp=zeros([3,6]);
    %% m
    for j=1:numfiles
        m_hat=m_hat+1/numfiles*(theta_hat{j}-J_truth);
    end
    %% s
    for j=1:numfiles
        tmp=tmp+ (theta_hat{j}-J_truth-m_hat).^2;
    end
    s_hat=(tmp/(numfiles-1)).^0.5;
    %% diameter
    z_star=1.96;
    diameter=2*z_star*s_hat/sqrt(numfiles);
    diameter_flat=[diameter(1,1);diameter(2,2);diameter(3,3);diameter(1,4);diameter(2,5);diameter(3,6)];
    diameter_list(:,T_index)=diameter_flat;
    T_index=T_index+1;
end


%% plot
load("datas\quadrotors\Trajecory_Task1_uniform_f02_M02.mat");
tiledlayout(3,2)
label_list={'Diameter $|J_{xx}|$','Diameter $|J_{yy}|$','Diameter $|J_{zz}|$',
    'Diameter $|J_{zz}-J_{yy}|$','Diameter $|J_{xx}-J_{zz}|$','Diameter $|J_{yy}-J_{xx}|$'};
for i=1:6
    nexttile
    plot(T_list,diameter_list(i,:));
    hold on
    grid on
    xlim([50 5000])
    xlabel('T')
    legend('OLS')
    ylabel(label_list(i),'Interpreter','latex')
end




