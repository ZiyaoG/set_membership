clear all;
close all;
%%
% load("datas\quadrotors\tg\Trajecory_01_Task1_tg_f05_M05.mat");
numfiles=2;
dt=0.002;
Traj=cell(1,numfiles);
for j=1:numfiles
    inputFolder='datas\quadrotors\tg';
    baseFileName = sprintf('Trajecory_%2.2d_Task1_tg_f05_M05.mat', j);
    fullFileName = fullfile(inputFolder, baseFileName);
    Traj{j}=load(fullFileName);
end

loose=1;
B=loose*[0.5;0.5;0.5]*dt;
% Diameter_J=zeros([size(start_point:step:end_point,2),1]);
T_list=[50 100 200 300 500 1000 2000];
diameter_list=zeros([6,size(T_list,2),numfiles]);

%% Find two J
for j=1:numfiles
    fprintf('------------------------------------------------------\n')
    fprintf('Trajectory %.0f: \n',j)
    J1_x = sdpvar(1,1);
    J1_y = sdpvar(1,1);
    J1_z = sdpvar(1,1);
    J1_zy = sdpvar(1,1);
    J1_xz = sdpvar(1,1);
    J1_yx = sdpvar(1,1);
    J2_x = sdpvar(1,1);
    J2_y = sdpvar(1,1);
    J2_z = sdpvar(1,1);
    J2_zy = sdpvar(1,1);
    J2_xz = sdpvar(1,1);
    J2_yx = sdpvar(1,1);
    theta1_star=[J1_x 0 0 J1_zy 0 0;
                 0 J1_y 0 0 J1_xz 0;
                 0 0 J1_z 0 0 J1_yx];
    theta2_star=[J2_x 0 0 J2_zy 0 0;
                 0 J2_y 0 0 J2_xz 0;
                 0 0 J2_z 0 0 J2_yx];
    
    Constraints=[];
    Constraints=[Constraints,-J1_x<=0,-J1_y<=0,-J1_z<=0,-J2_x<=0,-J2_y<=0,-J2_z<=0];
    T_index=1; % for logging
    
    for T=T_list
        for i=1:T
            LMI_j1=abs(theta1_star*Traj{j}.phi(:,i)-(Traj{j}.M(:,i)+Traj{j}.M_expl(:,i))*dt)-B<=0;
            LMI_j2=abs(theta2_star*Traj{j}.phi(:,i)-(Traj{j}.M(:,i)+Traj{j}.M_expl(:,i))*dt)-B<=0;
            Constraints=[Constraints,LMI_j1,LMI_j2];
        end
    
        Objective =  abs(J1_x-J2_x)+abs(J1_y-J2_y)+abs(J1_z-J2_z)+abs(J1_zy-J2_zy)+abs(J1_xz-J2_xz)+abs(J1_yx-J2_yx);
        options = sdpsettings('verbose', 0); % Choose a solver that is available and appropriate
        sol = optimize(Constraints, -Objective, options);
        
        if sol.problem == 0
            J1_flat=double([J1_x, J1_y, J1_z, J1_zy, J1_xz, J1_yx]);
            J2_flat=double([J2_x, J2_y, J2_z, J2_zy, J2_xz, J2_yx]);
            diameter_list(:,T_index,j)=abs(J1_flat-J2_flat);
            
            fprintf('Diameter at %.1f: [', T);
            fprintf('%g ', abs(J1_flat-J2_flat));
            fprintf(']\n');
        else
            disp('No solution found');
            return;
        end
        T_index=T_index+1;
    end
end

diameter_list_sum=zeros([6,size(T_list,2)]);
for i=1:numfiles
    diameter_list_sum=diameter_list_sum+diameter_list(:,:,i);
end
diameter_list_mean=diameter_list_sum/numfiles;
%% plot
tiledlayout(3,2)
label_list={'Diameter $|J_{xx}|$','Diameter $|J_{yy}|$','Diameter $|J_{zz}|$',
    'Diameter $|J_{zz}-J_{yy}|$','Diameter $|J_{xx}-J_{zz}|$','Diameter $|J_{yy}-J_{xx}|$'};
for i=1:6
    nexttile
    plot(T_list,diameter_list_mean(i,:));
    grid on
    % xlim(T_list)
    xlabel('T')
    legend('SM')
    ylabel(label_list(i),'Interpreter','latex')
    ytickformat('%,.2f')
end

save('datas\quadrotors\tg\SM12.mat','diameter_list','diameter_list_mean');