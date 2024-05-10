
clear all;
close all;
% load("datas\quadrotors\Trajecory_Task1_tg_f05_M05.mat");
% load("datas\quadrotors\Trajecory_Task1_test.mat");
load("datas\quadrotors\Trajecory_Task1_uniform_f05_M05.mat");

loose=1;
% B=loose*[1;1;1];
B=loose*[0.5;0.5;0.5];
start_point=100;
end_point=2000;
step=100;
Diameter_J=zeros([size(start_point:step:end_point,2),1]);

%% Find two J
J1_1 = sdpvar(1,1);
J1_2 = sdpvar(1,1);
J1_3 = sdpvar(1,1);
J2_1 = sdpvar(1,1);
J2_2 = sdpvar(1,1);
J2_3 = sdpvar(1,1);
Constraints=[];
 
Constraints=[Constraints,-J1_1<=0,-J1_2<=0,-J1_3<=0,-J2_1<=0,-J2_2<=0,-J2_3<=0];

for n=start_point:step:end_point
    for j=1:n
        LMI_j1=abs(diag([J1_1,J1_2,J1_3])/dt * (w(:,j+1) - w(:,j)) - M(:,j) + cross(w(:,j), diag([J1_1,J1_2,J1_3]) * w(:,j))) - B <= 0;
        LMI_j2=abs(diag([J2_1,J2_2,J2_3])/dt * (w(:,j+1) - w(:,j)) - M(:,j) + cross(w(:,j), diag([J2_1,J2_2,J2_3]) * w(:,j))) - B <= 0;
        Constraints=[Constraints,LMI_j1,LMI_j2];
    end

    Objective =  abs(J1_1-J2_1)+abs(J1_2-J2_2)+abs(J1_3-J2_3);
    options = sdpsettings('verbose', 0); % Choose a solver that is available and appropriate
    sol = optimize(Constraints, -Objective, options);
    
    if sol.problem == 0
        J1_value = double(diag([J1_1,J1_2,J1_3]));
        J2_value = double(diag([J2_1,J2_2,J2_3]));
        % disp('J1:');
        % disp(J1_value);
        % disp('J2:');
        % disp(J2_value);
        disp('--------------------------------------------')
        fprintf('Diameter at %.1f: max||J1-J2||_2=%f\n',n,norm(J1_value-J2_value))
        % disp(sum(abs(J1_value-J2_value),'all'));
    else
        disp('No solution found');
        return;
    end
    Diameter_J(n/step,1)=norm(J1_value-J2_value);
end
%%
% plot(dt*(i+(1:data_number)),Diameter_J)
plot(start_point:step:end_point,Diameter_J)
xlabel('Datapoints');
ylabel('Diameter of J: max||J1-J2||_1')
title('Estimation Bound of J')

% fname = sprintf('figs/quadrotors/Jbound_Task%0.1f_loosebound%0.1f_step%0.1f.png', flag.traj,loose,step);
% exportgraphics(gcf,fname)