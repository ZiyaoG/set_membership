
clear all;
close all;
% load("datas\quadrotors\Trajecory_Task1_tg_f05_M05.mat");
% load("datas\quadrotors\Trajecory_Task1_test.mat");
load("datas\quadrotors\Trajecory_Task1_uniform_f05_M05.mat");

loose=1;
% B=loose*[1;1;1];
B=loose*[0.5;0.5;0.5];
start_point=1;
end_point=2000;
step=100;
Diameter_J=zeros([size(start_point:step:end_point,2),1]);

%% Find two J
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
Constraints=[];
 
Constraints=[Constraints,-J1_x<=0,-J1_y<=0,-J1_z<=0,-J2_x<=0,-J2_y<=0,-J2_z<=0];

for n=start_point:step:end_point
    for j=1:n
        LMI_j1=abs(diag([J1_x,J1_y,J1_z])/dt * (w(:,j+1) - w(:,j)) - M(:,j) + cross(w(:,j), diag([J1_x,J1_y,J1_z]) * w(:,j))) - B <= 0;
        LMI_j2=abs(diag([J2_x,J2_y,J2_z])/dt * (w(:,j+1) - w(:,j)) - M(:,j) + cross(w(:,j), diag([J2_x,J2_y,J2_z]) * w(:,j))) - B <= 0;
        Constraints=[Constraints,LMI_j1,LMI_j2];
    end

    Objective =  abs(J1_x-J2_x)+abs(J1_y-J2_y)+abs(J1_z-J2_z);
    options = sdpsettings('verbose', 0); % Choose a solver that is available and appropriate
    sol = optimize(Constraints, -Objective, options);
    
    if sol.problem == 0
        J1_value = double(diag([J1_x,J1_y,J1_z]));
        J2_value = double(diag([J2_x,J2_y,J2_z]));
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
plot(start_point:step:end_point,Diameter_J)
xlabel('Datapoints');
ylabel('Diameter of J: max||J1-J2||_1')
title('Estimation Bound of J')

% fname = sprintf('figs/quadrotors/Jbound_Task%0.1f_loosebound%0.1f_step%0.1f.png', flag.traj,loose,step);
% exportgraphics(gcf,fname)