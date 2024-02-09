
clear all;
close all;
load("datas\quadrotors\Trajecory_Task1_constMdist05.mat");

loose=1;
B=loose*[0.5;0.5;0.5];
i=2000;
data_number=40;
step=11;
Diameter_J=zeros([data_number,1]);

%% Find two J
J1 = sdpvar(3,3,'diagonal');
J2 = sdpvar(3,3,'diagonal');
Constraints=[];
for l=1:3
    J1_min = -J1(l,l)<=0;
    J2_min = -J2(l,l)<=0;
    Constraints=[Constraints,J1_min,J2_min];
end
for n=1:data_number
    for j=step*(1:n)
        LMI_j1=abs(J1/dt * (w(:,i+j+1) - w(:,i+j)) - M(:,i+j) + cross(w(:,i+j), J1 * w(:,i+j))) - B <= 0;
        LMI_j2=abs(J2/dt * (w(:,i+j+1) - w(:,i+j)) - M(:,i+j) + cross(w(:,i+j), J2 * w(:,i+j))) - B <= 0;
        Constraints=[Constraints,LMI_j1,LMI_j2];
    end
    Distance=J1-J2;
    Distance_vec(1)=Distance(1,1);
    Distance_vec(2)=Distance(2,2);
    Distance_vec(3)=Distance(3,3);
    Objective = -sum(abs(Distance_vec));
    options = sdpsettings('verbose', 0, 'solver', 'mosek'); % Choose a solver that is available and appropriate
    sol = optimize(Constraints, Objective, options);
    
    if sol.problem == 0
        J1_value = double(J1);
        J2_value = double(J2);
        % disp('J1:');
        % disp(J1_value);
        % disp('J2:');
        % disp(J2_value);
        disp('--------------------------------------------')
        fprintf('Diameter at %.1f: max||J1-J2||_1=%f\n',n,sum(abs(J1_value-J2_value),'all'))
        % disp(sum(abs(J1_value-J2_value),'all'));
    else
        disp('No solution found');
        return;
    end
    Diameter_J(n,1)=sum(abs(J1_value-J2_value),'all');
end
%%
% plot(dt*(i+(1:data_number)),Diameter_J)
plot(i+(1:data_number),Diameter_J)
xlabel('Datapoints');
ylabel('Diameter of J: max||J1-J2||_1')
title('Estimation Bound of J')

fname = sprintf('figs/quadrotors/Jbound_Task%0.1f_loosebound%0.1f_step%0.1f.png', flag.traj,loose,step);
exportgraphics(gcf,fname)