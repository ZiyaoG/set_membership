clear all;
close all;
%% This script analyzes the data from 'uncertaintybound_convergence.m', show a few figures that point out noticable remarks in the estimation of setmembership
load('datas\test1.mat');

%% convergence plot
figure('Name','Convergence')
semilogy(T_series,volume_LSE)
xlabel('T')
ylabel('volume')
hold on
grid on
semilogy(T_series,volume_chisquare)

legend('LSE','chisquare')

%% Figure 1 shows that the estimation set from set membership doesn't always include the true solution
T_t=74;
ineq_LSE_T=set_LSE(:,:,T_t);
ineq_chisquare_T=set_chisquare(:,:,T_t);

figure('Name','incorrect estimate')
h1 = surf(A_mesh,B_mesh,0.8*double(ineq_LSE_T),'FaceColor','blue') ;
set(h1,'LineStyle','none');
xlabel('A');
ylabel('B');
hold on;
h3=plot3(A_star,B_star,1,'b*');
title('LSE Confidence Set p=0.8,T=350')
title(sprintf('LSE Confidence Set p=%d,T=%d', p, T_t));


h2= surf(A_mesh,B_mesh,0.7*double(ineq_chisquare_T),'FaceColor','green') ;
set(h2,'LineStyle','none');

map = [1.0 1.0 1.0
    1.0 0 0
    0 1.0 0
    0 0 1.0];
colormap(map)
% colorbar
alpha(0.5)
zlim([0.01,1.2])
view(0,90)

legend([h1,h2,h3], {'LSE', 'Chi-Squared','Ground Truth'});
xlim([-0.5,3])
ylim([-0.5,3])

%% Figure 2 summarizes all the estimated sets from chisquared and normalizes it. It shows that only about 80% of sets include the true value. 







%% Figure 3 shows that the estimation set could be very small or even doesn't exists, because its existence is probabilistic and only at p (80% in this dataset). I also directly used gaussian distribution to justify this. 
% zoom in figure from convergence plot
%% number of times when estimation fail, i.e., ground truth not in the set
succ_LSE = 0;
succ_Chisquared = 0;
for i=1:T
    succ_LSE= succ_LSE+set_LSE((1+(B_star-B(1))/pixel),(1+(A_star-A(1))/pixel),i);
    succ_Chisquared= succ_Chisquared+set_chisquare((1+(B_star-B(1))/pixel),(1+(A_star-A(1))/pixel),i);
end

prct_LSE=succ_LSE/T
prct_chisquared=succ_Chisquared/T




