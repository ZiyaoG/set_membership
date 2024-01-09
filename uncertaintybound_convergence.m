%% Ziyao Guo UIUC 
clear all;
close all;
%% Simulation parameters
T=350; % 5 250
p=0.8; % confidence level
delta=1-p;

% system parameters
A_star=0.8;
B_star=1;
mu_w=0;
mu_u=0;
% sigma_w=0.5;
sigma_w=1;
sigma_u=0.5;
w_max=1;
u_max=1;
x_0=0;

lambda=0.1;
L=sigma_w;
S=sqrt(A_star^2+B_star^2);

% for truncated Gaussian
gp=makedist('Normal','mu',mu_u,'sigma',sigma_u);
gp_trun=truncate(gp,-u_max,u_max);

pixel=0.01;
A=-2:pixel:4;
B=-5:pixel:7;
[A_mesh,B_mesh]=meshgrid(A,B);

X=zeros(T+1,size(A_star,1));
U=zeros(T,size(B_star,1));
W=zeros(T,size(B_star,1));
Z=zeros(T,(size(A_star,1)+size(B_star,1)));
% initial value 
X(1,:)=0;

% data collection
for i=1:T
    x_t=X(i,:);
    u_t=random(gp_trun);
    w_t=normrnd(mu_w,sigma_w);
    x_next = A_star*x_t+B_star*u_t+w_t;

    % logging data
    W(i,:)=w_t;
    U(i,:)=u_t;
    X(i+1,:)=x_next;
    Z(i,:)=[x_t',u_t'];
end
%% LSE confidence bound from Abbasi-Yadkori and Szepesvari [2011] 
ineq_LSE = calculate_LSE(X,Z,T,A_mesh,B_mesh,lambda,delta,L,S);

%% hoefding induced uncertainty set, using a recursive way
% ineq_hoefding_recur = calculate_hoeffding(X,U,T,A_mesh,B_mesh,sigma_w,p);

%% Chebyshev's inequality, with Yi=x_i+1 - Ax_i-Bu_i
ineq_chebyshev = calculate_chebyshev(X,U,T,A_mesh,B_mesh,sigma_w,p);

%% chisquare
ineq_chisquare = calculate_chisquare(X,U,T,A_mesh,B_mesh,p);

%% gaussian - theretical best of set membership
ineq_gaussian = calculate_gaussian(X,U,T,A_mesh,B_mesh,p);

%% plot
%---could be improved---------------------------------------------%
h1 = surf(A_mesh,B_mesh,double(ineq_LSE),'FaceColor','blue') ;
set(h1,'LineStyle','none');
xlabel('A');
ylabel('B');
hold on;
plot3(A_star,B_star,1,'b*');
title('LSE Confidence Set p=0.8,T=350')

% h2 = surf(A_mesh,B_mesh,0.5*ineq_hoefding_recur) ;
% set(h2,'LineStyle','none');

% h3 = surf(A_mesh,B_mesh,0.8*double(ineq_chebyshev)) ;
% set(h3,'LineStyle','none');

h4= surf(A_mesh,B_mesh,double(ineq_chisquare),'FaceColor','green') ;
set(h4,'LineStyle','none');

h5= surf(A_mesh,B_mesh,ineq_gaussian,'FaceColor','red') ;
set(h5,'LineStyle','none');

map = [1.0 1.0 1.0
    1.0 0 0
    0 1.0 0
    0 0 1.0];
colormap(map)
% colorbar

alpha(0.5)
legend([h1,h4,h5], {'LSE', 'Chi-Squared','Gaussian'});
zlim([0.01,1.2])
% view(0,90)
%% Convergence test
T_series=5:1:350;
length=size(T_series,2);
volume_LSE=zeros(length,1);
volume_chebyshev=zeros(length,1);
volume_hoeffding=zeros(length,1);
volume_chisquare=zeros(length,1);
% log all the estimated sets
set_chisquare=zeros(size(A_mesh));
set_LSE=zeros(size(A_mesh));

for t=1:length
    T=T_series(t); 

    % Initialization
    X=zeros(T+1,size(A_star,1));
    U=zeros(T,size(B_star,1));
    W=zeros(T,size(B_star,1));
    Z=zeros(T,(size(A_star,1)+size(B_star,1)));

    % Collect data
    for i=1:T
        x_t=X(i,:);
        u_t=random(gp_trun);
        w_t=normrnd(mu_w,sigma_w);
        x_next = A_star*x_t+B_star*u_t+w_t;

        % logging data
        W(i,:)=w_t;
        U(i,:)=u_t;
        X(i+1,:)=x_next;
        Z(i,:)=[x_t',u_t'];
    end

    ineq_LSE = calculate_LSE(X,Z,T,A_mesh,B_mesh,lambda,delta,L,S);
    ineq_chisquare = calculate_chisquare(X,U,T,A_mesh,B_mesh,p);
    ineq_chebyshev = calculate_chebyshev(X,U,T,A_mesh,B_mesh,sigma_w,p);
    % tic;
    % ineq_hoefding_recur = calculate_hoeffding(X,U,T,A_mesh,B_mesh,sigma_w,p);
    % toc;
    % log data
    volume_LSE(t)=sum(ineq_LSE,'all')*pixel^2;
    volume_chebyshev(t)=sum(ineq_chebyshev,'all')*pixel^2;
    set_LSE(:,:,i)=ineq_LSE;
    % volume_hoeffding(t)=sum(ineq_hoefding_recur,'all')*pixel^2;
    volume_chisquare(t)=sum(ineq_chisquare,'all')*pixel^2;
    set_chisquare(:,:,i)=ineq_chisquare;
    fprintf('Finished analysis for T=%d\n', T);
end

figure('Name','Convergence')
semilogy(T_series,volume_LSE)
xlabel('T')
ylabel('volume')
hold on
grid on
semilogy(T_series,volume_chebyshev)
% semilogy(T_series,volume_hoeffding)
semilogy(T_series,volume_chisquare)

% legend('LSE','Chebyshev','Hoeffding','chisquare')
legend('LSE','Chebyshev','chisquare')


% save('datas/test1')


