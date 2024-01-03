%% Copyright Ziyao Guo UIUC 
clear all;
close all;
%% Simulation parameters
T=10; % 5 250
p=0.8; % confidence level
delta=1-p;

% system parameters
A_star=0.8;
B_star=1;
mu_w=0;
mu_u=0;
sigma_w=0.5;
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

figure('Name','LSE')
h = pcolor(A_mesh,B_mesh,double(ineq_LSE)) ;
h.EdgeColor = 'none' ;
xlabel('A');
ylabel('B');
hold on;
plot(A_star,B_star,'*');
title('LSE Confidence Set 80%')

%% plot hoefding induced uncertainty set, using a recursive way

ineq_hoefding_recur = calculate_hoeffding(X,U,T,A_mesh,B_mesh,sigma_w,p);

figure('Name','Hoeffding')
h = pcolor(A_mesh,B_mesh,ineq_hoefding_recur) ;
h.EdgeColor = 'none' ;
xlabel('A');
ylabel('B');
hold on;
plot(A_star,B_star,'*');
title('Hoeffding Induced Set 80%')

%% Chebyshev's inequality, with Yi=x_i+1 - Ax_i-Bu_i
ineq_chebyshev = calculate_chebyshev(X,U,T,A_mesh,B_mesh,sigma_w,p);

figure('Name','Chebyshev')
h = pcolor(A_mesh,B_mesh,double(ineq_chebyshev)) ;
h.EdgeColor = 'none' ;
xlabel('A');
ylabel('B');
hold on;
plot(A_star,B_star,'*');
title('Chebyshev Induced Set 80%')

%% Convergence test
T_series=5:3:350;
length=size(T_series,2);
volume_LSE=zeros(length,1);
volume_chebyshev=zeros(length,1);
volume_hoeffding=zeros(length,1);

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
    ineq_chebyshev = calculate_chebyshev(X,U,T,A_mesh,B_mesh,sigma_w,p);
    tic;
    ineq_hoefding_recur = calculate_hoeffding(X,U,T,A_mesh,B_mesh,sigma_w,p);
    toc;
    % log data
    volume_LSE(t)=sum(ineq_LSE,'all')*pixel^2;
    volume_chebyshev(t)=sum(ineq_chebyshev,'all')*pixel^2;
    volume_hoeffding(t)=sum(ineq_hoefding_recur,'all')*pixel^2;
    fprintf('Finished analysis for T=%d\n', T);
    
end

figure('Name','Convergence')
semilogy(T_series,volume_LSE)
xlabel('T')
ylabel('volume')
hold on
grid on
semilogy(T_series,volume_chebyshev)
semilogy(T_series,volume_hoeffding)
legend('LSE','Chebyshev','Hoeffding')





