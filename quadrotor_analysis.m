% Analysis the datas from geometric controlled quadrototrs and use 
% set membership to give estimation of the m and J. 
clear all;
close all;
load("datas\quadrotors\sincos05.mat");
m_est.value=4:0.002:5;
msize=size(m_est.value);
m_est.flag=ones(msize);
Bound=3;
for i=1:msize(2)
    for j=1:(T/dt-1)
        m_tmp=abs(m_est.value(i)/dt*(v(:,j+1)-v(:,j))-m_est.value(i)*g*e3+f(j)*R(:,:,j)*e3);
        J_tmp=
        if sum(double(m_tmp>Bound))>0
            m_est.flag(i)=0;
            break;
        end
    end
end

%% plot
plot(m_est.value(1,:),m_est.flag(1,:),'o')
hold on
plot(m,1,'*')
legend('estimation','Ground Truth')

fname = sprintf('figs/quadrotors/Task%0.3f_Bound%0.3f.png', flag,Bound);
exportgraphics(gcf,fname)

