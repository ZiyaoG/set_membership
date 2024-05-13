clear all;
close all;
load("datas\quadrotors\tg\Trajecory_Task1_tg_f05_M05.mat");

%%
i=1000;

A=[w(1,i+1)-w(1,i) 0 0 dt*w(2,i)*w(3,i) 0 0;
    0 w(2,i+1)-w(2,i) 0 0 dt*w(3,i)*w(1,i) 0;
    0 0 w(3,i+1)-w(3,i) 0 0 dt*w(1,i)*w(2,i)];
b=M(:,i)+M_expl(:,i)+0.5*[1;1;1];

V=lcon2vert(A,b)