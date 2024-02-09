function [x_next,v_next,R_next,w_next]=quadrotor_dyn(xt,vt,Rt,wt,ft,Mt,f_dist_t,M_dist_t,dt,m,g,J,e3)
    x_next=xt+vt*dt;
    v_next=vt+(g*e3-ft/m*Rt*e3)*dt+f_dist_t/m*dt;
    R_next=Rt+(Rt*wedge(wt))*dt;
    w_next=wt+(-J^(-1)*cross(wt,J*wt)+J^(-1)*Mt+J^(-1)*M_dist_t)*dt;
end