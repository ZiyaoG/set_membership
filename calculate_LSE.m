function ineq_LSE = calculate_LSE(X,Z,T,A_mesh,B_mesh,lambda,delta,L,S)
%% This function generates the uncertainty bound using chebyshev inequality
    Theta_hat=((Z'*Z+lambda*eye(size(Z,2)))^(-1))*Z'*X(2:T+1,:);

    V_t=lambda*eye(size(Z,2));
    for i=1:T
        V_t=V_t+Z(i,:)'*Z(i,:);
    end
    beta_t=(size(X,2)*L*(2*log(det(V_t)^0.5*det(lambda*eye(size(Z,2)))^(-0.5)/delta))^0.5+lambda^0.5*S)^2;
    
    
    ineq_LSE=V_t(1,1)*(A_mesh-Theta_hat(1)).^2+(V_t(1,2)+V_t(2,1))*(A_mesh-Theta_hat(1)).*(B_mesh-Theta_hat(2))+V_t(2,2)*(B_mesh-Theta_hat(2)).^2<beta_t;

end