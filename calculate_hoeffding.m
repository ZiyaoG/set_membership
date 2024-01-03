function ineq_hoefding_recur = calculate_hoeffding(X,U,T,A_mesh,B_mesh,sigma_w,p)
%% This function generates the uncertainty bound using Hoeffding inequality
    T_recur=5:T;
    ineq_hoefding_recur=ones(size(A_mesh));
    for i=1:size(T_recur,2)
        T_curr=T_recur(i);
        X_curr=X(1:T_curr+1,:);
        U_curr=U(1:T_curr,:);
        % W_curr=W(1:T_curr,:);
    
        t=(-2*T_curr*sigma_w^2*log((1-p)/2))^0.5;
        Y_sum=zeros(size(A_mesh));
        for j=1:T_curr
            Y_j=X_curr(j+1,:)*ones(size(A_mesh))-A_mesh*X(j,:)-B_mesh*U_curr(j,:);
            Y_sum=Y_sum+Y_j;
            % w_j=W_curr(j,:);
        end
        ineq_hoefding_curr=double(abs(Y_sum)<t);
        ineq_hoefding_recur=ineq_hoefding_recur.*ineq_hoefding_curr;
    end

end