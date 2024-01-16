function ineq_gaussian_intersec = calculate_gaussian(X,U,T,A_mesh,B_mesh,sigma_w,p)
%% This function generates the uncertainty bound using chisquare inequality, assuming standard normal distribution
    ineq_gaussian_intersec=zeros(size(A_mesh));
    for i=1:T
        y=X(i+1,:)-A_mesh*X(i,:)-B_mesh*U(i,:);
        t=sigma_w*norminv((1+p)/2);
        ineq_gaussian=abs(y)<t;
        ineq_gaussian_intersec=ineq_gaussian_intersec+(1/(T))*double(ineq_gaussian);
    end
end