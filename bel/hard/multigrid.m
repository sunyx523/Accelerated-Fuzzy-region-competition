function [u_next] = multigrid(lambda,beta,N,method,k,g,m,I,dE)
%MULTIGRID Summary of this function goes here
%   Detailed explanation goes here
    dx = 1./m;
    u0 = zeros(m,m);
    for i = 1:m
       for j = 1:m
           if (sqrt((i - m./2).^2 + (j - m./2)^2) - m./4) < 0
                u0(j,i) = 1.0;
           else
                u0(j,i) = 0.0;
           end
       end
    end
    u = u0;
    u_ = u0;
    z = 4*N*beta*g/dx/dx;
    a = 2.*pi*sqrt(beta.*g);     %damping coeffience
    
    for num = 1:k
        utemp = u;
        [u,dt] = getU(u,u_,I,method,a,z,lambda,dx,beta,dE,g);
        u_ = utemp;
    end

    if m > 2
        m = m./2;
        dx = 1./m;  
        DxF = (u(:,[2:end end],:) - u)./dx;
        DyF = (u([2:end end],:,:) - u)./dx;
        c1 = sum(sum(u.*I))./sum(sum(u));
        c2 = sum(sum((1 - u).*I))./sum(sum(1 - u));
        Bx = beta.*DxF./(sqrt(1 + beta.^2.*(DxF.^2 + DyF.^2)));
        By = beta.*DyF./(sqrt(1 + beta.^2.*(DxF.^2 + DyF.^2)));
        Du = (Bx - Bx(:,[1 1:end - 1],:))./dx + (By - By([1 1:end - 1],:,:))./dx;
        r = dE -(- g.*Du + lambda.*((I - c1).^2 - (I - c2).^2));
        rc = r(1:2:end, 1:2:end);
        I = I(1:2:end, 1:2:end);

        ec = multigrid(lambda,beta,N,method,k,g,m,I,rc);
        e = zeros(m.*2);
        e(1:2:end,1:2:end) = ec;
        e(2:2:end - 1,:) = 0.5.*(e(3:2:end,:) + e(1:2:end - 2,:));
        e(:,2:2:end - 1) = 0.5.*(e(:, 3:2:end) + e(:,1:2:end - 2));
        u = max(min(u + e - u0,1),0);

    end
    u_ = u;
    for num = 1:k
        utemp = u;
        [u,dt] = getU(u,u_,I,method,a,z,lambda,dx,beta,dE,g);
        u_ = utemp;
    end
    u_next = u;    
end
    
