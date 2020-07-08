function [u] = multigrid2(u,u_,lambda,beta,N,method,k,g,m,I)
%MULTIGRID Summary of this function goes here
%   Detailed explanation goes here
    dx = 1./m;
    z = 4*N*beta*g/dx/dx;
    a = 2.*pi*sqrt(beta.*g);     %damping coeffience
    DxF = (u(:,[2:end end],:) - u)./dx;
    DyF = (u([2:end end],:,:) - u)./dx;
    c1 = sum(sum(u.*I))./sum(sum(u));
    c2 = sum(sum((1 - u).*I))./sum(sum(1 - u));
    Bx = beta.*DxF./(sqrt(1 + beta.^2.*(DxF.^2 + DyF.^2)));
    By = beta.*DyF./(sqrt(1 + beta.^2.*(DxF.^2 + DyF.^2)));
    Du = (Bx - Bx(:,[1 1:end - 1],:))./dx + (By - By([1 1:end - 1],:,:))./dx;
    dE= - g.*Du + lambda.*((I - c1).^2 - (I - c2).^2);
    
    utemp = u;
    [u,dt] = getU(u,u_,I,method,a,z,lambda,dx,beta,dE,g);
    u_ = utemp;

    
    if m > 2
        m = m./2;
        u = u(1:2:end, 1:2:end);
        u_ = u_(1:2:end, 1:2:end);
        I = I(1:2:end, 1:2:end);

        ec = multigrid2(u,u_,lambda,beta,N,method,k,g,m,I);
        e = zeros(m.*2);
        e(1:2:end,1:2:end) = ec;
        e(2:2:end - 1,:) = 0.5.*(e(3:2:end,:) + e(1:2:end - 2,:));
        e(:,2:2:end - 1) = 0.5.*(e(:,3:2:end) + e(:,1:2:end - 2));
        u = e;
    end
    
end
    
