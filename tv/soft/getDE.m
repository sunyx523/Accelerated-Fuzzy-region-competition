% getDE.m
% Calculate the current potential energy and gradient flow
% 
% getDu():Calculate the gradient descent of regularizatiion term
%
% Yuxin Sun
% 1498290038@qq.com
% Georgia Tech
% 2019.9.5

function [e,dE] = getDE(u,I,lambda,dx,g,h)

%     r1 = u;
%     r1(u < 0) = 1;
%     r1(u >= 0) = 0;
%     r2 = u;
%     r2(u > 1) = 1;
%     r2(u <= 1) = 0; 
%     h1 = r1.*(-u);
%     h2 = r2.*(u - 1);
    
%     r1 = -u;
%     r1(r1 < 0) = 0;
%     r2 = u - 1;
%     r2(r2 < 0) = 0;
%     h1 = 0.5.*r1.^2;
%     h2 = 0.5.*r2.^2;

    c1 = sum(sum(u.*I))./sum(sum(u));
    c2 = sum(sum((1 - u).*I))./sum(sum(1 - u));
    dE = - g.*getDu(u,dx) + lambda.*((I - c1).^2 - (I - c2).^2);
%    dE = - g.*getDu(u,dx) + lambda.*((I - c1).^2 - (I - c2).^2) + h.*(r2 - r1);
   
    DxF = (u(:,[2:end end],:) - u)./dx;
    DyF = (u([2:end end],:,:) - u)./dx; 
    grad = sqrt(DxF.^2 + DyF.^2);
    e = g.*sum(sum(grad)) + lambda.*sum(sum(u.*(I - c1).^2)) + lambda.*sum(sum((1 - u).*(I - c2).^2));
%    e = lambda.*sum(sum(u.*(I - c1).^2)) + lambda.*sum(sum((1 - u).*(I - c2).^2));
%    e = g.*sum(sum(grad)) + lambda.*sum(sum(u.*(I - c1).^2)) + lambda.*sum(sum((1 - u).*(I - c2).^2)) + h.*sum(sum(h1 + h2));
    
end