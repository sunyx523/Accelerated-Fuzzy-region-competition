% getDE.m
% Calculate the current potential energy and gradient flow
% 
% getDu():Calculate the gradient descent of regularizatiion term
%
% Yuxin Sun
% 1498290038@qq.com
% Georgia Tech
% 2019.9.5

function [e,dE] = getDE(u,I,lambda,dx,g)

    c1 = sum(sum(u.*I))./sum(sum(u));
    c2 = sum(sum((1 - u).*I))./sum(sum(1 - u));
    dE = - g.*getDu(u,dx) + lambda.*((I - c1).^2 - (I - c2).^2);
   
    DxF = (u(:,[2:end end],:) - u)./dx;
    DyF = (u([2:end end],:,:) - u)./dx; 
    grad = sqrt(DxF.^2 + DyF.^2);
    e = g.*sum(sum(grad)) + lambda.*sum(sum(u.*(I - c1).^2)) + lambda.*sum(sum((1 - u).*(I - c2).^2));
   
end