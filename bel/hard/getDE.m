% getDE.m
% Calculate the current potential energy and gradient flow
% 
% getDu():Calculate the gradient descent of regularizatiion term
%
% Yuxin Sun
% 1498290038@qq.com
% Georgia Tech
% 2019.9.5
function [e,dE] = getDE(u,I,lambda,dx,beta,g)
   
    c1 = sum(sum(u.*I))./sum(sum(u));
    c2 = sum(sum((1 - u).*I))./sum(sum(1 - u));
    DxF = (u(:,[2:end end],:) - u)./dx;
    DyF = (u([2:end end],:,:) - u)./dx;
    Bx = beta.*DxF./(sqrt(1 + beta.^2.*(DxF.^2 + DyF.^2)));
    By = beta.*DyF./(sqrt(1 + beta.^2.*(DxF.^2 + DyF.^2)));
    Du = (Bx - Bx(:,[1 1:end-1],:))./dx + (By - By([1 1:end-1],:,:))./dx;
%    dE = - g.*getDu(u,dx,beta) + lambda.*((I - c1).^2 - (I - c2).^2);
    dE = - g.*Du + lambda.*((I - c1).^2 - (I - c2).^2);
  
    e = g./beta.*sum(sum(sqrt(1 + beta.^2.*(DxF.^2 + DyF.^2)))) + lambda.*sum(sum(u.*(I - c1).^2)) + lambda.*sum(sum((1 - u).*(I - c2).^2)); 

end