% getDu.m
% Calculate the gradient descent of Beltramiregularizatiion term
%
% Yuxin Sun
% 1498290038@qq.com
% Georgia Tech
% 2019.9.5
function Du = getDu(u,dx,beta)

    DxF = (u(:,[2:end end],:) - u)./dx;
    DyF = (u([2:end end],:,:) - u)./dx;
    Bx = beta.*DxF./(sqrt(1 + beta.^2.*(DxF.^2 + DyF.^2)));
    By = beta.*DyF./(sqrt(1 + beta.^2.*(DxF.^2 + DyF.^2)));
    Du = (Bx - Bx(:,[1 1:end-1],:))./dx + (By - By([1 1:end-1],:,:))./dx;

end