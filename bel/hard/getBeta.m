%getBeta.m
%compute the maximum beta for each iteration
%
% Yuxin Sun
% 1498290038@qq.com
% Georgia Tech
% 2019.9.5
function beta = getBeta(u,dx,b)

    Dx = (u(:,[2:end end],:) - u)./dx;
    Dy = (u([2:end end],:,:) - u)./dx;
    grad = Dx.^2 + Dy.^2;
    gradm = min(grad(:));
    beta = b./(sqrt(1 + b.^2.*gradm));
    
end