function a = getDamp(u,dx,b,g,m,n)
%GETDAMP Summary of this function goes here
%   Detailed explanation goes here
    Dx = (u(:,[2:end end],:) - u)./dx;
    Dy = (u([2:end end],:,:) - u)./dx;
    grad = Dx.^2 + Dy.^2;
    dom = sqrt(1 + b.^2.*grad);
    domAvg = m.*n./sum(sum(dom(:)));
    a = 2.*pi.*sqrt(b.*g/m/n).*sqrt(domAvg);
    

end

