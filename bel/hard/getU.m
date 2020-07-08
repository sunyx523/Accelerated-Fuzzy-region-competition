% getU.m
% Update level set using acceleration scheme
%
% Yuxin Sun
% 1498290038@qq.com
% Georgia Tech
% 2019.9.5
function [u,dt] = getU(u,u_,I,method,a,z,lambda,dx,beta,dE,g)
    switch method
        case 0
            dt = 2/z;
            u = max(min(u - dt*dE,1),0);
        case 1
            dt = sqrt(4/z + (a/z)^2) + a/z;
            u = max(min(u + (u - u_)/(1 + a*dt) - dt^2*dE/(1 + a*dt),1),0);
            u(1,:) = 0;
            u(:,1) = 0;
            u(size(u,1),:) = 0;
            u(:,size(u,2)) = 0;
        case 2
            dt = 2/sqrt(z);
            u = max(min(u + (2 - a*dt)*(u - u_)/(2 + a*dt) - 2*dt^2*dE/(2 + a*dt),1),0);
        case 3
            dt = 2/sqrt(3*z);
            v = max(min(u + (2 - a*dt)*(u - u_)/(2 + a*dt),1),0);
            [~,dE] = getDE(u,I,lambda,dx,beta,g);
            u = max(min(v - 2*dt^2*dE/(2 + a*dt),1),0);
     end
end