% main_acc_bel.m
% PDE acceleration for TV regularziation of segmentation
% The segmentation scheme we used here is called fuzzy region competition from:
% Mory-Ardon, 2007, "Fuzzy Region Competition: A Convex Two-Phase Segmentation Framework"
%
% getDE(): Calculate the current potential energy and gradient flow
% getU(): Update level set using acceleration scheme
%
% Yuxin Sun
% 1498290038@qq.com
% Georgia Tech
% 2019.9.5
clear all

image = imread('../../cameraman.png');
image = double(image);
I = (image - min(image(:)))./(max(image(:)) - min(image(:)));
load('../../cameraman_noisy.mat') %For noisy image input
[m,n] = size(I);

%Initial level set
u0 = zeros(m,n);
for i = 1:n
    for j = 1:m
        if (sqrt((i - n./2).^2 + (j - m./2)^2.*n./m) - n./4) < 0
            u0(j,i) = 1.0;
        else
            u0(j,i) = 0.0;
        end
    end
end

%parameters
lambda = 0.05;   %Fidelity weight
dx = 1;          %Grid point distance
tmax = 10;
N = 2;           %Dimension
g = 1;           %Regularization weight
method = 1;
%0:gradient descent
%1:1-order
%2:2-order
%3:semi-implicit

tic

for k = [7]                  %times of actual damping coefficient to theortical best damping coefficient
for b = [2]                          %beta in beltrami regularization
for method = [1]
    
    figure
    a = k.*2.*sqrt(b.*g.*(pi.^2)./m./n);     %damping coeffience
    E = zeros(tmax,1);              %Initial potential energy function
    u = u0;                          %u for level set(primal variable)
    u_ = u;                          %u_ for previous level set
    ind = 0;
    t1 = clock;
    
    while (ind < tmax)

        %beta = getBeta(u,dx,b);
        beta = b;
        utemp = u;
        z = 4*N*beta*g/dx/dx;
        [e,dE] = getDE(u,I,lambda,dx,beta,g);
        [u,dt] = getU(u,u_,I,method,a,z,lambda,dx,beta,dE,g);
        u_ = utemp;
        t2 = clock;
        
        if(etime(t2,t1) > ind)
            E(ind + 1) = e;
            ind = ind + 1;
%             surf(u);
%             imshow(uint8(I.*(max(image(:)) - min(image(:))) + min(image(:)))); colormap(gray); hold on
%             contour(u, [0.5 0.5], 'r');hold off
%             title(sprintf('Contour at Level-Set 0 at Time = %d ,lambda=%d,k=%d,b=%d', floor(etime(t2,t1)),lambda,k,b));
            drawnow;
        end
        
    end
    
    %automatically save the energy
    name = strcat('t_hard_b =',num2str(b),',t=',num2str(tmax),',lambda=',num2str(lambda),',k=',num2str(k),'.mat');
    save(name,'E');

end
end
end
toc

    

