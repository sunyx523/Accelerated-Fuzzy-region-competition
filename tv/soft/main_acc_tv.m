% main_acc_tv.m
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

image = imread('../cameraman.png');
image = double(image);
I = (image - min(image(:)))./(max(image(:)) - min(image(:)));
load('../cameraman_noisy.mat') %For noisy image input
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

lambda = 1;       %Fidelity weight
g = 1;            %Regularization weight
N = 2;            %Dimension
Q = 1./255;       %Quantization level
dx = 1;           %Grid point distance
num = 1;          %Iteration index
itmax = 10000;    %Maximum iteration
h = 1;
method = 1;
dumax = 1;
%0:gradient descent
%1:1-order
%2:2-order
%3:semi-implicit

u = u0;
u_ = u;

tic
for h = [1]
for Q = [1./255]
for method = [1]  
for a = [0.5]    %damping coeffience
    
    figure
    E = zeros(itmax,1);            %Initla potential energy function
    u = u0;                        %u for level set(primal variable)
    u_ = u;                        %u_ for previous level set
    du  =100;
    
    while  (du > dumax && num < itmax)    
  %  while num < itmax
        
        utemp = u;
        z = 4.*sqrt(N).*g./Q./dx;
        [e,dE] = getDE(u,I,lambda,dx,g,h);
        [u] = getU(u,u_,I,method,a,z,lambda,dx,dE,g);
        E(num) = e;
        u_ = utemp;

%        visualize contour
        if(mod(num,20) == 0)
            imshow(uint8(I.*(max(image(:)) - min(image(:))) + min(image(:)))); colormap(gray); hold on
            contour(u, [0.5 0.5], 'r');hold off
            title(sprintf('Clean = %d ,a=%d,m=%d,Q=%d',num,a,method,Q));
            drawnow;
        end

        num = num + 1;
    end
    
    %automatically save the energy
    num = 1;
    name = strcat('qua_Q =',num2str(1./Q),',i=',num2str(itmax),',a=',num2str(a),'h=',num2str(h),',method=',num2str(method),'.mat');
    save(name,'E');
    du = max(abs(u(:) - u_(:)));
    
    if(method == 0)
        break;
    end

end
end
end
end
toc
