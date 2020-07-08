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

image = imread('../../brain.jpg');
image = rgb2gray(image);
image = double(image);
I = (image - min(image(:)))./(max(image(:)) - min(image(:)));
load('../../cameraman_noisy.mat') %For noisy image input
I = imresize(I,[512 512]);
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

lambda = 500;       %Fidelity weight
g = 1;            %Regularization weight
N = 2;            %Dimension
Q = 1./255;       %Quantization level
dx = 1/m;           %Grid point distance
tmax = 1000;      %Maximum time
method = 1;
dumax = 1;
%0:gradient descent
%1:1-order
%2:2-order
%3:semi-implicit

u = u0;
u_ = u;

tic
for Q = [1./255]
for method = [1]  
for a = [1]    %damping coeffience
    
    figure
    E = zeros(tmax,1);            %Initla potential energy function
    DU = zeros(tmax,1);
    u = u0;                        %u for level set(primal variable)
    u_ = u;                        %u_ for previous level set
    ind = 0;
    du  =0.1;
    t1 = clock;
    
%    while  (du > dumax && num < itmax)    
    while ind < tmax/10
        
        utemp = u;
        z = 4.*sqrt(N).*g./Q./dx;
        [e,dE] = getDE(u,I,lambda,dx,g);
        [u] = getU(u,u_,I,method,a,z,lambda,dx,dE,g);
        u_ = utemp;
        t2 = clock;

        if(etime(t2,t1) > ind)
            E(int32(ind*10 + 1)) = e;
            DU(int32(ind*10 + 1)) = du;
            ind = ind + 0.1;
        end
%        visualize contour
%         if(mod(int32(ind*10),10) == 0)
%             imshow(uint8(I.*(max(image(:)) - min(image(:))) + min(image(:)))); colormap(gray); hold on
%             contour(u, [0.5 0.5], 'r');hold off
%             title(sprintf('Contour at Level-Set 0 at Time = %d ,a=%d, lambda=%d', ind, a, lambda));
%             drawnow;
%         end
        du = max(abs(u(:) - u_(:)));
    end
    
    %automatically save the energy
    ind = 1;
    name = strcat('qua_Q =',num2str(1./Q),',i=',num2str(tmax),',a=',num2str(a),',method=',num2str(method),'.mat');
    save(name,'E');

    
end
end
end
toc
