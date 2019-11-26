% main_pd_tv.m
% Chambolle-Pock primal-dual for TV regularization
% The segmentation scheme we used here is called fuzzy region competition from:
% Mory-Ardon, 2007, "Fuzzy Region Competition: A Convex Two-Phase Segmentation Framework"
% The convex energy is minimized for instance by using Chambolle-Pock primal-dual algorithm
% Chambolle-Pock, 2011, "A First-Order Primal-Dual Algorithm for Convex Problems with Applications to Imaging"
%
% Yuxin Sun
% 1498290038@qq.com
% Georgia Tech
% 2019.9.5
clear all

image = imread('../cameraman.png');
image = double(image);
%I = (image - min(image(:)))./(max(image(:)) - min(image(:)));
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

px = zeros(m,n); %Initial projector(dual variable) on x direction
py = zeros(m,n); %Initial projector(dual variable) on y dirextion

lambda = 1;      %Fidelity weight
dx = 1;          %Grid point distance
itmax = 10000;   %Maximum iterations
num = 1;         %Iteration index
dumax = 0.01;
L2 = 8;

tic
for n = [4]
for tau = [1]             %tau parameter from the Paper
    
    figure
    sigma = n./L2./tau;   %sigma parameter from the paper
    sigma = sqrt(1./L2);
    tau = sqrt(1./L2);
    E = zeros(itmax,1);   %Initla potential energy function

    u = u0;               %u for level set(primal variable)
    u_ = u;               %u_ for previous level set
    num = 1;
    du = 100;
    while  (du > dumax && num < itmax)
%    while num < itmax
        
        utemp = u;
        c1 = sum(sum(u.*I))/sum(sum(u));
        c2 = sum(sum((1 - u).*I))/sum(sum(1 - u));
        r = lambda.*((I - c1).^2 - (I - c2).^2);

        Dux = (u(:,[2:end end],:) - u)./dx;
        Duy = (u([2:end end],:,:) - u)./dx;
        E(num) = sum(sum(sqrt(Dux.^2 + Duy.^2))) + lambda.*sum(sum(u.*(I - c1).^2)) + lambda.*sum(sum((1 - u).*(I - c2).^2));

        %update dual variable
        Px = px + sigma.*Dux;
        Py = py + sigma.*Duy;
        px = Px./max(1,sqrt(Px.^2 + Py.^2));
        py = Py./max(1,sqrt(Px.^2 + Py.^2));
        
       %update primal variable
        Dpx = (px - px(:,[1 1:end-1],:))./dx;
        Dpy = (py - py([1 1:end-1],:,:))./dx;    
        u = max(min(u - tau.*r + tau.*(Dpx + Dpy),1),0);
        u_ = utemp;

        %visualize the contour
        if(mod(num,20) == 0)
            imshow(uint8(I.*(max(image(:)) - min(image(:))) + min(image(:)))); colormap(gray); hold on
            contour(u, [0.5 0.5], 'r');hold off
            title(sprintf('Contour at Level-Set 0 at Iteration = %d ,lambda=%d,n=%d,theta=%d',num,lambda,n,tau));
            drawnow;
        end

        num = num + 1;
        du = max(abs(u(:) - u_(:)));

    end
    
    %automatically save the energy
    name = strcat('pd_i=',num2str(itmax),',n=',num2str(n),',theta=',num2str(tau),'.mat');
    save(name,'E');
end
end
toc

% c = u.*c1 + (1 - u).*c2;
% figure
% imshow(c);


    

