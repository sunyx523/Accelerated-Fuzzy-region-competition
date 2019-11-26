% main_pd_bel.m
% Chambolle-Pock primal-dual for Beltrami regularization
% The segmentation scheme we used here is called fuzzy region competition from:
% Mory-Ardon, 2007, "Fuzzy Region Competition: A Convex Two-Phase Segmentation Framework"
% The convex energy is minimized for instance by using Beltrami version Chambolle-Pock  primal-dual algorithm
% Zosso-Bustin, 2014, "A Primal-Dual Projected Gradient Algorithm for Efficient Beltrami Regularization"
%
% getBeta(): An alternative choice for finding best beta for each
%            iteration, but computational cost will increase.
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


%parameters
lambda = 0.05;      %Fidelity weight
dx = 1;             %Grid point distance
itmax = 5000;      %Maximum iterations
dumax = 0.001;
L2 = 8;


tic
for b = [0.5]          %beta in beltrami regularization

    figure
    sigma = sqrt(1./L2/b); %sigma parameter from the paper
    tau = sqrt(1./L2/b);
    E = zeros(itmax,1);  %Initla potential energy function

    u = u0;              %u for level set(primal variable)
    u_ = u;              %u_ for previous level set
    num = 1;             %Iteration index
    du = 100;
    
%     while  (du > dumax && num < itmax)
    while   num < itmax
        
        beta = b;
        utemp = u;
        %beta = getBeta(u,dx,b);   
        c1 = sum(sum(u.*I))/sum(sum(u));
        c2 = sum(sum((1 - u).*I))/sum(sum(1 - u));
        r = lambda.*((I - c1).^2 - (I - c2).^2);

        Dux = (u(:,[2:end end],:) - u)./dx;
        Duy = (u([2:end end],:,:) - u)./dx;
        E(num) = 1./beta.*sum(sum(sqrt(1 + beta.^2.*(Dux.^2 + Duy.^2)))) + lambda.*sum(sum(u.*(I - c1).^2)) + lambda.*sum(sum((1 - u).*(I - c2).^2));
        
        %update dual variable
        Den = real(sqrt(beta.^2 - px.^2 - py.^2));
        Den(Den < 0) = 0;
        Px = px - sigma.*px + sigma.*Dux.*Den;
        Py = py - sigma.*py + sigma.*Duy.*Den;
        px = beta.*Px./max(beta,sqrt(Px.^2 + Py.^2));
        py = beta.*Py./max(beta,sqrt(Px.^2 + Py.^2));
        
        %update primal variable
        Dpx = (px - px(:,[1 1:end-1],:))./dx;
        Dpy = (py - py([1 1:end-1],:,:))./dx;    
        u = max(min(u - tau.*r + tau.*(Dpx + Dpy),1),0);
        u_ = utemp;

        %visualize the contour
        if(mod(num,100) == 0)
            imshow(uint8(I.*(max(image(:)) - min(image(:))) + min(image(:)))); colormap(gray); hold on
            contour(u, [0.5 0.5], 'r');hold off
            title(sprintf('Contour at Level-Set 0 at Iteration = %d', num));
            drawnow;
        end

        num = num + 1;
        du = max(abs(u(:) - u_(:)));

    end
    
    %automatically save the energy

    name = strcat('pd_b =',num2str(b),',i=',num2str(tmax),',tau=',num2str(sigma),',theta=',num2str(tau),'.mat');
    save(name,'E');
end
toc

% c = u.*c1 + (1 - u).*c2;
% figure
% imshow(c);
