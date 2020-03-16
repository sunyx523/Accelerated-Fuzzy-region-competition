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
tmax = 1000;      %Maximum iterations
L2 = 6;
g = 1;

tic
for b = [2]          %beta in beltrami regularization
    
    figure
    a = 2.*sqrt(b.*g.*(pi.^2)/m/n);
    sigma = sqrt(a./L2/b); %sigma parameter from the paper
    tau = sqrt(1/L2/a/b);
    E = zeros(tmax,1);     %Initla potential energy function
    u = u0;                %u for level set(primal variable)
    u_ = u;                %u_ for previous level set
    ind = 0;
    t1 = clock;
 
    while (ind < tmax/10)        
        beta = b;
        utemp = u;
        %beta = getBeta(u,dx,b);   
        c1 = sum(sum(u.*I))/sum(sum(u));
        c2 = sum(sum((1 - u).*I))/sum(sum(1 - u));
        r = lambda.*((I - c1).^2 - (I - c2).^2);

        Dux = (u(:,[2:end end],:) - u)./dx;
        Duy = (u([2:end end],:,:) - u)./dx;
        e = 1./beta.*sum(sum(sqrt(1 + beta.^2.*(Dux.^2 + Duy.^2)))) + lambda.*sum(sum(u.*(I - c1).^2)) + lambda.*sum(sum((1 - u).*(I - c2).^2));
        t2 = clock;
        
        if(etime(t2,t1) > ind)
            E(int32(ind.*10 + 1)) = e;
            ind = ind + 0.1;
        end
        
%        visualize contour
%         if(mod(int32(ind*10),5) == 0)
% %             imshow(uint8(I.*(max(image(:)) - min(image(:))) + min(image(:)))); colormap(gray); hold on
% %             contour(u, [0.5 0.5], 'r');hold off
% %             title(sprintf('Primal dual Beltrami'));
%             imshow(u); hold on
%             drawnow;
%         end
        
        %update dual variable
        Den = real(sqrt(beta.^2 - px.^2 - py.^2));
        Den(Den < 0) = 0;
        Px = px - sigma.*beta.*px + beta.*sigma.*Dux.*Den;
        Py = py - sigma.*beta.*py + beta.*sigma.*Duy.*Den;
        px = beta.*Px./max(beta,sqrt(Px.^2 + Py.^2));
        py = beta.*Py./max(beta,sqrt(Px.^2 + Py.^2));
        
        %update primal variable
        Dpx = (px - px(:,[1 1:end-1],:))./dx;
        Dpy = (py - py([1 1:end-1],:,:))./dx;    
        u = max(min(u - tau.*r + tau.*(Dpx + Dpy),1),0);
        u_ = utemp;
        
    end
    
    %automatically save the energy

    name = strcat('t_pd_b =',num2str(b),',t=',num2str(tmax),',sigma=',num2str(sigma),',tau=',num2str(tau),'.mat');
    save(name,'E');
end
toc

% c = u.*c1 + (1 - u).*c2;
% figure
% imshow(c);
