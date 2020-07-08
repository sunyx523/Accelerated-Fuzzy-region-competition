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

px = zeros(m,n); %Initial projector(dual variable) on x direction
py = zeros(m,n); %Initial projector(dual variable) on y dirextion

lambda = 500;      %Fidelity weight
dx = 1/m;          %Grid point distance
tmax = 1000;       %Maximum time
ind = 1;           %Iteration index
dumax = 0.01;
L2 = 6;

tic
for tau = [1]            
    
    figure
    a = 1;
    sigma = dx*sqrt(a/L2); %sigma parameter from the paper
    tau = dx*sqrt(1/a/L2);
    E = zeros(tmax,1);   %Initla potential energy function
    DU = zeros(tmax,1);

    u = u0;               %u for level set(primal variable)
    u_ = u;               %u_ for previous level set
    ind = 0;
    du = 100;
    t1 = clock;
    
%    while  (du > dumax && ind < tmax)
    while ind < tmax/10
        
        utemp = u;
        c1 = sum(sum(u.*I))/sum(sum(u));
        c2 = sum(sum((1 - u).*I))/sum(sum(1 - u));
        r = lambda.*((I - c1).^2 - (I - c2).^2);

        Dux = (u(:,[2:end end],:) - u)./dx;
        Duy = (u([2:end end],:,:) - u)./dx;
        e = sum(sum(sqrt(Dux.^2 + Duy.^2))) + lambda.*sum(sum(u.*(I - c1).^2)) + lambda.*sum(sum((1 - u).*(I - c2).^2));
        t2 = clock;
        
        %visualize the contour
        if(etime(t2,t1) > ind)
            E(int32(ind.*10 + 1)) = e;
            DU(int32(ind.*10 + 1)) = du;
            ind = ind + 0.1;
        end
        
%        visualize contour
%         if(mod(int32(ind*10),10) == 0)
%             imshow(uint8(I.*(max(image(:)) - min(image(:))) + min(image(:)))); colormap(gray); hold on
%             contour(u, [0.5 0.5], 'r');hold off
%             title(sprintf('Primal dual TV'));
%             drawnow;
%         end    

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

        du = max(abs(u(:) - u_(:)));

    end
    
    %automatically save the energy
    name = strcat('pd_t=',num2str(tmax),',sigma=',num2str(sigma),',tau=',num2str(tau),'.mat');
    save(name,'E');
end
toc

% c = u.*c1 + (1 - u).*c2;
% figure
% imshow(c);


    

