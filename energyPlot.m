% energyPlot.m
% display loglog plot for energy functions
% visualize the result for potential energy variation
%
% Yuxin Sun
% 1498290038@qq.com
% Georgia Tech
% 2019.9.5

filename = 'bel/temp/';%file path
energyFolder = fullfile(strcat(filename));
energyInput = dir(fullfile(energyFolder,'*.mat'));
energyNum = length(energyInput);
load(strcat(filename,energyInput(1).name));
En = E;
le= '';

for i=1:energyNum;
    name = strcat(filename,energyInput(i).name);
    energy = load(name);
    energy = energy.E;
    En = [En energy];
    name = energyInput(i).name;
    name = name(1:length(name) - 4);
    title = char(length(name) + 3);
    temp = '''';
    title(1) = temp;
    title(2:length(name) + 1) = name;
    title(length(name) + 2) = temp;
    title(length(name) + 3) = ',';
    
    le= strcat(le,title);
end
En = En(:,2:energyNum + 1);
le = le(1:length(le) - 1);
figure
plot(En,'LineWidth',1);
le = str2num(le(:));
legend(le);
xlabel('CPU Time(s)');
ylabel('Potential Energy');
%legend('acceleration(k=5.5)','primaldual','quadratic(h=0.12)','quadratic(h=0.15)','quadratic(h=0.2)','quadratic(h=0.3)')
%legend('Hard-clip','Chambolle-Pock','h=0.5','h=1','h=3','h=5')