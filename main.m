%___________________________________________________________________%
%  Grey Wolf Optimizer (GWO) source codes version 1.0               %
%                                                                   %
%  Developed in MATLAB R2011b(7.13)                                 %
%                                                                   %
%  Author and programmer: Seyedali Mirjalili                        %
%                                                                   %
%         e-Mail: ali.mirjalili@gmail.com                           %
%                 seyedali.mirjalili@griffithuni.edu.au             %
%                                                                   %
%       Homepage: http://www.alimirjalili.com                       %
%                                                                   %
%   Main paper: S. Mirjalili, S. M. Mirjalili, A. Lewis             %
%               Grey Wolf Optimizer, Advances in Engineering        %
%               Software , in press,                                %
%               DOI: 10.1016/j.advengsoft.2013.12.007               %
%                                                                   %
%___________________________________________________________________%

% You can simply define your cost in a seperate file and load its handle to fobj 
% The initial parameters that you need are:
%__________________________________________
% fobj = @YourCostFunction
% dim = number of your variables
% Max_iteration = maximum number of generations
% SearchAgents_no = number of search agents
% lb=[lb1,lb2,...,lbn] where lbn is the lower bound of variable n
% ub=[ub1,ub2,...,ubn] where ubn is the upper bound of variable n
% If all the variables have equal lower bound you can just
% define lb and ub as two single number numbers

% To run GWO: [Best_score,Best_pos,GWO_cg_curve]=GWO(SearchAgents_no,Max_iteration,lb,ub,dim,fobj)
%__________________________________________
function main()
clear all 
clc
SearchAgents_no=20; % Number of search agents

% Name of the test function that can be from F1 to F23 (Table 1,2,3 in the paper)

Max_iteration=1000; % Maximum numbef of iterations

% Load details of the selected benchmark function
[fobj]=Get_Functions_details('cost');
M=20;%number of PMs
N=7;%number of VMs
res=get_res(M);
req=get_req(N);
tic
[Best_score,Best_pos,GWO_cg_curve]=GWO(SearchAgents_no,Max_iteration,fobj,M,N,res,req);
toc
% figure('Position',[500 500 660 290])
% %Draw search space
% subplot(1,2,1);
% func_plot(Function_name);
% title('Parameter space')
% xlabel('x_1');
% ylabel('x_2');
% zlabel([Function_name,'( x_1 , x_2 )'])

%Draw objective space
% subplot(1,2,2);
% semilogy(GWO_cg_curve,'Color','r')
% title('Objective space')
% xlabel('Iteration');
% ylabel('Best score obtained so far');

% axis tight
% grid on
% box on
% legend('GWO')

display(['The best solution obtained by GWO is : ', num2str(Best_pos)]);
display(['The best optimal value of the objective funciton found by GWO is : ', num2str(Best_score)]);
        
end
function req=get_req(N)
VM_RAM=[1, 2, 3, 4]; %GB
VM_BW=[100, 1000]; %Mb/s
VM_MIPS=[5000, 10000, 2*10000, 3*10000];
VM_Storage=[10, 20, 30, 40, 50]; %GB
VM_PE=[1, 2, 3, 4];
VM_Live_Time=[1, 3, 6, 12, 24, 48, 36, 72];%hours
for VMi=1:N
    req.VM(VMi).RAM=VM_RAM(floor(rand()*length(VM_RAM))+1);
    req.VM(VMi).BW=VM_BW(floor(rand()*length(VM_BW))+1);
    req.VM(VMi).MIPS=VM_MIPS(floor(rand()*length(VM_MIPS))+1);
    req.VM(VMi).Storage=VM_Storage(floor(rand()*length(VM_Storage))+1);
    req.VM(VMi).PE=VM_PE(floor(rand()*length(VM_PE))+1);
    req.VM(VMi).VM_LT=VM_Live_Time(floor(rand()*length(VM_Live_Time))+1);
end
end
function res=get_res(M)
PM_RAM=[4, 8, 16, 32]; %GB
PM_BW=[1000, 2000]; %Mb/s
PM_MIPS=[5000, 10000, 2*10000, 3*10000];
PM_Storage=[100, 200, 300, 400, 500,1000]; %GB
PM_PE=[8, 16, 32]; %GB

PM_PWI=[58.4,66,42.3,93.7,86,105,41.6];
PM_PWM=[222,247,113,135,117,169,113];

for pmi=1:M
    res.PM(pmi).RAM=PM_RAM(floor(rand()*length(PM_RAM))+1);
    res.PM(pmi).BW=PM_BW(floor(rand()*length(PM_BW))+1);
    res.PM(pmi).MIPS=PM_MIPS(floor(rand()*length(PM_MIPS))+1);
    res.PM(pmi).Storage=PM_Storage(floor(rand()*length(PM_Storage))+1);
    res.PM(pmi).PE=PM_PE(floor(rand()*length(PM_PE))+1);
    r=floor(rand()*length(PM_PWI))+1;
    res.PM(pmi).PWI=PM_PWI(r);
    res.PM(pmi).PWM=PM_PWM(r);
end
end



