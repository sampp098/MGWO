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

% Grey Wolf Optimizer
function [Alpha_score,Alpha_pos,Convergence_curve]=GWO(SearchAgents_no,Max_iter,fobj,M,N,res,req)
lb=1;
ub=M;%number of PMs

% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,N);
Alpha_score=inf; %change this to -inf for maximization problems

Beta_pos=zeros(1,N);
Beta_score=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,N);
Delta_score=inf; %change this to -inf for maximization problems

%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,ub,lb,N,M,res,req);

Convergence_curve=zeros(1,Max_iter);

l=0;% Loop counter



% Main loop
while l<Max_iter
    for i=1:size(Positions,1)
        
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=Positions(i,:)>ub;
        Flag4lb=Positions(i,:)<lb;
        Positions(i,:)=((Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb);
        %         disp('dooooooo');
        % Calculate objective function for each search agent
        fitness=fobj(Positions(i,:),M,N,res,req);
        % Update Alpha, Beta, and Delta
        if fitness<Alpha_score
            Alpha_score=fitness; % Update alpha
            Alpha_pos=Positions(i,:);
%             FPMs1=FullPMs(Alpha_pos,M,N,res,req)
            if l==1
                Alpha_score
                Beta_score=Alpha_score % init beta
                Beta_pos=Alpha_pos;
                Delta_score=Alpha_score % init delta
                Delta_pos=Alpha_pos;
            end
        end
        
        if fitness>Alpha_score && fitness<Beta_score
            Beta_score=fitness; % Update beta
            Beta_pos=Positions(i,:);
%             FPMs2=FullPMs(Beta_pos,M,N,res,req)
        end
        
        if fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score
            Delta_score=fitness; % Update delta
            Delta_pos=Positions(i,:);
%             FPMs3=FullPMs(Delta_pos,M,N,res,req)
        end
    end
    
    
    a=2-l*((2)/Max_iter); % a decreases linearly fron 2 to 0
    
    % Update the Position of search agents including omegas
    for i=1:size(Positions,1)
        for j=1:size(Positions,2)
            
            r1=rand(); % r1 is a random number in [0,1]
            r2=rand(); % r2 is a random number in [0,1]
            
            A1=2*a*r1-a; % Equation (3.3)
            C1=2*r2; % Equation (3.4)
            
            D_alpha=abs(C1*Alpha_pos(j)-Positions(i,j)); % Equation (3.5)-part 1
            X1=Alpha_pos(j)-A1*D_alpha; % Equation (3.6)-part 1
            
            r1=rand();
            r2=rand();
            
            A2=2*a*r1-a; % Equation (3.3)
            C2=2*r2; % Equation (3.4)
            
            D_beta=abs(C2*Beta_pos(j)-Positions(i,j)); % Equation (3.5)-part 2
            X2=Beta_pos(j)-A2*D_beta; % Equation (3.6)-part 2
            
            r1=rand();
            r2=rand();
            
            A3=2*a*r1-a; % Equation (3.3)
            C3=2*r2; % Equation (3.4)
            
            D_delta=abs(C3*Delta_pos(j)-Positions(i,j)); % Equation (3.5)-part 3
            X3=Delta_pos(j)-A3*D_delta; % Equation (3.5)-part 3
            
            Positions(i,j)=floor((X1+X2+X3)/3)+1;% Equation (3.7)
            %             disp('average');
        end
    end
    l=l+1;
    Convergence_curve(l)=Alpha_score;
    figure(1);
    plot(Convergence_curve);
    %hold on;
end
    cost2(Alpha_pos,M,N,res,req)
end
function FPMs=FullPMs(x,M,N,res,req)
w=toStruct(x,M,N);
c=1;
FPMs=[];

for i=1:M
    VM_num=w.PM(i).VMs(1);
    if(VM_num>0)
        PM=res.PM(i);
        VMs=req.VM(w.PM(i).VMs(2:end));
        u=get_PM_MIPS_used(VMs,PM);
        if u<0.95 &&  u>0.85
            FPMs(c)=i;
            c=c+1;
        end
    end
end

end
function power = cost2(x,M,N,res,req)
    w=toStruct(x,M,N);
    power=0;
    for i=1:M
       VM_num=w.PM(i).VMs(1);
       if(VM_num>0)
           PM=res.PM(i);
           VMs=req.VM(w.PM(i).VMs(2:end));
           
           %RAM CHECK
           ram_p=get_PM_RAM_used(VMs,PM)
           
           %BW CHECK
%            bw_p=get_PM_BW_used(VMs,PM)           
           %Storage CHECK
           Storage_p=get_PM_Storage_used(VMs,PM)           
           %MIPS CHECK
           MIPS_p=get_PM_MIPS_used(VMs,PM)
       end
    end

end

function RAM=get_PM_RAM_used(VMs,PM)
    S=0;
    for i=1:length(VMs)
        S=S+VMs(i).RAM;
    end
    RAM=S/PM.RAM;
end
% function BW=get_PM_BW_used(VMs,PM)
%     S=0;
%     for i=1:length(VMs)
%         S=S+VMs(i).BW;
%     end
%     BW=S/PM.BW;
% end
function Storage=get_PM_Storage_used(VMs,PM)
    S=0;
    for i=1:length(VMs)
        S=S+VMs(i).Storage;
    end
    Storage=S/PM.Storage;
end
function MIPS=get_PM_MIPS_used(VMs,PM)
    S=0;
    for i=1:length(VMs)
        S=S+(VMs(i).MIPS*VMs(i).PE);
    end
    MIPS=S/(PM.MIPS*PM.PE);
end
function w=toStruct(x,M,N)
w.PM=[];
for m=1:M
    w.PM(m).VMs(1)=0;
end
for vm=1:N
    m=x(vm);
    w.PM(m).VMs(1)=w.PM(m).VMs(1)+1;
    w.PM(m).VMs(w.PM(m).VMs(1)+1)=vm;
end
end