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
function [Alpha_power,Alpha_score,Alpha_pos,Convergence_curve]=GWO(SearchAgents_no,Max_iter,M,N,res,req,PMAX)
lb=1;
ub=M;%number of PMs
inf_num=0;
w_num=0;
% initialize alpha, beta, and delta_pos
Alpha_pos=zeros(1,N);
Alpha_score=inf; %change this to -inf for maximization problems
Alpha_powe=inf; %change this to -inf for maximization problems


Beta_pos=zeros(1,N);
Beta_score=inf; %change this to -inf for maximization problems
Beta_power=inf; %change this to -inf for maximization problems

Delta_pos=zeros(1,N);
Delta_score=inf; %change this to -inf for maximization problems
Delta_power=inf; %change this to -inf for maximization problems

%Initialize the positions of search agents
Positions=initialization(SearchAgents_no,ub,lb,N,M,res,req);

Convergence_curve=zeros(1,Max_iter);

l=0;% Loop counter

f_VMs=zeros(SearchAgents_no,1);
w_scr=zeros(SearchAgents_no,1);
power=zeros(SearchAgents_no,1);

Alpha_power=inf;
Beta_power=inf;
Delta_power=inf;
% Main loop
while l<Max_iter
    for i=1:size(Positions,1)
        orig_norm=0;
                if (orig_norm==1)
        % Return back the search agents that go beyond the boundaries of the search space
                    Flag4ub=Positions(i,:)>ub;
                    Flag4lb=Positions(i,:)<lb;
                    Positions(i,:)=((Positions(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb);
                else
                    idx1=find(Positions(i,:)>ub);
                    idx2=find(Positions(i,:)<lb);
                    idxs=[idx1 idx2];
                    pos=Positions(i,:);
                    pos(idxs)=[];
                    PMs=available_PMS(pos,M,length(pos),res,req,0.6);
                    Positions(i,idxs)=PMs(floor(rand(1,length(idxs))*length(PMs))+1);
                end

                
                
        %         disp('dooooooo');
        % Calculate objective function for each search agent
        
        [f_VM, power(i),fitness]=cost(Positions(i,:),M,N,res,req,PMAX);
        f_VMs(i)=f_VM;
%         w_num=w_num+1;
%         if fitness==inf
%             inf_num=inf_num+1
%         end
        w_scr(i)=fitness;
        % Update Alpha, Beta, and Delta
        if fitness<=Alpha_score
            
%             Delta_power=Beta_power; % Update delta
%             Delta_score=Beta_score; % Update delta
%             Delta_pos=Beta_pos;
%             %----------------------------------
%             Beta_power=Alpha_power; % Update alpha
%             Beta_score=Alpha_score; % Update beta
%             Beta_pos=Alpha_pos;
%             %----------------------------------
            Alpha_power=power(i); % Update alpha
            Alpha_score=fitness; % Update alpha
            Alpha_pos=Positions(i,:);
        end
        
        if fitness>Alpha_score && fitness<=Beta_score
%             Delta_power=Beta_power; % Update delta
%             Delta_score=Beta_score; % Update delta
%             Delta_pos=Beta_pos;
%             %----------------------------------
            Beta_power=power(i); % Update alpha
            Beta_score=fitness; % Update beta
            Beta_pos=Positions(i,:);
            
        end
        
        if fitness>Alpha_score && fitness>Beta_score && fitness<=Delta_score
            Delta_power=power(i); % Update delta
            Delta_score=fitness; % Update delta
            Delta_pos=Positions(i,:);
            %             FPMs3=FullPMs(Delta_pos,M,N,res,req)
        end
        
        %        ************************************
    end
    for i=1:size(Positions,1)
          [Positions(i,:),power(i),fitness]=crossover(w_scr,Positions,M,N,res,req,i,power(i),PMAX);
%             w_num=w_num+1;
%             if fitness==inf
%                 inf_num=inf_num+1
%             end
        if fitness<=Alpha_score
            
%             Delta_power=Beta_power; % Update delta
%             Delta_score=Beta_score; % Update delta
%             Delta_pos=Beta_pos;
%             %----------------------------------
%             Beta_power=Alpha_power; % Update alpha
%             Beta_score=Alpha_score; % Update beta
%             Beta_pos=Alpha_pos;
%             %----------------------------------
            Alpha_power=power(i); % Update alpha
            Alpha_score=fitness; % Update alpha
            Alpha_pos=Positions(i,:);
        end
        
        if fitness>Alpha_score && fitness<=Beta_score
%             Delta_power=Beta_power; % Update delta
%             Delta_score=Beta_score; % Update delta
%             Delta_pos=Beta_pos;
%             %----------------------------------
            Beta_power=power(i); % Update alpha
            Beta_score=fitness; % Update beta
            Beta_pos=Positions(i,:);
            
        end
        
        if fitness>Alpha_score && fitness>Beta_score && fitness<=Delta_score
            Delta_power=power(i); % Update delta
            Delta_score=fitness; % Update delta
            Delta_pos=Positions(i,:);
            %             FPMs3=FullPMs(Delta_pos,M,N,res,req)
        end
    end
    
    %     f_VMs
    % % %     size(f_VMs)
    %     pause
    
       a=2-l*((2)/Max_iter); % a decreases linearly fron 2 to 0
%      a=0.3;
    % Update the Position of search agents including omegas
    for i=1:size(Positions,1)
        chs=rand();
%         f_VMs(i)
%         pause
        if(chs>1 && f_VMs(i)>0)
            L=N-f_VMs(i)+1;
%             [pms] = available_PMS(Positions(i,1:f_VMs(i)),M,f_VMs(i),res,req,0.5);
            [pms] = available_PMS(Positions(i,:),M,N,res,req,0.2);
            if (pms(1)==-1)
                pms=1:M;
                f_PM=Positions(i,f_VMs(i));
                pms(f_PM)=[];
            end
            
            Positions(i,f_VMs(i):end)=pms(floor(rand(1,(L))*length(pms))+1);
        else
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
        
    end
    l=l+1;
    Convergence_curve(l)=Alpha_score;
    figure(1);
    plot(Convergence_curve);
    %hold on;
end

%     cost2(Alpha_pos,M,N,res,req)
    disp('valuable wolves');
    inf_num/w_num
end

function [pos,power,scr]=crossover(wscr,Positions,M,N,res,req,i,pw1,PMAX)
%     w=find(wscr~=inf);
w=find(wscr~=inf);
l=length(w);
if(l>1)
    pos_t=zeros(1,N);
    r1=floor(rand()*l)+1;
    p1=Positions(w(r1),:);
    w(r1)=[];
    r2=floor(rand()*(l-1))+1;
    p2=Positions(w(r2),:);
    
    s=int32(N/2);
    
    pos_t(1:s)=p1(1:s);
    pos_t(s+1:end)=p2(s+1:end);
    
    
    [f_VM, pw2,scr2] = cost(pos_t,M,N,res,req,PMAX);
    
    if(scr2<wscr(i))
        pos=pos_t;
        power=pw2;
        scr=scr2;
    else
        pos=Positions(i,:);
        power=pw1;
        scr=wscr(i);
    end
    
else
    pos=Positions(i,:);
    power=pw1;
    scr=wscr(i);
end

end
function [PM] = available_PMS(x,M,N,res,req,ul)
res2=res;
%     w=toStruct(x,M,N);
done=0;
f_VM=0;
PM=-1;
c=1;
for n=1:N
    pm=x(n);
    res2.PM(pm).MIPS=res2.PM(pm).MIPS-req.VM(n).MIPS;
    res2.PM(pm).RAM=res2.PM(pm).RAM-req.VM(n).RAM;
    res2.PM(pm).Storage=res2.PM(pm).Storage-req.VM(n).Storage;
end
for i=1:M
    ru=1-(res.PM(i).RAM-res2.PM(i).RAM/res.PM(i).RAM);
    su=1-(res2.PM(i).Storage/res.PM(i).Storage);
    mu=1-(res2.PM(i).MIPS/res.PM(i).MIPS);
    if(ru<ul && su<ul && mu<ul)
%         disp('like---------------------');
        PM(c)=i;
        c=c+1;
    end
end
end
