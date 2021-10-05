function FFD()
clear all
clc

% SearchAgents_no=20; % Number of search agents
% Max_iteration=1000; % Maximum numbef of iterations
% [fobj]=Get_Functions_details('cost');
gwo_FFD=1;
offline=1;
str='data\30\30-100\';
SearchAgents_no=60;
Max_iteration=1000;
if offline==0
    M=40;%number of PMs
    N=40;%number of VMs
    res=get_res(M);
    req=get_req(N);
    res.storage=0;
    res.ram=0;
    res.mips=0;
    for i=1:M
        res.PM(i).MIPS=res.PM(i).MIPS*res.PM(i).PE;
        res.PM(i).PE=1;
        
        res.mips=res.mips+res.PM(i).MIPS;
        res.storage=res.storage+res.PM(i).Storage;
        res.ram=res.ram+res.PM(i).RAM;
        
    end
    req.storage=0;
    req.ram=0;
    req.mips=0;
    for i=1:N
        req.VM(i).MIPS=req.VM(i).MIPS*req.VM(i).PE;
        req.VM(i).PE=1;
        
        req.mips=req.mips+req.VM(i).MIPS;
        req.storage=req.storage+req.VM(i).Storage;
        req.ram=req.ram+req.VM(i).RAM;
    end
    
    save('res.mat','res');
    save('req.mat','req');
elseif offline==1
    load(strcat(str,'res.mat'));
    load(strcat(str,'req.mat'));
    N=length(req.VM);
    M=length(res.PM);
end
PMAX=0;
for i=1:M
    PMAX=PMAX+res.PM(i).PWM;
end
PMAX=PMAX*2;
if gwo_FFD==1
    FFDO(N,M,res,req,PMAX);
    FFO(N,M,res,req,PMAX);
    FFIO(N,M,res,req,PMAX);
    
    GA(SearchAgents_no,Max_iteration,M,N,res,req,PMAX);
else
    
    tic
    [Best_power,Best_score,Best_pos,GWO_cg_curve]=GWO(SearchAgents_no,Max_iteration,M,N,res,req,PMAX);
    toc
    display(['The best solution obtained by GWO is : ', num2str(Best_pos)]);
    display(['The best optimal score of the objective funciton found by GWO is : ', num2str(Best_score)]);
    display(['The best optimal power of the objective funciton found by GWO is : ', num2str(Best_power)]);
end

end
function fit=isFit(PM,VM)
fit=0;
if (PM.RAM-VM.RAM>0)
    if (PM.MIPS-VM.MIPS>0)
        if (PM.Storage-VM.Storage>0)
            fit= 1;
        end
    end
end
end
function req=get_req(N)
VM_RAM=[2,2,4,4,4,4,8,8,8,8,16,16]; %GB
% VM_BW=[100, 1000]; %Mb/s
VM_MIPS=[40000,40000,40000,40000,40000,40000,40000,40000,40000,40000,40000,40000];
VM_Storage=[10,40,20,40,20,40,40,60,50,100,50,100]; %GB
VM_PE=[1,1,2,2,4,4,4,4,8,8,8,8];
VM_Live_Time=[1, 3, 6, 12, 24, 48, 36, 72];%hours

for VMi=1:N
    r=floor(rand()*length(VM_RAM))+1;scr=r;
    req.VM(VMi).RAM=VM_RAM(r);
    %     r=floor(rand()*length(VM_BW))+1;scr=scr+r;
    %     req.VM(VMi).BW=VM_BW(r);
    scr=scr+r;
    req.VM(VMi).MIPS=VM_MIPS(r);
    scr=scr+r;
    req.VM(VMi).Storage=VM_Storage(r);
    scr=scr+r;
    req.VM(VMi).PE=VM_PE(r);
    
%     r=floor(rand()*length(VM_Live_Time))+1;
%     req.VM(VMi).LT=VM_Live_Time(r);
    req.VM(VMi).SCR=scr;
end
end
function res=get_res(M)
PM_RAM=[32,32,64,128,256]; %GB
% PM_BW=[1000, 2000]; %Mb/s
PM_MIPS=[40000,40000,40000,40000,40000];
PM_Storage=[1000,2000,5000,5000,10000]; %GB
PM_PE=[16, 36, 64, 72, 72];

PM_PWI=[20,30,40,60,60];
PM_PWM=[130,195,260,390,390];

for pmi=1:M
    r=floor(rand()*length(PM_RAM))+1;
    scr=r;
    res.PM(pmi).RAM=PM_RAM(r);
    %     r=floor(rand()*length(PM_BW))+1;scr=r+scr;
    %     res.PM(pmi).BW=PM_BW(r);
    scr=r+scr;
    res.PM(pmi).MIPS=PM_MIPS(r);
    scr=r+scr;
    res.PM(pmi).Storage=PM_Storage(r);
    scr=r+scr;
    res.PM(pmi).PE=PM_PE(r);
    res.PM(pmi).PWI=PM_PWI(r);
    res.PM(pmi).PWM=PM_PWM(r);
    
    res.PM(pmi).SCR=scr;
end
end
function position=FFDO(N,M,res,req,PMAX)
%sort
VM_SCRs=[req.VM(:).SCR];
[B1,sVM]=sort(VM_SCRs);

PM_SCRs=[res.PM(:).SCR];
[B2,sPM]=sort(PM_SCRs);


res2=res;
%placement
position=zeros(1,N);
done=0;
for i=1:N
    for j=1:M
        done=0;
        if isFit(res2.PM(sPM(j)),req.VM(sVM(i)))==1
            position(sVM(i))=sPM(j);
            res2.PM(sPM(j)).MIPS=res2.PM(sPM(j)).MIPS-req.VM(sVM(i)).MIPS;
            res2.PM(sPM(j)).RAM=res2.PM(sPM(j)).RAM-req.VM(sVM(i)).RAM;
            res2.PM(sPM(j)).Storage=res2.PM(sPM(j)).Storage-req.VM(sVM(i)).Storage;
            %                 res2.PM(sPM(j)).BW=res2.PM(sPM(j)).BW-req.VM(sVM(i)).BW;
            done=1;
            break;
        end
    end
    if done==0;
        break;
    end
end
if done==1
    power=0;
    for i=1:M
        u=1-(res2.PM(i).MIPS/res.PM(i).MIPS);
        if u>0
            power=(res.PM(i).PWI+u*(res.PM(i).PWM-res.PM(i).PWI))+power;
        end
    end
%     int32(power)
    %disp(sprintf('power is:%d',power));
    %          position
    disp('*******************************************');
    [f_VM, power,scr]=cost(position,M,N,res,req,PMAX);
    display(['The best solution obtained by FFD is : ', num2str(position)]);
    display(['The best optimal score of the objective funciton found by FFD is : ', num2str(scr)]);
    display(['The best optimal power of the objective funciton found by FFD is : ', num2str(power)]);
else
    disp('failed');
end

end

function position=FFO(N,M,res,req,PMAX)
%sort
% VM_SCRs=[req.VM(:).SCR];
% [B1,sVM]=sort(VM_SCRs);
% 
% PM_SCRs=[res.PM(:).SCR];
% [B2,sPM]=sort(PM_SCRs);


res2=res;
%placement
position=zeros(1,N);
done=0;
for i=1:N
    for j=1:M
        done=0;
        if isFit(res2.PM(j),req.VM(i))==1
            position(i)=j;
            res2.PM(j).MIPS=res2.PM(j).MIPS-req.VM(i).MIPS;
            res2.PM(j).RAM=res2.PM(j).RAM-req.VM(i).RAM;
            res2.PM(j).Storage=res2.PM(j).Storage-req.VM(i).Storage;
            %                 res2.PM(sPM(j)).BW=res2.PM(sPM(j)).BW-req.VM(sVM(i)).BW;
            done=1;
            break;
        end
    end
    if done==0;
        break;
    end
end
if done==1
    power=0;
    for i=1:M
        u=1-(res2.PM(i).MIPS/res.PM(i).MIPS);
        if u>0
            power=(res.PM(i).PWI+u*(res.PM(i).PWM-res.PM(i).PWI))+power;
        end
    end
    %int32(power)
    %disp(sprintf('power is:%d',power));
    %          position
    disp('*******************************************');
    [f_VM, power,scr]=cost(position,M,N,res,req,PMAX);
    save('position.mat','position');
    display(['The best solution obtained by FF is : ', num2str(position)]);
    display(['The best optimal score of the objective funciton found by FF is : ', num2str(scr)]);
    display(['The best optimal power of the objective funciton found by FF is : ', num2str(power)]);
else
    disp('failed');
end

end

function position=FFIO(N,M,res,req,PMAX)
%sort
VM_SCRs=[req.VM(:).SCR];
[B1,sVM]=sort(VM_SCRs);

PM_SCRs=[res.PM(:).SCR];
[B2,sPM]=sort(PM_SCRs);
sPM_t=sPM(end:-1:1);
sPM=sPM_t;

res2=res;
%placement
position=zeros(1,N);
done=0;
for i=1:N
    for j=1:M
        done=0;
        if isFit(res2.PM(sPM(j)),req.VM(sVM(i)))==1
            position(sVM(i))=sPM(j);
            res2.PM(sPM(j)).MIPS=res2.PM(sPM(j)).MIPS-req.VM(sVM(i)).MIPS;
            res2.PM(sPM(j)).RAM=res2.PM(sPM(j)).RAM-req.VM(sVM(i)).RAM;
            res2.PM(sPM(j)).Storage=res2.PM(sPM(j)).Storage-req.VM(sVM(i)).Storage;
            %                 res2.PM(sPM(j)).BW=res2.PM(sPM(j)).BW-req.VM(sVM(i)).BW;
            done=1;
            break;
        end
    end
    if done==0;
        break;
    end
end
if done==1
    power=0;
    for i=1:M
        u=1-(res2.PM(i).MIPS/res.PM(i).MIPS);
        if u>0
            power=(res.PM(i).PWI+u*(res.PM(i).PWM-res.PM(i).PWI))+power;
        end
    end
%     int32(power)
    %disp(sprintf('power is:%d',power));
    %          position
    disp('*******************************************');
    [f_VM, power,scr]=cost(position,M,N,res,req,PMAX);
    display(['The best solution obtained by FFI is : ', num2str(position)]);
    display(['The best optimal score of the objective funciton found by FFI is : ', num2str(scr)]);
    display(['The best optimal power of the objective funciton found by FFI is : ', num2str(power)]);
else
    disp('failed');
end

end

function    crosspop=crossover(crosspop,pop,ncross,M,N,res,req,PMAX)

f=[pop.scr];
f=f./sum(f);
f=cumsum(f);


for n=1:2:ncross

    
    i1=find(rand<=f,1,'first');
    i2=find(rand<=f,1,'first');
    
    
p1=pop(i1).par;
p2=pop(i2).par;



R=rand(size(p1));

o1=floor((p1.*R)+(p2.*(1-R)))+1;
%SAM..back to the range
Flag4ub=o1(:)>M;
Flag4lb=o1(:)<1;
o1(:)=((o1(:).*(~(Flag4ub+Flag4lb)))+M.*Flag4ub+1.*Flag4lb);

o2=floor((p2.*R)+(p1.*(1-R)))+1;
%SAM..back to the range
Flag4ub=o2(:)>M;
Flag4lb=o2(:)<1;
o2(:)=((o2(:).*(~(Flag4ub+Flag4lb)))+M.*Flag4ub+1.*Flag4lb);

crosspop(n).par=o1;%[vm,pop(i).scr,pop(i).scr]=cost(pop(i).par,M,N,res,req)
[vm,crosspop(n).power,crosspop(n).scr]=cost(o1,M,N,res,req,PMAX);


crosspop(n+1).par=o2;
%crosspop(n+1).scr=fitness(o2);
[vm,crosspop(n+1).power,crosspop(n+1).scr]=cost(o2,M,N,res,req,PMAX);
end

end

function  mutpop=mutation(mutpop,pop,nmut,npop,lb,ub,nvar,M,N,res,req,PMAX)

for n=1:nmut

i=randi([1 npop]);     
sol=pop(i).par;




j=randi([1 nvar]);

d=0.1*unifrnd(-1,1)*(ub(j)-lb(j));


sol(j)=sol(j)+d;

sol=min(sol,ub);
sol=max(sol,lb);



mutpop(n).par=sol;
[vm,mutpop(n).power,mutpop(n).scr]=cost(sol,M,N,res,req,PMAX);

end



end
function GA(SearchAgents_no,Max_iteration,M,N,res,req,PMAX)

%% paramters setting



nvar=N;    % number of variable



npop=SearchAgents_no;         % number of population

maxiter=Max_iteration;      % max of iteration


pc=0.1;                  % percent of crossover
ncross=2*round(npop*pc/2);   % number of cross over offspring

pm=1-pc;                 % percent of mutation
nmut=round(npop*pm);     % number of mutation offsprig

lb=ones(1,nvar); 
ub=M*ones(1,nvar); 


%% initialization
tic
empty.par=[];
empty.scr=[];
empty.power=[];

pop=repmat(empty,npop,1);


for i=1:npop
    
   pop(i).par=floor(lb+rand(1,nvar).*ub);
   [vm,pop(i).power,pop(i).scr]=cost(pop(i).par,M,N,res,req,PMAX);
   
end




%% main loop


BEST=zeros(maxiter,1);
MEAN=zeros(maxiter,1);
maxiter
for iter=1:maxiter


   % crossover
   
   crosspop=repmat(empty,ncross,1);
   crosspop=crossover(crosspop,pop,ncross,M,N,res,req,PMAX);
   
   
   
   
   % mutation
   mutpop=repmat(empty,nmut,1);
   mutpop=mutation(mutpop,pop,nmut,npop,lb,ub,nvar,M,N,res,req,PMAX);
   
   
   % merged
  [pop]=[pop;crosspop;mutpop];


  % select
  [value,index]=sort([pop.scr]);
  pop=pop(index);
  pop=pop(1:npop);

 gpop=pop(1);   % global pop



 BEST(iter)=gpop.scr;
 MEAN(iter)=mean([pop.scr]);



%  disp([ ' Iter = '  num2str(iter)  ' BEST = '  num2str(BEST(iter))]);


  if iter>400 && BEST(iter)==BEST(iter-400)
      break;
  end

end
%% results

%disp(' ')
%disp([ ' Best par = '  num2str(int32(gpop.par))])
%disp([ ' Best fitness = '  num2str(gpop.scr)])
%disp([ ' Time = '  num2str(toc)])

disp('*******************************************');
display(['The best solution obtained by GA is : ' num2str(int32(gpop.par))]);
display(['The best optimal score of the objective funciton found by GA is : ' num2str(gpop.scr)]);
display(['The best optimal power of the objective funciton found by GA is : ' num2str(gpop.power)]);

score=zeros(1,maxiter);
score(:)=BEST(iter);
score(1:iter)=BEST(1:iter);
figure(1)
plot(score,'b','LineWidth',2)
hold on
%plot(MEAN(1:iter),'b','LineWidth',2)


xlabel('Iteration')
ylabel(' Fitness')

%legend('BEST','MEAN')

%title('GA for TSP')
end

%%

function RAM=get_PM_RAM_used(VMs,PM)
S=0;
for i=1:length(VMs)
    S=S+VMs(i).RAM;
end
RAM=S/PM.RAM;
end
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


