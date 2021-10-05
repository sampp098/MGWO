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

% This function initialize the first population of search agents
function Positions=initialization(SearchAgents_no,ub,lb,N,M,res,req)
Boundary_no= size(ub,2); % numnber of boundaries
Boundary_no=2;
Positions=zeros(SearchAgents_no,N);
% If the boundaries of all variables are equal and user enter a signle
% number for both ub and lb
if Boundary_no==1
    disp('Boundary_no==1');
    Positions=floor(rand(SearchAgents_no,N).*M)+1;
end

% If each variable has a different lb and ub
if Boundary_no>1
    disp('Boundary_no>1');
    Positions(1,:)=FFDO(N,M,res,req);
    Positions(2,:)=FFO(N,M,res,req);
    Positions(3,:)=FFIO(N,M,res,req);
    for j=4:SearchAgents_no
        vms=1:N;
        c=1;
        for i=1:N
            r=floor(rand().*length(vms))+1;
            Positions(j,vms(r))=c;
            vms(r)=[];
            c=c+1;
            c=mod(c,M)+1;
        end
    end
%     Positions
%     pause
end
end
function position=FFDO(N,M,res,req)
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
%     if done==1
%         power=0;
%         for i=1:M
%             u=1-(res2.PM(i).MIPS/res.PM(i).MIPS);
%             if u>0
%                 power=(res.PM(i).PWI+u*(res.PM(i).PWM-res.PM(i).PWI))+power;
%             end
%         end
%         int32(power)
%         %disp(sprintf('power is:%d',power));
% %         position
%     else
%         disp('failed');
%     end

end

function position=FFO(N,M,res,req,res2)
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
% if done==1
%     power=0;
%     for i=1:M
%         u=1-(res2.PM(i).MIPS/res.PM(i).MIPS);
%         if u>0
%             power=(res.PM(i).PWI+u*(res.PM(i).PWM-res.PM(i).PWI))+power;
%         end
%     end
%     int32(power)
%     %disp(sprintf('power is:%d',power));
%     %          position
%     cost(position,M,N,res,req)
% else
%     disp('failed');
% end

end

function position=FFIO(N,M,res,req,res2)
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
% if done==1
%     power=0;
%     for i=1:M
%         u=1-(res2.PM(i).MIPS/res.PM(i).MIPS);
%         if u>0
%             power=(res.PM(i).PWI+u*(res.PM(i).PWM-res.PM(i).PWI))+power;
%         end
%     end
%     int32(power)
%     %disp(sprintf('power is:%d',power));
%     %          position
%     cost(position,M,N,res,req)
% else
%     disp('failed');
% end

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