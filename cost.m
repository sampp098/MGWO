function [f_VM, power,scr] = cost(x,M,N,res,req,PMAX)
x=int32(x);
res2=res;
%     w=toStruct(x,M,N);
done=0;
f_VM=0;
for n=1:N
    pm=x(n);
    res2.PM(pm).MIPS=res2.PM(pm).MIPS-req.VM(n).MIPS;
    res2.PM(pm).RAM=res2.PM(pm).RAM-req.VM(n).RAM;
    res2.PM(pm).Storage=res2.PM(pm).Storage-req.VM(n).Storage;
    if (res2.PM(pm).MIPS>=0)
        if (res2.PM(pm).RAM>=0)
            if (res2.PM(pm).Storage>=0)
                done =1;
            else
                done =0;
                f_VM=n;
%                 disp('st');
                break;
            end
        else
            done =0;
            f_VM=n;
            
%             disp('ram');
            break;
        end
    else
        done =0;
        f_VM=n;
%         disp('cp');
        break;
    end
end
scr=0;
if done==1
    power=0;
    for i=1:M
        u=1-(res2.PM(i).MIPS/res.PM(i).MIPS);
        ru=res2.PM(pm).RAM/res.PM(pm).RAM;
        if u>0
            p=(res.PM(i).PWI+u*(res.PM(i).PWM-res.PM(i).PWI));
            power=p+power;
            s=p*1/ru;
            scr=scr+s;
        end
    end
    
else
    scr=PMAX;
    power=PMAX;
end

end