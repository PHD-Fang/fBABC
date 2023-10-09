% function [Res,cg_curve1] = BRO_Fun(N,Dim,maxiter,model)
function  [minIter,notUpdateTimes,GAP, Convergence , Best_Obj , Best_Fit ,Best_Sol] = BRO_Fun(N,Dim,maxiter,model,initPos,CrossType)
global NFE;
MaxFault = 4;
thetama = 1;
thetami = 0.1;
% CrossType = 3;
Res = [];
NumRun = 1;
tic
%% 
% cg_curve1 = zeros(1,maxiter);
Hit = 0;
Gap = 0;
EVL = [];
ListNFE = [];
NFE = 0;
R=1;
Shrink = ceil(log10(maxiter));
Step = round(maxiter/Shrink);
%%
soldier.xy = [];
soldier.Fit = [];
soldier.Fault  = zeros(1,N);
flag.Fit = 0;
soldier.xy = [];
flag.xy = [];
Eval = @MaxCut;
%     ShrinkMinMax = XMinMax;

%% Initialization
% soldier.xy = initPos;%randi([0 1],N,Dim);

for i=1:N
    soldier.xy(i,:) = initPos(i,:);
    soldier.Fit(i) = Eval(model,soldier.xy(i,:));
end
initFit = soldier.Fit;

[flag.Fit,indx] = max(soldier.Fit);
flag.xy = soldier.xy(indx,:);

%%
STDdim = zeros(Dim,1);
cg_curve = zeros(1,maxiter);
cg_curve(1) = flag.Fit;
for iter = 2:maxiter
    %         theta = thetama-((thetama-thetami)/maxiter)*iter;
    theta = Dim-((Dim-1)/maxiter)*iter;
    a = theta/Dim;
    
    for i = 1:N
        [~,Cn] =  Dsimilarity(soldier.xy(i,:),soldier.xy,i);
        dam = i;
        Vic = Cn;
        if (soldier.Fit(Cn) < soldier.Fit(i))
            dam = Cn;
            Vic = i;
        end

        if (soldier.Fault(dam) < MaxFault)
%                 [soldier.xy(dam,:) soldier.Fit(dam)]= MyCrossOverFcn(soldier.xy(dam,:),flag.xy,Dim,CrossType,Eval,model);
            [soldier.xy(dam,:)]= MyCrossOverFcn(soldier.xy(dam,:),flag.xy,Dim,CrossType,Eval,model);
            %%%%% mutation for Dam
            if rand<0.3
                [soldier.xy(dam,:)] = mutation(soldier.xy(dam,:),soldier.Fit(dam),Eval,model);
            end
            %%%%%%%%%%%%%%%%%
            soldier.Fit(dam) =  Eval(model,soldier.xy(dam,:));
            soldier.Fault(dam)= soldier.Fault(dam) + 1;
            soldier.Fault(Vic) = 0;
        else
            b = randperm(Dim,ceil(a*Dim));
            soldier.xy(dam,b) = randi([0 1],1,numel(b));
            soldier.Fault(dam) = 0;
            soldier.Fit(dam) = Eval(model,soldier.xy(dam,:));
        end
%             [soldier.xy(Vic,:),soldier.Fit(Vic)] = mutation(soldier.xy(Vic,:),soldier.Fit(Vic),Eval,model);
            
            %%%%% mutation for Vic
            if rand<0.3
                [soldier.xy(Vic,:)] = mutation(soldier.xy(Vic,:),soldier.Fit(Vic),Eval,model);
            end
            soldier.Fit(Vic) =  Eval(model,soldier.xy(Vic,:));
%             [flag.xy,flag.Fit] = mutation(flag.xy,flag.Fit,Eval,model);

        if(soldier.Fit(dam)>flag.Fit)
            flag.Fit = soldier.Fit(dam);
            flag.xy = soldier.xy(dam,:);
        end
        if(soldier.Fit(Vic)>flag.Fit)
            flag.Fit = soldier.Fit(Vic);
            flag.xy = soldier.xy(Vic,:);
        end
    end
%     disp([num2str(iter),'   ',num2str(flag.Fit)]);
%     flag.xy
    %%%%% mutation for Flag
    [temp_xy] = mutation(flag.xy,flag.Fit,Eval,model);
    temp_Fit = Eval(model,flag.xy);
    if temp_Fit>flag.Fit
        flag.Fit = temp_Fit;
        flag.xy = temp_xy;
    end 
    cg_curve(iter) = flag.Fit;
    %%%%%%%%%
end
ListNFE(R) = NFE;

Res.fit(R) = flag.Fit;
MeanTim = toc/NumRun;
MeanBest = mean(Res.fit);
Best = max(Res.fit);
Worst = min(Res.fit);
StdBest = std(Res.fit);

Res.cg_curve = cg_curve;
Res.MeanTim =  MeanTim;
Res.MeanBest = MeanBest;
Res.Best = Best;
Res.Worst = Worst;
Res.StdBest = StdBest;
% Res.NFE = mean(ListNFE);
Res.Dim = Dim;
Res.sol = flag.xy;
Convergence = cg_curve;
Best_Obj = flag.Fit;
Best_Fit = CalculateFitness(Best_Obj);
Best_Sol = flag.xy;
GAP = 0;
minIter=3000;
notUpdateTimes = 0;
end
