%This code helps figure out how code runtime will scale with processing
%speed and core count.

clear all;
close all;
load('thread_timer.mat');

%Check data format.
if strcmp(strData(2),'nThreads') && ...
   strcmp(strData(4),'complete-run') && ...
   strcmp(strData(10),'ode-parallel-integration') && ...
   strcmp(strData(11),'ode-propagateLR') && ...
   strcmp(strData(16),'wv-calc')
    %Okay.
else
    error('Where is the right data?');
end

threadVec = sort(unique(tData(:,2)));
threadBestIdx = 0*threadVec;
nThread = length(threadVec);
for k = 1:nThread
    kBestIdx = find(tData(:,2) == threadVec(k));
    kBestIdx = kBestIdx(tData(kBestIdx,4)==min(tData(kBestIdx,4)));
    threadBestIdx(k) = kBestIdx;
end

znT = zeros(nThread,1);
threadTotTime = znT;
threadEquilInputTime = znT;
threadFourFitTime = znT;
threadODEIntegrateTime = znT;
threadODEPropagateTime = znT;
threadVacTime = znT;
for k = 1:nThread
    threadTotTime(k) = tData(threadBestIdx(k),4);
    threadEquilInputTime(k) = tData(threadBestIdx(k),5);
    threadFourFitTime(k) = tData(threadBestIdx(k),7);
    threadODEIntegrateTime(k) = tData(threadBestIdx(k),10);
    threadODEPropagateTime(k) = tData(threadBestIdx(k),11);
    threadVacTime(k) = tData(threadBestIdx(k),16);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Build up the total time--see decompositions work.
%Should really take max(proxyTotTimeNoVac,proxyTotTimeVac)
proxyTotTimeNoVac = threadEquilInputTime   + ...
                    threadFourFitTime      + ...
                    threadODEIntegrateTime + ...
                    threadODEPropagateTime;
P = polyfit(threadTotTime,proxyTotTimeNoVac,1);
extraTimeOffset = P(2);

proxyTotTimeVac =  threadEquilInputTime   + ...
                   threadFourFitTime      + ...
                   threadVacTime;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now get a model for each timer
nModel = 100;
modelThreadVec = (1:nModel)';
modelEquilInputTime = mean(threadEquilInputTime)*ones(nModel,1); %Never really changes
modelFourFitTime = min(threadFourFitTime)*ones(nModel,1); %Bottoms out due to false sharing.
modelODEPropagateTime = mean(threadODEPropagateTime)*ones(nModel,1); %Barely changes with # of intervals.

iPt = 4; %This starting point makes the following two models agree very well on the available data.
f = @(b,x) b(1).*exp(b(2).*x) + b(3);
nrmrsd = @(b) norm(threadODEIntegrateTime(iPt:end) - f(b,threadVec(iPt:end)));
B0 = rand(3,1);
[B,~] = fminsearch(nrmrsd, B0);
ef = fit(threadVec(iPt:end),threadODEIntegrateTime(iPt:end),'exp2');
modelODEIntegrateTime = B(1)*exp(B(2)*modelThreadVec)+B(3);
%figure;
%plot(threadVec,threadODEIntegrateTime,'.','MarkerSize',20);
%hold on;
%plot(threadVec,ef.a*exp(ef.b*threadVec) + ef.c*exp(ef.d*threadVec));
%plot(threadVec,B(1)*exp(B(2)*threadVec)+B(3));

iPt = 4; %This starting point makes the following two models agree very well on the available data.
fVac = @(b,x) b(1)./(exp(b(2).*x)) + b(3);
nrmrsdVac = @(b) norm(threadVacTime(iPt:end) - fVac(b,threadVec(iPt:end)));
B0 = rand(3,1);
[B,~] = fminsearch(nrmrsdVac, B0);
ef = fit(threadVec(iPt:end),threadVacTime(iPt:end),'exp2');
modelVacTime = B(1)./(exp(B(2)*modelThreadVec))+B(3);
%figure;
%plot(threadVec,threadVacTime,'.','MarkerSize',20);
%hold on;
%plot(threadVec,ef.a*exp(ef.b*threadVec) + ef.c*exp(ef.d*threadVec));
%plot(threadVec,B(1)./exp(B(2)*threadVec)+B(3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Now model different processor speeds.
refProcSpeed = 2.7; %Intel(R) Xeon(R) CPU E5-2680 v4 @ 2.40GHz (but runs at an average of 2.7GHz, maxes out at 3.3)
modelProcSpeed = .1:.1:4.4;
modelThreadNum = 1:44;
[P,T] = meshgrid(modelProcSpeed,modelThreadNum);
projectedRunTime = extraTimeOffset          + ...
                   modelEquilInputTime(T)   + ...
                   modelFourFitTime(T)      + ...
                   max(modelVacTime(T), modelODEPropagateTime(T) + modelODEIntegrateTime(T));
projectedRunTime = projectedRunTime * refProcSpeed ./ P;

for iProc = 1:9
    if iProc == 1
        procStr = 'Intel® Xeon® Platinum 8156 Processor';
        nCores = 4;
        nThreads = 8;
        fBase = 3.6;
        fTurbo = 3.7;
    elseif iProc == 2
        procStr = 'Intel® Xeon® Platinum 8158 Processor';
        nCores = 12;
        nThreads = 24;
        fBase = 3;
        fTurbo = 3.7;
    elseif iProc == 3
        procStr = 'Intel® Xeon® Platinum 8160 Processor';
        nCores = 24;
        nThreads = 48;
        fBase = 2.1;
        fTurbo = 3.7;
    elseif iProc == 4
        procStr = 'Intel® Xeon® Platinum 8168 Processor';
        nCores = 24;
        nThreads = 48;
        fBase = 2.7;
        fTurbo = 3.7;
    elseif iProc == 5
        procStr = 'Intel® Xeon® Platinum 8180 Processor';
        nCores = 28;
        nThreads = 56;
        fBase = 2.5;
        fTurbo = 3.8;
    elseif iProc == 6
        procStr = 'INTEL® CORE i7-7820X';
        nCores = 8;
        nThreads = 16;
        fBase = 3.6;
        fTurbo = 4.3;
    elseif iProc == 7
        procStr = 'INTEL® CORE i7-8700K';
        nCores = 6;
        nThreads = 12;
        fBase = 3.7;
        fTurbo = 4.7;
    elseif iProc == 8
        procStr = 'INTEL® CORE i9-7940X';
        nCores = 14;
        nThreads = 28;
        fBase = 3.1;
        fTurbo = 4.3;
    elseif iProc == 9
        procStr = 'INTEL® CORE i7-6900K';
        nCores = 8;
        nThreads = 16;
        fBase = 3.2;
        fTurbo = 3.7;
    end
    
    disp([procStr ': ' num2str(nCores) ' cores [' num2str(fBase) 'GHz] ', num2str(round(projectedRunTime(nCores,10*fBase)*1000)) 'ms']);
end