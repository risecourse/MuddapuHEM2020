function [Brate] = BurstMeasure_firings(firings,Nneur)
%% CREDITS
% Created by
% Vignayanandam R. Muddapu (Ph.D. scholar)
% C/o Prof. V. Srinivasa Chakravarthy
% Indian Institute of Technology Madras
% India

% van Elburg, R. A. J., & van Ooyen, A. (2004). A new measure for bursting.
% Neurocomputing, 58–60, 497–502.
% https://doi.org/10.1016/j.neucom.2004.01.086

%%
Brate=[];dis1=[];dis2=[];
for neur=1:Nneur
    temptime=firings((firings(:,2)==neur));
    
    for i=1:numel(temptime)-2
        dis1(i) = temptime(i+1)-temptime(i);
        dis2(i) = temptime(i+2)-temptime(i);
    end
    Brate(neur)=(2*var(dis1)-var(dis2))/(2*((sum(dis1)/numel(dis1))^2));
%     disp(neur)
end
end