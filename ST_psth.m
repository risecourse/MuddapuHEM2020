function [temp1]=ST_psth(matrix)

% Extracting spike times from neuronal firing data
% matrix = firings (column1->Time point; column2->Neuron ID)
% temp1 = spike times for each neuron (row->neuron ID; column->timepoints)

neur=unique(matrix(:,2));
temp1=[];
for i=1:numel(neur)
temp=matrix((matrix(:,2)==neur(i)));
temp1=[temp1;temp];
% temp=temp0(temp0>start & temp0<stop);
% temp1{i}=temp;
end
% sr=[];
% for j=1:numel(temp1)
%     sr1=psth(temp1(:,:),500,10000,1,100001,0);
%     sr=[sr sr1];
% %     pause(0.01)
% end

% figure(2)
% imagesc(sr')
% colorbar

end