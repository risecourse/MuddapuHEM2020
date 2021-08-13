%% CREDITS
% Created by
% Vignayanandam R. Muddapu (Ph.D. scholar)
% C/o Prof. V. Srinivasa Chakravarthy
% Indian Institute of Technology Madras
% India

% RUN script for Glutamate-induced toxicity model
% Included:
% Energy deficiency
% Glutamate inhibition therapy
% Dopamine restoration therapy
% Calcium channel blockers therapy
% Calcium-binding protein therapy
% Apoptotic signal blocker therapy

%%
clc;clear;close all;
% plotting parameters
skip_plots = 12;    % the higher the number, the fewer the cells whose
                    % membrane potential traces are plotted. And keep in
                    % mind that we only recorded the ongoing voltages for
                    % a row of cells of each type. If there is any
                    % topological pattern of activity in this model, we are
                    % likely missing it. You may want to look at a random
                    % sampling of cells instead of the top row, if it seems
                    % like they could be different.

% Duration of simulation
dur=50000; % ms

% time step of simulation
dt=0.1; %ms


% Recording simulation time
time=clock;curdate=time(3);curmonth=time(2);

% Percentage of cells under energy deficiency
peren=[100];

% Weight of STN-->SNc
wstsn=[0.3];
scfa=0.00001;% scaling factor

% Apoptotic threshold
Aapopthr=[0.5];

nR=1; %0-autoreceptors 1-without autoreceptors
camtthr=[0.02]; % Perhaps this is "calcium concentration threshold after which mitochondria-induced apoptosis gets initiated"

numtrial=1;
if gpuDeviceCount("available")>0
    gpuon=1;% Code runs on GPU
else
    gpuon=0;%  Code runs on CPU only

end
% Percentage of cell loss at which therapy is initiated
cl=[25];

% Glutamate inhibition therapy
gion=0;
gi_dose=[0]; %all are supposed to be at .3, putting it at 0 so there is no effectiveness

% Dopamine restoration therapy
dron=0;
dr_dose=[0];

% Calcium channel blockers therapy
ccbon=0;
ccb_dose=[0];

% Calcium-binding protein therapy
cbdon=0;
cbd_dose=[0];

% Apoptotic signal blocker therapy
asbon=0;
asb_dose=[0];

scfa1=deci2str(scfa);
durr=deci2str(dur/1000);

for i=1:numel(peren)
    for j=1:numel(wstsn)
        for k=1:numel(Aapopthr)
            for l=1:numel(camtthr)
                for m=1:numel(cl)
                    for n4=1:numel(gi_dose)
                        for n3=1:numel(dr_dose)
                            for n2=1:numel(ccb_dose)
                                for n1=1:numel(cbd_dose)
                                    for n=1:numel(asb_dose)
                                        % create a struct that stores
                                        % parameter values for easy
                                        % documentation inside the saved
                                        % data *.mat files:
                                        params.ApoptosisThresh = Aapopthr(k);
                                        params.Weight_STN_SNc = wstsn(j);
                                        params.Under_energy_def = peren(i);
                                        params.CaThresh_MT_Apop = camtthr(l);
                                        params.CellLossLevel_Therapy = cl(m);
                                        params.ApopBlock_Therapy_Dose = asb_dose(n);
                                        
                                        wwstsn=deci2str(wstsn(j));
                                        peren1=deci2str(peren(i));
                                        apopthr=deci2str(Aapopthr(k));
                                        
                                        camtthr1=deci2str(camtthr(l));
                                        cl1=deci2str(cl(m));
                                        asb_dose1=deci2str(asb_dose(n));

                                        filename=strcat('H_t',num2str(apopthr),'_camtt',num2str(camtthr1),'_PD',num2str(peren1),'_Rstnsnc',num2str(wwstsn),'_std1_scfa',num2str(scfa1),'_CL',num2str(cl1),'%_',num2str(numtrial),'_',num2str(durr),'sec');

                                        
                                        %disp(filename)
                                        [deda,dDA,kid,simtime,srnd,VtrajectorySTN,VtrajectoryGPE,VtrajectorySNC,Ca_trajectorySNC,mCatrajectorySNC,Ttime,Nstn,erCatrajectorySNC,Ca_er,Ca_mt,er_catrajectorySNC,mt_catrajectorySNC,GLCe,GLCn,GLCn_trajectory,GLCe_trajectory]=MAIN_HEM_model(dt,dur,peren(i),wstsn(j),scfa,Aapopthr(k),camtthr(l),cl(m),gion,gi_dose(n4),dron,dr_dose(n3),ccbon,ccb_dose(n2),cbdon,cbd_dose(n1),asbon,asb_dose(n),gpuon);

                                        %                         [out_cell,simtime,srnd]=IS_bSNc_iSTN_GPe_gpu_long_CL_cell(dur,peren(i),wstsn(j),scfa,Aapopthr(k),camtthr(l),cl(m),gi_dose(n),gpuon);
                                        %                         [deda,dDA,kid,clp,simtime,srnd]=IS_bSNc_iSTN_GPe_gpu_long_CL_pattern(dur,peren(i),wstsn(j),scfa,Aapopthr(k),camtthr(l),cl(m),gi_dose(n),gpuon);
                                        parsave_CL(filename,deda,dDA,kid,simtime,srnd,params);
                                        %                         parsave_CL_cell(filename,out_cell,simtime,srnd);
                                        %                         parsave_CL_pattern(filename,out_cell,simtime,srnd);
                                        disp(sprintf(['%%The file results are in: ',filename, '.mat']))
                                        disp(sprintf(['\n%%You can load them in with the command:']))
                                        disp(['results=load(''',filename,''');'])

                                        % Here are some examples of how to
                                        % print the cell voltages:
                                        h=figure();
                                        for q=1:skip_plots:size(VtrajectorySTN,2)
                                            plot(dt:dt:dur, VtrajectorySTN(:,q,1))
                                            hold on
                                        end 
                                        xlabel('Time (ms)')
                                        ylabel('Membrane Potential (mV)')
                                        title({'STN Cells',filename},Interpreter='None')
                                        h2=figure();
                                        for q=1:skip_plots*skip_plots:size(VtrajectoryGPE,2)
                                            plot(dt:dt:dur, VtrajectoryGPE(:,q))
                                            hold on % hold on means, keep the current lines 
                                                    % on the figure and add the newest plot line on top
                                                    % You can also just
                                                    % plot one cell, or
                                                    % introduce subplots or
                                                    % new figures in the
                                                    % for loop so that each
                                                    % cell's v is plotted
                                                    % separately
                                        end 
                                        xlabel('Time (ms)')
                                        ylabel('Membrane Potential (mV)')
                                        title({'GPE Cells',filename},Interpreter='None')
                                        
                                        h3=figure();
                                        for q=1:size(VtrajectorySNC,2)
                                            plot(dt:dt:dur, VtrajectorySNC(:,q))
                                            hold on
                                        end 
                                        xlabel('Time (ms)')
                                        ylabel('Membrane Potential (mV)')
                                        title({'SNC Cells',filename},Interpreter='None')   
                                        
                                         h4=figure();
                                        for q=1:size(Ca_trajectorySNC,2)
                                            plot(dt:dt:dur, Ca_trajectorySNC(:,q))
                                            hold on
                                        end 
                                        xlabel('Time (ms)')
                                        ylabel('Internal Ca2+ concentration (mM)')
                                        title({'[Ca2+]i',filename},Interpreter='None')   
                                        
                                        % We can also calculate an array
                                        % such that the number of elements
                                        % in the array equals the number of
                                        % actual step counts. Then,
                                        % each time the array contains a
                                        % value of 1 we know the cell
                                        % spiked during that spike raster.
                                        
                                        % Here's some code that gets the
                                        % spikes from SNC cell #3. We can
                                        % create a vector with one element
                                        % for each time step, where it is 1
                                        % if the cell reaches spike
                                        % threshold during that time step
                                        % and 0 otherwise. This 'spikes'
                                        % vector we can convolve with a
                                        % "stress trajectory" kernel to
                                        % come up with a continuous
                                        % estimation of stress level in
                                        % this cell throughout the
                                        % simulation. We can also create a
                                        % list of times at which this cell
                                        % spiked 'spktimes' if useful for
                                        % plotting or other analyses:
                                        spikes = getspikes(VtrajectorySNC(:,1),-30);
                                        timevec = dt:dt:dur;
                                        spktimes=timevec(spikes==1);
                                        
                                        % For the convolution, you would
                                        % need something that describes the
                                        % trajectory of stress over time in
                                        % response to a spike in the cell.
                                        % This will be your kernel.
                                        % The exact function you could
                                        % decide amongst yourselves what to
                                        % use. And, if you get the stress
                                        % level from the detailed cell and
                                        % deconvolve it with the spiking
                                        % pattern of that cell, that would
                                        % be one way to get a kernel.
                                        % Here I just used an equation
                                        % similar to the one we use to
                                        % model synaptic conductance over
                                        % time. This example is not
                                        % appropriately scaled for use in
                                        % this equation; you would need to
                                        % look at the units and range of
                                        % the other stress measures in this
                                        % model code to find a scale that
                                        % makes sense if you wanted to use
                                        % this equation:
                                        % 
                                        t = 0:1800; % assuming stress
                                        % response to a single spike
                                        % follows a time course lasting
                                        % around 200 ms, our kernel will be
                                        % around 10x that length
                                        %
                                        stress_max = 5; % this is an
                                        % arbitrary example; I haven't
                                        % looked at the stress measure in
                                        % this code to find its scale. Same
                                        % with the time courses below:
                                        trise = 10; % ms
                                        tfall = 1010; % ms
                                        %
                                        stress_traj_mt = mt_catrajectorySNC;
                                        stress_traj_er = er_catrajectorySNC;
                                        %
                                        % Then we convolve the
                                        % stereotypical stress response
                                        % with the vector of spikes to get
                                        % the stress level over time. See
                                        % this site for help with the
                                        % convolution:
                                        % https://www.mathworks.com/help/matlab/ref/conv.html
                                        stress_mt = conv(VtrajectorySNC(:,1),spktimes);
                                        %stress_er = conv(spktimes,stress_traj_er);
                                        %
                                        z = 2;
                                        alpha = 0.05;
                                        w = linspace(1,50001,500000);
                                        w = w.';
%                                         polytool(w/10,mt_catrajectorySNC(:,1),z,alpha);
                                        
                                        
                                        stressData = table(w,VtrajectorySNC(:,1),mt_catrajectorySNC(:,1));
                                        stressDataLinear_mt = table(mt_catrajectorySNC(:,1),w,VtrajectorySNC(:,1));
                                        stressDataLinear_er = table(er_catrajectorySNC(:,1),w,VtrajectorySNC(:,1));
                                        
                                        %x = [w,Var3];
                                        a = [w mt_catrajectorySNC(:,1) + ones(size(mt_catrajectorySNC(:,1))) * 80];

                                        
%                                         mdl = fitlm(stressDataLinear,'Var1~w+Var3');
                                        %%%quadratic regression%%%
%                                         stepWiseMdl_mt = fitlm([w VtrajectorySNC(:,1)],mt_catrajectorySNC(:,1),'quadratic');
%                                         stepWiseMdl_er = fitlm([w VtrajectorySNC(:,1)],er_catrajectorySNC(:,1),'linear');
                                        Mdl_mt = fitrensemble([w VtrajectorySNC(:,1)],mt_catrajectorySNC(:,1));
                                        Mdl_er = fitrensemble([w VtrajectorySNC(:,1)],er_catrajectorySNC(:,1));
                                        %%%stepwise predicitons%%%
                                        XSNC = [w VtrajectorySNC(:,1)];
                                        XSTN = [w VtrajectorySTN(:,1)];
                                        XGPE = [w VtrajectoryGPE(:,1)];
                                        
                                        pred_mt_SNC = predict(Mdl_mt, XSNC);
                                        pred_mt_STN = predict(Mdl_mt, XSTN);
                                        pred_mt_GPE = predict(Mdl_mt, XGPE);
                                        pred_er_SNC = (predict(Mdl_er, XSNC));
                                        pred_er_STN = predict(Mdl_er, XSTN);
                                        pred_er_GPE = predict(Mdl_er, XGPE);
                                        
                                        
                                        
                                        figure;plot(pred_mt_SNC);
                                        title('Predicted Calcium levels in SNC Mitochondria');
                                        xlabel('Time (ms)');
                                        ylabel('Ca2+ Concentration (mM)');
                                        figure;plot(pred_mt_STN);
                                        title('Predicted Calcium levels in STN Mitochondria');
                                        xlabel('Time (ms)');
                                        ylabel('Ca2+ Concentration (mM)');
                                        figure;plot(pred_mt_GPE);
                                        title('Predicted Calcium levels in GPE Mitochondria');
                                        xlabel('Time (ms)');
                                        ylabel('Ca2+ Concentration (mM)');
                                        figure;plot(pred_er_SNC);
                                        title('Predicted Calcium levels in SNC Endoplasmic Reticulum');
                                        xlabel('Time (ms)');
                                        ylabel('Ca2+ Concentration (mM)');
                                        figure;plot(pred_er_STN);
                                        title('Predicted Calcium levels in STN Endoplasmic Reticulum');
                                        xlabel('Time (ms)');
                                        ylabel('Ca2+ Concentration (mM)');
                                        figure;plot(pred_er_GPE);
                                        title('Predicted Calcium levels in GPE Endoplasmic Reticulum');
                                        xlabel('Time (ms)');
                                        ylabel('Ca2+ Concentration (mM)');
                                        figure;plot(mt_catrajectorySNC);
                                        title('Actual Calcium levels in SNC Mitochondria');
                                        xlabel('Time (ms)');
                                        ylabel('Ca2+ Concentration (mM)');
                                        figure;plot(er_catrajectorySNC);
                                        title('Actual Calcium levels in SNC Endroplasmic Reticulum');
                                        xlabel('Time (ms)');
                                        ylabel('Ca2+ Concentration (mM)');
                                        
                                        
                                            
                                        %%%ridge fit%%%
                                        X1 = [w VtrajectorySNC(:,1)];
                                        d = x2fx(X1,'interaction');
                                        d(:,1) = [];
                                        k = 0:1e-5:5e-3;
                                        b = ridge(er_catrajectorySNC(:,1),d,k);
                                        
%                                         n = length(mt_catrajectorySNC(:,1));
%                                         rng('default');
%                                         a = cvpartition(n,'HoldOut',0.3);
%                                         c_train = training(a,1);
%                                         c_test = ~c_train;
%                                         k_1 = 5;
%                                         b_1 = ridge(mt_catrajectorySNC(:,1)(c_train),X1(c_train,:),k,0);
                                      

                                        


                                        
                                        
                                        x = table(w,VtrajectorySNC(:,1));
                                        y = er_catrajectorySNC(:,1);
                                        Mdl = fitrnet(x,y);


                                        
                                        %beta0 = randn(nVars,1);
                                        
                                        %double sigmoid - 0.5*( M1/(1+exp(-k1(x-a))) + M2/(1+exp(-k2(x-b)))) Where M1 and M2 are the asymptotics (the plateaus), a and b are the x values for the centers of the slopes and k1 and k2 are connected to the steepness of the curves.
                                        modelfun = @(b,x) b(1).*x{:,1}.^b(2) ...
                                             + b(3).* x{:,2}.^b(4) + b(5);
                                        
%                                         modelfun = @(b,x)b(1)/(1+exp(-b(2)*(b(3)-x(1,:))))
%                                         modelfun = @(b,x)(-.5)*(b(1)/(1+exp(-b(2)*(x(1,:)-b(3))))+...
%                                             b(4)/(1+exp(-b(5)*(x(2,:)-b(6)))));
                                        MaxFunEvals = 50;
                                        beta0 = [0 0 0 0 0];  
%                                         mdl1 = fitnlm(stressData,modelfun,beta0);
                                        SSECF = @(b) sum((i(:) - modelfun(b,x)).^2);
                                        [b_estd, SSE1] = fminsearch(SSECF, beta0);
                                        i_fit1 = modelfun(b_estd,stressData);
                                        
%                                         
%   
%  
%                                         
% %                                         load reaction;
%                                         ds = dataset({independent,xn(2,:),xn(3,:)},{dependent,yn});
%                                         i_fcn = @(b,x) b(1).*x(:,1).^b(2) ...
%                                             + b(3).* x(:,2).^b(4) + b(5);
%                                         x = [w(:)  VtrajectorySNC(:,1)];
%                                         b0 = rand(4,1);
%                                         SSECF = @(b) sum((i(:) - i_fcn(b,x)).^2);                           % Sum-Squared-Error Cost Function
%                                         [b_estd, SSE] = fminsearch(SSECF, b0);                              % Estimate Parameters
%                                         i_fit = i_fcn(b_estd,x);  
                                        




                                        
                                        % And we can shift it so that the
                                        % effects take place after the
                                        % spike time:
                                        %stress_mt = [zeros(1,round(length(stress_traj_mt)/2)) stress_mt(1:end-round(length(stress_traj_mt)/2))];
                                        %
                                        % OR take the full convolution
                                        % instead of the same size and only
                                        % use the first N indices of the
                                        % stress vector, where N is the
                                        % length of timevec. Ex:
                                        % figure;plot(timevec,stress(1:length(timevec)))
                                        %
                                        % And to plot with spikes overlaid:
                                        % for x =1:length(spikes)
                                        %     if spikes(x)==1
                                        %         hold on
                                        %         plot([timevec(x) timevec(x)],ylim, 'r-')
                                        %     end
                                        % end
                                        
                                        %load train
                                        %sound(y,Fs)
                                   end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end


% close all
% clear;clc;
% matlabpool close
toc
threshold_vec_m = zeros(Ttime, 8);
threshold_vec_er = zeros(Ttime, 8);
threshold_vec_in = zeros(Ttime, 8);
% for i = 1:5
%     threshold_vec(CatrajectorySNC(k, i)>0.00215)=1;
% end

for k=1:Ttime
    for t=1:8
        if mt_catrajectorySNC(k,t) > .0215
            threshold_vec_m(k,t) = 1;
        end
        if er_catrajectorySNC(k,t) > .00215
            threshold_vec_er(k,t) = 1;
        end
    end
end
%threshold_vec_m(mt_catrajectorySNC(:) > 0.0215)=1;
%threshold_vec_er(er_catrajectorySNC(:) > 0.00215)=1;
threshold_vec_in(Ca_trajectorySNC(:) > 0.02)=1;


function spikes = getspikes(v,threshold)
    spikes = zeros(size(v));
    s = find(v(1:end-1)<=threshold & v(2:end)>threshold);
    spikes(s) = 1;
end



