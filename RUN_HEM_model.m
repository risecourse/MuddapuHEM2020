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

% Duration of simulation
dur=2000;

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
camtthr=[0.02];

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
gi_dose=[0.3];

% Dopamine restoration therapy
dron=0;
dr_dose=[0.3];

% Calcium channel blockers therapy
ccbon=0;
ccb_dose=[0.3];

% Calcium-binding protein therapy
cbdon=0;
cbd_dose=[0.3];

% Apoptotic signal blocker therapy
asbon=0;
asb_dose=[0.3];

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
                                        
                                        wwstsn=deci2str(wstsn(j));
                                        peren1=deci2str(peren(i));
                                        apopthr=deci2str(Aapopthr(k));
                                        camtthr1=deci2str(camtthr(l));
                                        cl1=deci2str(cl(m));
                                        asb_dose1=deci2str(asb_dose(n));
                                        %                         filename=strcat('HybMod_SSG_GDA_stnlatRS_DAr0-5--15e-6__nRRP10_sncstn1_Rstnsnc5-5_',num2str(durr),'msec_',num2str(curdate),'_',num2str(curmonth));
                                        %                         filename=strcat('H_MAp-1_V0_0-1_t',num2str(apopthr),'_CA1_PD',num2str(peren1),'_allr_rA0-1_lamR10e-6-5e-6_Fexp4-87_DA1_2e-5_Rstnsnc',num2str(wstsn),'_std1_scfa',num2str(scfa1),'_1_',num2str(durr),'sec_',num2str(curdate),'_',num2str(curmonth));
                                        filename=strcat('H_t',num2str(apopthr),'_camtt',num2str(camtthr1),'_PD',num2str(peren1),'_Rstnsnc',num2str(wwstsn),'_std1_scfa',num2str(scfa1),'_CL',num2str(cl1),'%_',num2str(numtrial),'_',num2str(durr),'sec');
                                        %                         filename=strcat('H_t',num2str(apopthr),'_camtt',num2str(camtthr1),'_PD',num2str(peren1),'_Rstnsnc',num2str(wwstsn),'_std1_scfa',num2str(scfa1),'_1_',num2str(durr),'sec');
                                        
                                        %disp(filename)
                                        [deda,dDA,kid,simtime,srnd]=MAIN_HEM_model(dur,peren(i),wstsn(j),scfa,Aapopthr(k),camtthr(l),cl(m),gion,gi_dose(n4),dron,dr_dose(n3),ccbon,ccb_dose(n2),cbdon,cbd_dose(n1),asbon,asb_dose(n),gpuon);
                                        %                         [out_cell,simtime,srnd]=IS_bSNc_iSTN_GPe_gpu_long_CL_cell(dur,peren(i),wstsn(j),scfa,Aapopthr(k),camtthr(l),cl(m),gi_dose(n),gpuon);
                                        %                         [deda,dDA,kid,clp,simtime,srnd]=IS_bSNc_iSTN_GPe_gpu_long_CL_pattern(dur,peren(i),wstsn(j),scfa,Aapopthr(k),camtthr(l),cl(m),gi_dose(n),gpuon);
                                        parsave_CL(filename,deda,dDA,kid,simtime,srnd);
                                        %                         parsave_CL_cell(filename,out_cell,simtime,srnd);
                                        %                         parsave_CL_pattern(filename,out_cell,simtime,srnd);
                                        disp(sprintf(['%%The file results are in: ',filename, '.mat']))
                                        disp(sprintf(['\n%%You can load them in with the command:']))
                                        disp(['results=load(''',filename,''');'])
                                        disp(sprintf(['\n%%You could plot them with the command:\nfigure;plot(results.dDA,results.deda)']))
                                        disp(sprintf('\n%%I don''t know what these variables mean yet, having not read the paper...'))
                                        %                         plot_save_long_CL;
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
