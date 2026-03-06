%aggregate_datasets
clear all
addpath('Palamedes')
es_values=[0 .2 .5 .8 1.2] ;
for tp=1:2
    for es_idx=5
        es=es_values(es_idx);
        dataset=[];
        for idx=1:100
            mu_alpha_intercept=normrnd(-7.89,0.699);
            mu_beta_intercept=normrnd(-2.28,0.0334);
            mu_lambda_intercept=normrnd(-4.74,0.506);
            tau_b_alpha_intercept=0;
            while tau_b_alpha_intercept<=0
                tau_b_alpha_intercept=normrnd(7.89,0.646);
            end
            tau_b_beta_intercept=0;
            while tau_b_beta_intercept<=0
                tau_b_beta_intercept=normrnd(0.145,0.0587); 
            end
            tau_b_lambda_intercept=0;
            while tau_b_lambda_intercept<=0
                tau_b_lambda_intercept=normrnd(0.448,0.340);
            end
            tau_w_alpha_intercept=0;
            while tau_w_alpha_intercept<=0
                tau_w_alpha_intercept=normrnd(7.32,0.429); 
            end
            tau_w_beta_intercept=0;
            while tau_w_beta_intercept<=0
                tau_w_beta_intercept=normrnd(0.333,0.0446);
            end
            tau_w_lambda_intercept=0;
            while tau_w_lambda_intercept<=0
                tau_w_lambda_intercept=normrnd(1.92,0.325);
            end
            tau_b_alpha_effect=0;
            while tau_b_alpha_effect<=0
                tau_b_alpha_effect=normrnd(7.89,0.646)/2; 
            end
            tau_b_beta_effect=0;
            while tau_b_beta_effect<=0
                tau_b_beta_effect=normrnd(0.145,0.0587)/2;
            end
            
            if tp==1
                mu_alpha_effect=es*sqrt(2*tau_b_alpha_intercept^2+2*tau_w_alpha_intercept^2+tau_b_alpha_effect^2);
                mu_beta_effect=0;            
            else
                mu_alpha_effect=0;
                mu_beta_effect=es*sqrt(2*tau_b_beta_intercept^2+2*tau_w_beta_intercept^2+tau_b_beta_effect^2);  
            end
            dataset.target_parameter=tp;
            dataset.effect_size=es;
            dataset.idx=idx;
            dataset.mu_alpha_int=mu_alpha_intercept;
            dataset.mu_beta_int=mu_beta_intercept;
            dataset.mu_lambda_int=mu_lambda_intercept;
            dataset.tau_b_alpha_int=tau_b_alpha_intercept;
            dataset.tau_b_beta_int=tau_b_beta_intercept;
            dataset.tau_b_lambda_int=tau_b_lambda_intercept;
            dataset.tau_w_alpha_int=tau_w_alpha_intercept;
            dataset.tau_w_beta_int=tau_w_beta_intercept;
            dataset.tau_w_lambda_int=tau_w_lambda_intercept;
            dataset.mu_alpha_effect=mu_alpha_effect;
            dataset.mu_beta_effect=mu_beta_effect;
            dataset.tau_b_alpha_effect=tau_b_alpha_effect;
            dataset.tau_b_beta_effect=tau_b_beta_effect;
            
            for p=1:120
                tic
                disp(['Target: ',num2str(tp),' - Effect Size: ',num2str(es),' - Dataset: ',num2str(idx),' - Participant: ',num2str(p)])
            
                participant_alpha_intercept=normrnd(mu_alpha_intercept,tau_b_alpha_intercept);
                participant_beta_intercept=normrnd(mu_beta_intercept,tau_b_beta_intercept);
                participant_lambda_intercept=normrnd(mu_lambda_intercept,tau_b_lambda_intercept);
                
                true_alpha_control=normrnd(participant_alpha_intercept,tau_w_alpha_intercept);
                true_beta_control=normrnd(participant_beta_intercept,tau_w_beta_intercept);
                true_lambda_control=normrnd(participant_lambda_intercept,tau_w_lambda_intercept);
                
                true_alpha_treatment=normrnd(participant_alpha_intercept,tau_w_alpha_intercept)+normrnd(mu_alpha_effect,tau_b_alpha_effect);
                true_beta_treatment=normrnd(participant_beta_intercept,tau_w_beta_intercept)+normrnd(mu_beta_effect,tau_b_beta_effect);
                true_lambda_treatment=normrnd(participant_lambda_intercept,tau_w_lambda_intercept);
                 
                ac=true_alpha_control;
                bc=exp(true_beta_control);
                lc=.5/(1+exp(-true_lambda_control));
                at=true_alpha_treatment;
                bt=exp(true_beta_treatment);
                lt=.5/(1+exp(-true_lambda_treatment));
                
                PM_control=PAL_AMPM_setupPM(...
                    'stimRange',-50.5:1:50.5,...
                    'priorAlphaRange',-50.5:1:50.5,...
                    'priorBetaRange',log10(1./[20:-.1:.1]),...
                    'priorLambdaRange',0.02,...
                    'gammaEQlambda',true,...
                    'PF',@PAL_CumulativeNormal,...
                    'numTrials',150 ...
                    );
                PM_treatment=PM_control;
                
                for t=1:150
                    x=PM_control.xCurrent;
                    y=binornd(1,lc+(1-2*lc)*normcdf(x,ac,1/bc));
                    PM_control=PAL_AMPM_updatePM(PM_control,y);
                    
                    x=PM_treatment.xCurrent;
                    y=binornd(1,lt+(1-2*lt)*normcdf(x,at,1/bt));
                    PM_treatment=PAL_AMPM_updatePM(PM_treatment,y);
                end
                
                dataset.participant(p).alpha_control=ac;
                dataset.participant(p).beta_control=bc;
                dataset.participant(p).lambda_control=lc;
                dataset.participant(p).alpha_treatment=at;
                dataset.participant(p).beta_treatment=bt;
                dataset.participant(p).lambda_treatment=lt;
                dataset.participant(p).alpha_estimate_control=PM_control.threshold;
                dataset.participant(p).alpha_se_control=PM_control.seThreshold;
                dataset.participant(p).beta_estimate_control=PM_control.slope;
                dataset.participant(p).beta_se_control=PM_control.seSlope;
                dataset.participant(p).stimuli_control=PM_control.x(1:end-1)';
                dataset.participant(p).response_control=PM_control.response';
                dataset.participant(p).alpha_estimate_treatment=PM_treatment.threshold;
                dataset.participant(p).alpha_se_treatment=PM_treatment.seThreshold;
                dataset.participant(p).beta_estimate_treatment=PM_treatment.slope;
                dataset.participant(p).beta_se_treatment=PM_treatment.seSlope;
                dataset.participant(p).stimuli_treatment=PM_treatment.x(1:end-1)';
                dataset.participant(p).response_treatment=PM_treatment.response';
                
                toc
            end
            if tp==1         
                save(sprintf('datasets/PI/PI_Threshold_%1.1f_120_150_%i.mat',es,idx),'dataset')
            else
                save(sprintf('datasets/PI/PI_Slope_%1.1f_120_150_%i.mat',es,idx),'dataset')
            end
        end
    end
end

