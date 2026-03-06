%Script used to simulate data for the parameter and model recovery analysis
%This script uses functions forom the PALAMEDES toolbox (Kingdom and Prins)
%Author: Arthur S. Courtin
%Licence:MIT
%% Set parameters
n_dataset=100;
n_participant=50;
n_trial=80;
n_dp=n_participant*n_trial;

%% Initialize the thresholding object
addpath('Palamedes')

PM_o=PAL_AMPM_setupPM(...
    'stimRange',-50.5:1:50.5,...
    'priorAlphaRange',-50.5:1:50.5,...
    'priorBetaRange',log10(1./[20:-.1:.1]),...
    'priorLambdaRange',0.02,...
    'gammaEQlambda',true,...
    'PF',@PAL_CumulativeNormal,...
    'numTrials',n_trial...
    );
%% Simulate data for the Gaussian PF
gaussian=table();
for idx=1:n_dataset

    mu_alpha=normrnd(0,10);
    tau_alpha=0;
    while tau_alpha<=0
        tau_alpha=abs(normrnd(0,10));
    end

    mu_beta=normrnd(-1.3,.5);
    tau_beta=0;
    while tau_beta<=0
        tau_beta=abs(normrnd(0,.5));
    end
    
    mu_lambda=normrnd(-4,1);
    tau_lambda=0;
    while tau_lambda<=0
        tau_lambda=abs(normrnd(0,1));
    end

    gaussian.idx((1:n_dp)+(idx-1)*n_dp)=idx;
    gaussian.mu_alpha((1:n_dp)+(idx-1)*n_dp)=mu_alpha;
    gaussian.mu_beta((1:n_dp)+(idx-1)*n_dp)=mu_beta;
    gaussian.mu_lambda((1:n_dp)+(idx-1)*n_dp)=mu_lambda;
    gaussian.tau_alpha((1:n_dp)+(idx-1)*n_dp)=tau_alpha;
    gaussian.tau_beta((1:n_dp)+(idx-1)*n_dp)=tau_beta;
    gaussian.tau_lambda((1:n_dp)+(idx-1)*n_dp)=tau_lambda;

    for p=1:n_participant
        tic
        disp(['Dataset: ',num2str(idx),' - Participant: ',num2str(p)])

        alpha=normrnd(mu_alpha,tau_alpha);
        beta=exp(normrnd(mu_beta,tau_beta));
        lambda=.5/(1+exp(-normrnd(mu_lambda,tau_lambda)));

        PM=PM_o;
        
        for t=1:n_trial
            x=PM.xCurrent;
            y=binornd(1,lambda+(1-2*lambda)*normcdf(x,alpha,1/beta));
            PM=PAL_AMPM_updatePM(PM,y);
        end

        gaussian.participant((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=p;
        gaussian.alpha((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=alpha;
        gaussian.beta((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=beta;
        gaussian.lambda((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=lambda;
        gaussian.alpha_estimate((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=PM.threshold(end);
        gaussian.alpha_se((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=PM.seThreshold(end);
        gaussian.beta_estimate((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=PM.slope(end);
        gaussian.beta_se((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=PM.seSlope(end);
        gaussian.stimuli((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=PM.x(1:end-1)';
        gaussian.response((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=PM.response';

        toc
    end
end
writetable(gaussian,'datasets/Recovery/R_gaussian.csv');

%% Simulate data for the Gaussian PF
gumbel=table();
for idx=1:n_dataset

    mu_alpha=normrnd(0,10);
    tau_alpha=0;
    while tau_alpha<=0
        tau_alpha=abs(normrnd(0,10));
    end

    mu_beta=normrnd(-1.9,.5);
    tau_beta=0;
    while tau_beta<=0
        tau_beta=abs(normrnd(0,.5));
    end
    
    mu_lambda=normrnd(-4,1);
    tau_lambda=0;
    while tau_lambda<=0
        tau_lambda=abs(normrnd(0,1));
    end

    gumbel.idx((1:n_dp)+(idx-1)*n_dp)=idx;
    gumbel.mu_alpha((1:n_dp)+(idx-1)*n_dp)=mu_alpha;
    gumbel.mu_beta((1:n_dp)+(idx-1)*n_dp)=mu_beta;
    gumbel.mu_lambda((1:n_dp)+(idx-1)*n_dp)=mu_lambda;
    gumbel.tau_alpha((1:n_dp)+(idx-1)*n_dp)=tau_alpha;
    gumbel.tau_beta((1:n_dp)+(idx-1)*n_dp)=tau_beta;
    gumbel.tau_lambda((1:n_dp)+(idx-1)*n_dp)=tau_lambda;

    for p=1:n_participant
        tic
        disp(['Dataset: ',num2str(idx),' - Participant: ',num2str(p)])

        alpha=normrnd(mu_alpha,tau_alpha);
        beta=exp(normrnd(mu_beta,tau_beta));
        lambda=.5/(1+exp(-normrnd(mu_lambda,tau_lambda)));

        PM=PM_o;
        
        for t=1:n_trial
            x=PM.xCurrent;
            y=binornd(1,lambda+(1-2*lambda)*(1-exp(-1*10^(beta*(x-alpha)))));
            PM=PAL_AMPM_updatePM(PM,y);
        end
        
        gumbel.participant((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=p;
        gumbel.alpha((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=alpha;
        gumbel.beta((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=beta;
        gumbel.lambda((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=lambda;
        gumbel.alpha_estimate((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=PM.threshold(end);
        gumbel.alpha_se((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=PM.seThreshold(end);
        gumbel.beta_estimate((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=PM.slope(end);
        gumbel.beta_se((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=PM.seSlope(end);
        gumbel.stimuli((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=PM.x(1:end-1)';
        gumbel.response((1:n_trial)+(p-1)*n_trial+(idx-1)*n_dp)=PM.response';

        toc
    end
end
writetable(gumbel,'datasets/Recovery/R_gumbel.csv');
