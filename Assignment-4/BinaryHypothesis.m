N=10000;
A=1;
f = 1/sqrt(2);
SNR_db_list=-20:5:20;
PFE_list=0:0.1:1;

% Signal formation
time=1:1:N;
sig=A*cos(2*pi*f*time);
es=(1/N)*mtimes(sig,transpose(sig));

%% Plotting- Computing for multiple values of SNR figure
% plot(SNR_list,Pe_list)
title('Pd vs Pfa ')
hold on
legendList=[];
i=1;
for snr_val=SNR_db_list
    Pd=detectionProbability(es,snr_val,PFE_list);
    plot(PFE_list,Pd)
    i=i+1;
end 
disp(legendList)
hold off
xlabel('Pfa')
ylabel('Pd')
xlim([0 1])
ylim([0 1])
%lgn=legend(legendList);
lgn=legend('SNR= -20dB','SNR= -15dB','SNR= -10dB','SNR= -5dB','SNR=0dB','SNR =5dB','SNR =10dB','SNR=15dB','SNR=20dB');
lgn.Location='southeast';

%% Functions Defination 
function s=signal(A,f,n)
    s=vpa(A*cos(2*pi*f*n));
end 
function w = noise(sigma_val)
    w=normrnd(0,sigma_val);
end
%% Pd and Pfe calculation 
function Pd=detectionProbability(energy,snr_db,PFE_list)
    sigma=vpa(1/db2pow(snr_db));
    term = sqrt(energy*sigma);
    Pd=[];
    i=1;
    for pfe=PFE_list
        gamma=qfuncinv(term*pfe);
        Pd(i)=qfunc(gamma-energy/term);
        i=i+1;
    end
end