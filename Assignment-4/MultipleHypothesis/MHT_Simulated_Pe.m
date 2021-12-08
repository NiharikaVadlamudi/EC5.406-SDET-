A=1;
f=1;
M=8;
SNR_db=(-20:2:20);
SNR_list=vpa(db2pow(SNR_db));
sigma_list=reciprocal(SNR_list);
k_list=1:1:M;

% Iterating through sigma values .
xaxis=[];
yaxis=[];
i=1;
N=1000;
for sig=SNR_db
    sig_val=1./db2pow(sig);
    xaxis(i)=sig;
    yaxis(i)=peCalculation(N,sig_val);
    i=i+1;
end 

%% Plotting 
figure
semilogy(xaxis,yaxis)
% plot(xaxis,yaxis)
title('Pe vs SNR ')
xlabel('ENR - 10log(1/sigma^2))')
ylabel('Pe')

%% Function Definations
function pe = peCalculation(N,sigma_val)
    ntimes=1:1:N;
    countError=0;
    countTrue=0;
    for n=ntimes 
        [gd_n,pred_n]=randGeneratorSignal(1,sigma_val);
        if(gd_n~=pred_n)
            countError=countError+1;
        else
            countTrue=countTrue+1;
        end
    end 
    accuracy=countTrue/N;
    pe=1-accuracy;
    fprintf('Error is : %f \n',pe)
   
end

function [gd,pred]=randGeneratorSignal(A,sigma_val)
    M=8;
    r = randi([1,M],1);
    s_r=[A*cos(2*(2*r-1)*pi/M) A*sin(2*(2*r-1)*pi/M)];
    L=size(s_r);
    w=normrnd(0,1,L);
    x_r=s_r+w;    
    % Assignment of values 
    gd=r;
    pred=detector(x_r,sigma_val);
end
function index = detector(x,sigma_val)
    A=1;
    s1=h_s(1,1,sigma_val);
    s2=h_s(1,2,sigma_val);
    s3=h_s(1,3,sigma_val);
    s4=h_s(1,4,sigma_val);
    s5=h_s(1,5,sigma_val);
    s6=h_s(1,6,sigma_val);
    s7=h_s(1,7,sigma_val);
    s8=h_s(1,8,sigma_val);
    
    % Check if the si signal's size is equal to x .
    if(length(s1)~=length(x))
        return;
    end
    s=[s1;s2;s3;s4;s5;s6;s7;s8];
    L=length(s);
    scores=[];
    for i=1:1:8
        s_i=s(i,:);
        scores(i)=vpa(test_statistics(x,s_i));
    end
    % Finding maximum case . 
    [score,ind]=max(scores);
    index=ind;
end 
function h_s=h_s(A,k,sigma_val)
    h_s=signal(A,k)+noise(sigma_val);
end

function s=signal(A,k)
    s=[A*cos(2*(2*k-1)*pi/8) A*sin(2*(2*k-1)*pi/8)];
end 
function w = noise(sigma_val)
    w=[normrnd(0,sigma_val) normrnd(0,sigma_val)];
end

function test=test_statistics(x,s)
    A=1;
    test=mtimes(x,transpose(s))-0.5^(A*A);
end 
function rec=reciprocal(t)
    rec=vpa(1./t);
end



















% %% Assignment 2 Multiple Hypothesis Testing.
% %% Note : Single Observation 
% 
% A=1;
% f=1;
% M=8;
% N=10;
% SNR_db=vpa(-20:2:20);
% SNR_list=vpa(db2pow(SNR_db));
% sigma_list=reciprocal(SNR_list);
% k_list=1:1:M;
% 
% for sig=sigma_list
%     disp(sig)
%     peCalculation(1,sig)
% end 
% 
% %% Function Definations
% function pe = peCalculation(N,sigma_val)
%     ntimes=1:1:N;
%     countError=0;
%     countTrue=0;
%     for n=ntimes 
%         [gd_n,pred_n]=randGeneratorSignal(1,sigma_val);
%         if(gd_n~=pred_n)
%             countError=countError+1;
%         else
%             countTrue=countTrue+1;
%         end
%     end 
%     accuracy=countTrue/N;
%     fprintf('Error is : \n')
%     error=1-accuracy
% end
% 
% function [gd,pred]=randGeneratorSignal(A,sigma_val)
%     M=8;
%     r = randi([1,M],1);
%     s_r=[A*cos(2*(2*r-1)*pi/8) A*sin(2*(2*r-1)*pi/8)];
%     L=size(s_r);
%     w=normrnd(0,r*0.005,L);
%     x_r=s_r+w;
%     
%     % Assignment of values 
%     gd=r;
%     pred=detector(x_r,sigma_val);
% end
% 
% function index = detector(x,sigma_val)
%     A=1;
%     s1=h_s(1,1,sigma_val);
%     s2=h_s(1,2,sigma_val);
%     s3=h_s(1,3,sigma_val);
%     s4=h_s(1,4,sigma_val);
%     s5=h_s(1,5,sigma_val);
%     s6=h_s(1,6,sigma_val);
%     s7=h_s(1,7,sigma_val);
%     s8=h_s(1,8,sigma_val);
%     
%     % Check if the si signal's size is equal to x .
%     if(length(s1)~=length(x))
%         return;
%     end
%     s=[s1;s2;s3;s4;s5;s6;s7;s8];
%     L=length(s);
%     scores=[];
%     for i=1:1:8
%         s_i=s(i,:);
%         scores(i)=vpa(test_statistics(x,s_i));
%     end
%     % Finding maximum case . 
%     [score,ind]=max(scores);
%     index=ind;
% end 
% 
% function h_i=hypothesis(A,k,sigma_val)
%     h_i=signal(A,k)+noise(sigma_val);
% end
% 
% function s=signal(A,k)
%     s=[A*cos(2*(2*k-1)*pi/8);A*sin(2*(2*k-1)*pi/8)];
% end 
% function w = noise(sigma_val)
%     w=[normrnd(0,sigma_val) normrnd(0,sigma_val)];
% end
% 
% function test=test_statistics(x,s)
%     test=mtimes(x,transpose(s));
% end 
% function rec=reciprocal(t)
%     rec=vpa(1./t);
% end