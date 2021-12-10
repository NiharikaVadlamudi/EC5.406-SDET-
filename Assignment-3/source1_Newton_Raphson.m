%% Script for Source 1 Localisation 

%%Global Variables.
p = 1 ; 
SNR=80;
sigma=1/db2pow(SNR);
beta=20/sigma*(log(10));
beta_sq=(20./sqrt(sigma)*log(10));
eps=0.005;
% Grid size of 10 x 10 .
radius=50;
% Fixing source 1 
s1_pos_gd=[4;1];

%Fixed Entity for a given experiment 
% snr_db=1:5:50;
% sigma_list=1./db2pow(snr_db);
% L=length(sigma_list);
% sigma=sigma_list(i);
% beta=20/sigma*(log(10));
% beta_sq=(20./sqrt(sigma)*log(10));

M=10;
meanErrorList=[];
N_list=[40,60,80,100,120];
Rpos=generateRandomPosition(M,radius);
L=length(N_list);
for i=1:1:L
    N=N_list(i);
    thetha_sum=[0;0];
    thetha_init=[0;0];
    for n=1:1:N
        yR_n=recieverPower(M,sigma,radius);
        thetha_opt=newtonRaphson(M,Rpos,thetha_init,eps,yR_n,beta,beta_sq,p,radius);
        thetha_sum=thetha_sum+thetha_opt;
    end
    thetha_final=(1/N)*(thetha_sum);
    meansqerr=norm(thetha_final-s1_pos_gd);
    meanErrorList(i)=meansqerr;
    fprintf('Iter : %d , MSE : %f Thetha ->[%d %d]\n:',i,meansqerr,int16(thetha_final))
    fprintf('\n')
end 

%% Plotting- Computing for multiple values of SNR figure
f1 = figure();
set(f1, 'Visible', 'off');
title('mse vs snr ')
plot(N_list,meanErrorList)
xlabel('N')
ylabel('MSE')
saveas(gcf,'res1.png')

%% Helper Functions 

% Source 1 
function f1=s1(p)
    f1=p;
end

% g(x,r) function 
function res=g(r,l)
    rx=r(1);
    ry=r(2);
    lx=l(1);
    ly=l(2);
    res=(rx-lx).^2+(ry-ly).^2;
end 

% f function 
function fval = ft(p,r,l)
    fval=10*log10(s1(p)./g(r,l)); 
end 

% Random Position Generator - Rx nodes.
function R=generateRandomPosition(M,radius) 
    R_x=randi([0 radius],1,M)+ randn([1 M]);
    R_y=randi([0 radius],1,M)+ randn([1 M]);
    R=[R_x;R_y];
end 

% Random Signal Value generator for Rx 
function yR=recieverPower(M,sigma,radius)
    w=sigma*randn([1 M]);
    dterm=radius.^2/2*randn([1 M])+0.25*radius.^2;
    p=1;
    for i=1:1:M
        yR(i)=10*log10(p*1/dterm(i))+w(i)*(sigma)*rand(1);
    end
    yR=transpose(yR);
end 

% d(lnP)/dx function 
function val = dlnpx(r,y,l,factor,p)
    rx=r(1);
    lx=l(1);
    num= (y-ft(p,r,l))*(rx-lx);
    den=g(r,l);
    val=factor*(num./den);
end 

% d(lnP)/dy function for s1(n) signal . 
function val = dlnpy(r,y,l,factor,p)
    ry=r(2);
    ly=l(2);
    num=(y-ft(p,r,l))*(ry-ly);
    den=g(r,l);
    val=factor*(num./den);
end 

% @S1First derivative of log-likelihood function [ d(ln(p)/dx d(ln(p)/dy ]
function p1=lnpvec(M,R,Y,l,factor,p)
    fx=0;
    fy=0;
    for i=1:1:M
        r_i=R(:,i);
        y_i=Y(i);
        fx=fx+dlnpx(r_i,y_i,l,factor,p);
        fy=fy+dlnpy(r_i,y_i,l,factor,p);
    end 
    p1=[fx;fy];
end 

% Single Fisher Matrix . 
function I_m=fisher_m(r,l,factor)
 rx=r(1);
 ry=r(2);

 lx=l(1);
 ly=l(2);
 
 g_sq= g(r,l).^2;

 % Matrix Components 
 I_a=(rx-lx).^2/g_sq;
 I_b=(rx-lx)*(ry-ly)/g_sq;
 I_c=I_b;
 I_d=(ry-ly).^2/g_sq;
 
 I_m=factor*[[I_a I_b];[I_c I_d]];   
end 

% Complete Fisher Matrix - I for M observations .
function f=I(M,R,l,factor)
    f=[[0 0];[0 0]];
    for i=1:1:M
        r_i=R(:,i);
        f=f+fisher_m(r_i,l,factor);
    end 
    % Perform Inverse of the matrix
%     f=(1/M)*f;
    f=inv(f);
end 

% Newton Raphson Algorithm.
function thetha_opt=newtonRaphson(M,R,thetha_init,eps,Y,beta,beta_sq,p,radius)
    thetha_old=thetha_init;
    factor1=beta;
    factor2=beta_sq;
    err=inf;
    count=0;
    while(err>=eps)
        Q=mtimes(I(M,R,thetha_old,factor1),lnpvec(M,R,Y,thetha_old,factor2,p));
        if(isinf(Q)|isnan(Q))
            thetha_new=thetha_old;
            err=0;
            fprintf('Breaking due to Q-inf\n')
            break
        end
        % Update the vector. 
        thetha_new=thetha_old+Q;
        err=norm(thetha_new-thetha_old);
        thetha_old=thetha_new;
        count=count+1;
    end 
    % Now , the while loop has broken .
    thetha_opt=thetha_new;
    fprintf('Iterations Taken : %d\n',count)
end 





