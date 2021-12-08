%% Assignment 2 Multiple Hypothesis Testing.
%% Note : Single Observation 
A=1;
f=1;
M=8;

SNR_db=(-20:2:5);
% Iterating through sigma values .
xaxis=[];
yaxis=[];
i=1;
N=1000;
for sig=SNR_db
    k=1;
    sig_val=1./db2pow(sig);
    xaxis(i)=sig;
    yaxis(i)=1-vpa(pcIntegral(sig_val,1));
    i=i+1;
end 

%% Plotting 
figure
semilogy(xaxis,yaxis)
% plot(xaxis,yaxis)
title('Theoretical Pe vs SNR ')
xlabel('ENR - 10log(1/sigma^2))')
ylabel('Theoretical Pe')

%% Function Definations
function pc=pcIntegral(sigma_val,k)
    hi=hypothesis(k,sigma_val);
    [xl0,xh0,yl0,yh0]=bounds(k);
    pc=abs(vpa(integral2(hi,xl0,xh0,yl0,yh0)));
end 

function [xl,xh,yl,yh]=bounds(k)
    A=1;
    thetha=2*(2*k-1)*pi/8;
    thetha_l=thetha-pi/8;
    thetha_h=thetha+pi/8;
    % Cartesian Coordinate Bounds 
    xl=A*cos(thetha_l);
    xh=A*cos(thetha_h);
    yl=A*sin(thetha_l);
    yh=A*sin(thetha_h);


end 

function hi=hypothesis(k,sigma_val)
    A=1;
    s=signal(A,k);
    hi=@(x,y) normpdf(x,s(1),sigma_val).*normpdf(y,s(2),sigma_val);
end

function s=signal(A,k)
    s=[A*cos(2*(2*k-1)*pi/8);A*sin(2*(2*k-1)*pi/8)];
end 
function w = noise(sigma_val)
    w=[normrnd(0,sigma_val) normrnd(0,sigma_val)];
end

function test=test_statistics(x,s)
    test=mtimes(x,transpose(s));
end 
function rec=reciprocal(t)
    rec=1./t;
end