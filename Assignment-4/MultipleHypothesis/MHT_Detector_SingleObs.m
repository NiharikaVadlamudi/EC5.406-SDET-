%% Assignment 2 Multiple Hypothesis Testing. 
A=1;
f=1;

% Testing signal , close to signal k=2 + added noise . 
sm=[cos(2*(2*4-1)*pi/8) sin(2*(2*4-1)*pi/8)];
L=size(sm);
w=normrnd(0,0.05,L);
x_test=sm+w;

%Result 
res=detector(x_test);
fprintf("Input Signal(x) matches Hypothes H_%d ; i.e S_%d\n",res,res)




%% Function Definations
function index = detector(x)
    A=1;
    s1=signal(1,1);
    s2=signal(1,2);
    s3=signal(1,3);
    s4=signal(1,4);
    s5=signal(1,5);
    s6=signal(1,6);
    s7=signal(1,7);
    s8=signal(1,8);

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

function s=signal(A,k)
    s=[A*cos(2*(2*k-1)*pi/8) A*sin(2*(2*k-1)*pi/8)];
end 

function w = noise(sigma_val)
    w=[normrnd(0,sigma_val) normrnd(0,sigma_val)];
end

function test=test_statistics(x,s)
    test=mtimes(x,transpose(s));
end 
