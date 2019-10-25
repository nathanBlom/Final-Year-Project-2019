%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear all;close all;
%rng('shuffle') % For reproducibility
N = 10^2;
Pin_avg = [-5:5:40];
Pin_avg_dB = 10.^(Pin_avg/10);
repeats = 10^5;
U = zeros(repeats,N);
for j = 1:length(Pin_avg_dB)
    for  k = 1:repeats

        Bmax = 200*Pin_avg_dB(j);
        B=zeros(1,N);
        Pout = zeros(1,N);


        gauss1 = normrnd(0,1,[1,N]);
        gauss2 = normrnd(0,1,[1,N]);
        gamma = gauss1.^2+gauss2.^2;

        Pin = exprnd(Pin_avg_dB(j),[1,N]);

        for i =1:N

            if (i == 1)
                Pout(i) = min([0,Pin_avg_dB(j)]); % Battery is initially emtpy hence 0
            else
                Pout(i) = min([B(i-1),Pin_avg_dB(j)]);
            end

            U(k,i) = log2(1+(Pout(i)*gamma(i)));

            if (i == 1)
                B(i) = 0+Pin(i)-Pout(i);  % Battery is initially emtpy hence 0
            else
                B(i) = B(i-1)+Pin(i)-Pout(i);
            end

            if(B(i)>= Bmax)
                B(i) = Bmax;     % if condition is never met. 
            end
        end
    end
    Ro = 1;
    U_logical = U<Ro;
    U_bar = mean(U_logical(:,2:end),2);
    final_average(j) = mean(U_bar);
    %final_results = sum(sum(U_logical,1))/(repeats* N)
end
semilogy(Pin_avg,final_average,'o-')
grid on
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc;clear all;
%rng('shuffle') % For reproducibility
N = 10^3;
Pin_avg = [-5:5:40];
Pin_avg_dB = 10.^(Pin_avg/10);
repeats = 10^4;
U = zeros(repeats,N);
for j = 1:length(Pin_avg_dB)
    for  k = 1:repeats

        Bmax = 200*Pin_avg_dB(j);
        B=zeros(1,N);
        Pout = zeros(1,N);


        gauss1 = normrnd(0,1,[1,N]);
        gauss2 = normrnd(0,1,[1,N]);
        gamma = gauss1.^2+gauss2.^2;

        Pin = exprnd(Pin_avg_dB(j),[1,N]);

        for i =1:N

            if (i == 1)
                Pout(i) = min([0,Pin_avg_dB(j)]); % Battery is initially emtpy hence 0
            else
                Pout(i) = min([B(i-1),Pin_avg_dB(j)]);
            end

            U(k,i) = log2(1+(Pout(i)*gamma(i)));

            if (i == 1)
                B(i) = 0+Pin(i)-Pout(i);  % Battery is initially emtpy hence 0
            else
                B(i) = B(i-1)+Pin(i)-Pout(i);
            end

            if(B(i)>= Bmax)
                B(i) = Bmax;     % if condition is never met. 
            end
        end
    end
    Ro = 1;
    U_logical = U<Ro;
    U_bar = mean(U_logical(:,2:end),2);
    final_average(j) = mean(U_bar);
    %final_results = sum(sum(U_logical,1))/(repeats* N)
end
semilogy(Pin_avg,final_average,'o-')
grid on
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear all;
%rng('shuffle') % For reproducibility
N = 10^4;
Pin_avg = [-5:5:40];
Pin_avg_dB = 10.^(Pin_avg/10);
repeats = 10^4;
U = zeros(repeats,N);
for j = 1:length(Pin_avg_dB)
    for  k = 1:repeats

        Bmax = 200*Pin_avg_dB(j);
        B=zeros(1,N);
        Pout = zeros(1,N);


        gauss1 = normrnd(0,1,[1,N]);
        gauss2 = normrnd(0,1,[1,N]);
        gamma = gauss1.^2+gauss2.^2;

        Pin = exprnd(Pin_avg_dB(j),[1,N]);

        for i =1:N

            if (i == 1)
                Pout(i) = min([0,Pin_avg_dB(j)]); % Battery is initially emtpy hence 0
            else
                Pout(i) = min([B(i-1),Pin_avg_dB(j)]);
            end

            U(k,i) = log2(1+(Pout(i)*gamma(i)));

            if (i == 1)
                B(i) = 0+Pin(i)-Pout(i);  % Battery is initially emtpy hence 0
            else
                B(i) = B(i-1)+Pin(i)-Pout(i);
            end

            if(B(i)>= Bmax)
                B(i) = Bmax;     % if condition is never met. 
            end
        end
    end
    Ro = 1;
    U_logical = U<Ro;
    U_bar = mean(U_logical(:,2:end),2);
    final_average(j) = mean(U_bar);
    %final_results = sum(sum(U_logical,1))/(repeats* N)
end
semilogy(Pin_avg,final_average,'o-')
grid on
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 10^7;
Pin_avg = [-5:5:40];
Pin_avg_dB = 10.^(Pin_avg/10);
U_bar = zeros(size(Pin_avg_dB));

for j = 1:length(Pin_avg_dB)
    
    Bmax = inf;
    B=zeros(1,N);
    Pout = zeros(1,N);
    U = zeros(1,N);
    
    gauss1 = normrnd(0,1,[1,N]);
    gauss2 = normrnd(0,1,[1,N]);
    gamma = gauss1.^2+gauss2.^2;

    Pin = exprnd(Pin_avg_dB(j),[1,N]);
    %x = linspace(min(Pin),max(Pin),50);
    %hist(Pin,x),title('Exponential Distibution of Pin')

    for i =1:N

        Pout(i) = Pin_avg_dB(j);

        U(i) = log2(1+Pout(i)*gamma(i));

        if (i == 1)
            B(i) = 0+Pin(i)-Pout(i);
        else
            B(i) = B(i-1)+Pin(i)-Pout(i);
        end

        if(B(i)>= Bmax)
            B(i) = Bmax;
        end
    end

    Ro = 1;
    U_logical = U<Ro;
    U_bar(j) = mean(U_logical);
end
semilogy(Pin_avg,U_bar,'o-')
grid on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xlabel('$\bar{P_{in}}$ (in dB)','Interpreter','Latex')
ylabel('$\bar{U}$ (Outage probability)','Interpreter','Latex')
title('Figure 1 - Outage Probability for Point-to-Point EH and non-EH')

h = legend('EH, N = ${10}^{2}$','EH, N = ${10}^{3}$', 'EH, N = ${10}^{4}$','non-EH, N $\rightarrow\infty$');
set(h,'interpreter','Latex')
legend('Location','northeast')
% Completed 30/7/2019