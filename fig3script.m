clc;clear all;close all;
warning('off','all')

Pin_avg = [-25:5:20];
Pin_avg_dB = 10.^(Pin_avg/10);
lambda_vec = zeros(1,length(Pin_avg));

syms x lambda Pc epsilon
f = ((1./Pin_avg_dB).*exp(-x./Pin_avg_dB));
first_integral = int(Pc*f,0,Inf);
 
second_integral = int((epsilon/lambda)*f,0,Inf);

 % numerical integration
for i=1:length(Pin_avg_dB)
    
    epsilon = 5;
    fun = @(x) (epsilon.*(1/Pin_avg_dB(i))*exp(-x./Pin_avg_dB(i))./x);
    third_integral = integral(fun,0,Inf);


    epsilon = 5; Pc = 10^(-25/10);
    lambda_vec(i) = epsilon/(Pin_avg_dB(i)-Pc+epsilon*third_integral);
end
lambda_vec;
%%%%%
N = 10^7; Pc = 10^(-25/10); epsilon = 5;
U_bar = zeros(size(Pin_avg_dB));
P_d = zeros(size(Pin_avg_dB));

for j = 1:length(Pin_avg_dB)
    lambda = lambda_vec(j);
    Bmax = inf;
    B=zeros(1,N);
    Pout = zeros(1,N);
    U = zeros(1,N);
    
    gauss1 = normrnd(0,1,[1,N]);
    gauss2 = normrnd(0,1,[1,N]);
    gamma = gauss1.^2+gauss2.^2;

    Pin = exprnd(Pin_avg_dB(j),[1,N]);

    for i =1:N
        if ((gamma(i))>lambda)
            P_d(i) = Pc + epsilon*((1/lambda)-(1/gamma(i)));
        else
            P_d(i) = 0;
        end

        if (i == 1)
                Pout(i) = min([0,P_d(i)]); % Battery is initially emtpy hence 0
        else
                Pout(i) = min([B(i-1),P_d(i)]);
        end
        
        if ((Pout(i)-Pc))<0
            
            U(i) = log2(1+(1/epsilon)*(0)*gamma(i));
        else
            U(i) = log2(1+(1/epsilon)*(Pout(i)-Pc)*gamma(i));
        end


        if (i == 1)
            B(i) = 0+Pin(i)-Pout(i);
        else
            B(i) = B(i-1)+Pin(i)-Pout(i);
        end

        if(B(i)>= Bmax)
            B(i) = Bmax;
        end
    end

    U_bar(j) = mean(U);
end
semilogy(Pin_avg,U_bar,'k-')
grid on

hold on

%%% 2Pin
N = 10^4;Pc = 10^(-25/10); epsilon = 5;
repeats = 10^3;
U = zeros(repeats,N);
for j = 1:length(Pin_avg_dB)
    lambda = lambda_vec(j);
    for  k = 1:repeats

        Bmax = 2*Pin_avg_dB(j);
        B=zeros(1,N);
        Pout = zeros(1,N);

        gauss1 = normrnd(0,1,[1,N]);
        gauss2 = normrnd(0,1,[1,N]);
        gamma = gauss1.^2+gauss2.^2;

        Pin = exprnd(Pin_avg_dB(j),[1,N]);

        for i =1:N
            
            if ((gamma(i))>lambda)
                P_d(i) = Pc + epsilon*((1/lambda)-(1/gamma(i)));
            else
                P_d(i) = 0;
            end

            if (i == 1)
                Pout(i) = min([0,P_d(i)]); % Battery is initially emtpy hence 0
            else
                Pout(i) = min([B(i-1),P_d(i)]);
            end
        
            if ((Pout(i)-Pc))<0
                 U(k,i) = log2(1+(1/epsilon)*(0)*gamma(i));
            else
                 U(k,i) = log2(1+(1/epsilon)*(Pout(i)-Pc)*gamma(i));
            end

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
    final_average(j) = mean(mean(U));
    %final_results = sum(sum(U_logical,1))/(repeats* N)
end
semilogy(Pin_avg,final_average,'ro-')
grid on


%%% 20Pin
N = 10^4;Pc = 10^(-25/10); epsilon = 5;
repeats = 10^3;
U = zeros(repeats,N);
for j = 1:length(Pin_avg_dB)
    lambda = lambda_vec(j);
    for  k = 1:repeats

        Bmax = 20*Pin_avg_dB(j);
        B=zeros(1,N);
        Pout = zeros(1,N);

        gauss1 = normrnd(0,1,[1,N]);
        gauss2 = normrnd(0,1,[1,N]);
        gamma = gauss1.^2+gauss2.^2;

        Pin = exprnd(Pin_avg_dB(j),[1,N]);

        for i =1:N
            
            if ((gamma(i))>lambda)
                P_d(i) = Pc + epsilon*((1/lambda)-(1/gamma(i)));
            else
                P_d(i) = 0;
            end

            if (i == 1)
                Pout(i) = min([0,P_d(i)]); % Battery is initially emtpy hence 0
            else
                Pout(i) = min([B(i-1),P_d(i)]);
            end
        
            if ((Pout(i)-Pc))<0
                 U(k,i) = log2(1+(1/epsilon)*(0)*gamma(i));
            else
                 U(k,i) = log2(1+(1/epsilon)*(Pout(i)-Pc)*gamma(i));
            end

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
    final_average(j) = mean(mean(U));
    %final_results = sum(sum(U_logical,1))/(repeats* N)
end
semilogy(Pin_avg,final_average,'mx-')
grid on

%%% 200Pin

N = 10^4;Pc = 10^(-25/10); epsilon = 5;
repeats = 10^3;
U = zeros(repeats,N);
for j = 1:length(Pin_avg_dB)
    lambda = lambda_vec(j);
    for  k = 1:repeats

        Bmax = 200*Pin_avg_dB(j);
        B=zeros(1,N);
        Pout = zeros(1,N);

        gauss1 = normrnd(0,1,[1,N]);
        gauss2 = normrnd(0,1,[1,N]);
        gamma = gauss1.^2+gauss2.^2;

        Pin = exprnd(Pin_avg_dB(j),[1,N]);

        for i =1:N
            
            if ((gamma(i))>lambda)
                P_d(i) = Pc + epsilon*((1/lambda)-(1/gamma(i)));
            else
                P_d(i) = 0;
            end

            if (i == 1)
                Pout(i) = min([0,P_d(i)]); % Battery is initially emtpy hence 0
            else
                Pout(i) = min([B(i-1),P_d(i)]);
            end
        
            if ((Pout(i)-Pc))<0
                 U(k,i) = log2(1+(1/epsilon)*(0)*gamma(i));
            else
                 U(k,i) = log2(1+(1/epsilon)*(Pout(i)-Pc)*gamma(i));
            end

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
    final_average(j) = mean(mean(U));
    %final_results = sum(sum(U_logical,1))/(repeats* N)
end
semilogy(Pin_avg,final_average,'bs-')
grid on

title('Figure 3 - Average data rate for Point-to-Point EH and equivalent non-EH systems for various battery capacities B_{max}')
xlabel('$\bar{P_{in}}$ (in dB)','Interpreter','Latex')
ylabel('$\bar{U}$ (Rate in bits/symb)','Interpreter','Latex')

h = legend('non-EH, N $\rightarrow\infty$','$B_{max}$=2$\bar{P_{in}}$, N=${10}^{4}$','$B_{max}$=20$\bar{P_{in}}$, N=${10}^{4}$','$B_{max}$=200$\bar{P_{in}}$, N=${10}^{4}$');
set(h,'interpreter','Latex')
legend('Location','northwest')