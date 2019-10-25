clc;clear all;close all;
format long g

Pin_avg = [-20:5:20];
Pin_avg_dB = 10.^(Pin_avg/10);
lambda_vec = zeros(1,length(Pin_avg));

for i=1:length(Pin_avg)
   
    %closed form
    syms x lambda
    f = ((1/Pin_avg_dB(i))*exp(-x/Pin_avg_dB(i))/lambda);
    first_integral = int(f,0,Inf);

    % numerical integration
    fun = @(x) ((1/Pin_avg_dB(i))*exp(-x./Pin_avg_dB(i))./x);
    second_integral = integral(fun,0,Inf);

    % solve for lambda

%     syms lambda
%     expected_value = first_integral - second_integral;
%     a = solve( expected_value == Pin_avg_dB)               % unsure

      lambda_vec(i) = 1/(Pin_avg_dB(i) + second_integral);
end

%%%%%%%%%%%
N = 10^2;Pc = 0; epsilon = 1;
repeats = 10^2;
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
semilogy(Pin_avg,final_average,'ro-')
grid on
hold on

%%%%%%%%%%%%%%%%%%
N = 10^3;Pc = 0; epsilon = 1;
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
semilogy(Pin_avg,final_average,'gx-')
grid on
hold on
%%%%%%%%%%%

N = 10^4;Pc = 0; epsilon = 1;
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
hold on

%%%%%%%%%%%%%%%%%%


N = 10^7; Pc = 0; epsilon = 1;
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

xlabel('$\bar{P_{in}}$ (in dB)','Interpreter','Latex')
ylabel('$\bar{U}$ (Outage probability)','Interpreter','Latex')
title('Figure 2 - Average data rate for Point-to-Point EH and equivalent non-EH systems for various N')

h = legend('EH, N = ${10}^{2}$','EH, N = ${10}^{3}$', 'EH, N = ${10}^{4}$','non-EH, N $\rightarrow\infty$');
set(h,'interpreter','Latex')
legend('Location','northwest')
