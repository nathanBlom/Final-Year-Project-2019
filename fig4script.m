clc;clear all;close all;
format long g

Pin_avg = 10;
gamma1 = 1;
gamma2 = 10;
alpha = linspace(0,1,10^5);

R1 = log2(1 + (alpha.*Pin_avg.*gamma1)./((1-alpha).*Pin_avg.*gamma1 + 1));

R2 = log2(1 + (1-alpha).*Pin_avg*gamma2);

plot(R1,R2,'r','LineWidth',2)

%%%%%%%%%%%%
N = 10^5;
Bmax = 200*Pin_avg;
alpha = linspace(0,1,200);
R1_avgs = zeros(1,length(alpha));
R2_avgs = zeros(1,length(alpha));
repeats = 30;
% Assume battery is initially empty
for k=1:length(alpha)
    R1_repeats = zeros(1,repeats);
    R2_repeats = zeros(1,repeats);
    for j = 1:repeats
        B = zeros(1,N);
        Pin = exprnd(Pin_avg,[1,N]);        % should this change for every alpha?
        Pout_1 = zeros(1,N);
        Pout_2 = zeros(1,N);
        R1 = zeros(1,N);
        R2 = zeros(1,N);

        for i = 1:N
            if i == 1
                Pout_1(i) = 0;
                Pout_2(i) = 0;
            else
                if B(i-1) >= alpha(k)*Pin_avg
                    Pout_1(i) = alpha(k)*Pin_avg;
                else
                    Pout_1(i) = B(i-1);
                end

                if B(i-1) >= Pin_avg
                    Pout_2(i) = (1-alpha(k))*Pin_avg;
                else
                    Pout_2(i) = B(i-1) - Pout_1(i);
                end
            end

            %Update battery capacity
            if i == 1
                B(i) = 0 + Pin(i) - Pout_1(i) - Pout_2(i);
            else
                B(i) = B(i-1) + Pin(i) - Pout_1(i) - Pout_2(i);
            end
            
            if(B(i)>= Bmax)
                B(i) = Bmax;     % if condition is never met. 
            end
            
            R1(i) = log2(1 + (Pout_1(i)*gamma1)/(Pout_2(i)*gamma1 + 1));
            R2(i) = log2(1 + Pout_2(i)*gamma2);
        end
        R1_repeats(j) = mean(R1);
        R2_repeats(j) = mean(R2);
    end
    R1_avgs(k) = mean(R1_repeats);
    R2_avgs(k) = mean(R2_repeats);
end
        
    
hold on
plot(R1_avgs,R2_avgs,'m-.')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 10^4;
Bmax = 200*Pin_avg;
alpha = linspace(0,1,200);
R1_avgs = zeros(1,length(alpha));
R2_avgs = zeros(1,length(alpha));
repeats = 100;
% Assume battery is initially empty
for k=1:length(alpha)
    R1_repeats = zeros(1,repeats);
    R2_repeats = zeros(1,repeats);
    for j = 1:repeats
        B = zeros(1,N);
        Pin = exprnd(Pin_avg,[1,N]);        % should this change for every alpha?
        Pout_1 = zeros(1,N);
        Pout_2 = zeros(1,N);
        R1 = zeros(1,N);
        R2 = zeros(1,N);

        for i = 1:N
            if i == 1
                Pout_1(i) = 0;
                Pout_2(i) = 0;
            else
                if B(i-1) >= alpha(k)*Pin_avg
                    Pout_1(i) = alpha(k)*Pin_avg;
                else
                    Pout_1(i) = B(i-1);
                end

                if B(i-1) >= Pin_avg
                    Pout_2(i) = (1-alpha(k))*Pin_avg;
                else
                    Pout_2(i) = B(i-1) - Pout_1(i);
                end
            end

            %Update battery capacity
            if i == 1
                B(i) = 0 + Pin(i) - Pout_1(i) - Pout_2(i);
            else
                B(i) = B(i-1) + Pin(i) - Pout_1(i) - Pout_2(i);
            end
            
            if(B(i)>= Bmax)
                B(i) = Bmax;     % if condition is never met. 
            end
            
            R1(i) = log2(1 + (Pout_1(i)*gamma1)/(Pout_2(i)*gamma1 + 1));
            R2(i) = log2(1 + Pout_2(i)*gamma2);
        end
        R1_repeats(j) = mean(R1);
        R2_repeats(j) = mean(R2);
    end
    R1_avgs(k) = mean(R1_repeats);
    R2_avgs(k) = mean(R2_repeats);
end
        
    
hold on
plot(R1_avgs,R2_avgs,'m-.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 10^3;
Bmax = 200*Pin_avg;
alpha = linspace(0,1,200);
R1_avgs = zeros(1,length(alpha));
R2_avgs = zeros(1,length(alpha));
repeats = 1000;
% Assume battery is initially empty
for k=1:length(alpha)
    R1_repeats = zeros(1,repeats);
    R2_repeats = zeros(1,repeats);
    for j = 1:repeats
        B = zeros(1,N);
        Pin = exprnd(Pin_avg,[1,N]);        % should this change for every alpha?
        Pout_1 = zeros(1,N);
        Pout_2 = zeros(1,N);
        R1 = zeros(1,N);
        R2 = zeros(1,N);

        for i = 1:N
            if i == 1
                Pout_1(i) = 0;
                Pout_2(i) = 0;
            else
                if B(i-1) >= alpha(k)*Pin_avg
                    Pout_1(i) = alpha(k)*Pin_avg;
                else
                    Pout_1(i) = B(i-1);
                end

                if B(i-1) >= Pin_avg
                    Pout_2(i) = (1-alpha(k))*Pin_avg;
                else
                    Pout_2(i) = B(i-1) - Pout_1(i);
                end
            end

            %Update battery capacity
            if i == 1
                B(i) = 0 + Pin(i) - Pout_1(i) - Pout_2(i);
            else
                B(i) = B(i-1) + Pin(i) - Pout_1(i) - Pout_2(i);
            end
            
            if(B(i)>= Bmax)
                B(i) = Bmax;     % if condition is never met. 
            end
            
            R1(i) = log2(1 + (Pout_1(i)*gamma1)/(Pout_2(i)*gamma1 + 1));
            R2(i) = log2(1 + Pout_2(i)*gamma2);
        end
        R1_repeats(j) = mean(R1);
        R2_repeats(j) = mean(R2);
    end
    R1_avgs(k) = mean(R1_repeats);
    R2_avgs(k) = mean(R2_repeats);
end
        
    
hold on
plot(R1_avgs,R2_avgs,'m-.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 10^2;
Bmax = 200*Pin_avg;
alpha = linspace(0,1,200);
R1_avgs = zeros(1,length(alpha));
R2_avgs = zeros(1,length(alpha));
repeats = 2000;
% Assume battery is initially empty
for k=1:length(alpha)
    R1_repeats = zeros(1,repeats);
    R2_repeats = zeros(1,repeats);
    for j = 1:repeats
        B = zeros(1,N);
        Pin = exprnd(Pin_avg,[1,N]);        % should this change for every alpha?
        Pout_1 = zeros(1,N);
        Pout_2 = zeros(1,N);
        R1 = zeros(1,N);
        R2 = zeros(1,N);

        for i = 1:N
            if i == 1
                Pout_1(i) = 0;
                Pout_2(i) = 0;
            else
                if B(i-1) >= alpha(k)*Pin_avg
                    Pout_1(i) = alpha(k)*Pin_avg;
                else
                    Pout_1(i) = B(i-1);
                end

                if B(i-1) >= Pin_avg
                    Pout_2(i) = (1-alpha(k))*Pin_avg;
                else
                    Pout_2(i) = B(i-1) - Pout_1(i);
                end
            end

            %Update battery capacity
            if i == 1
                B(i) = 0 + Pin(i) - Pout_1(i) - Pout_2(i);
            else
                B(i) = B(i-1) + Pin(i) - Pout_1(i) - Pout_2(i);
            end
            
            if(B(i)>= Bmax)
                B(i) = Bmax;     % if condition is never met. 
            end
            
            R1(i) = log2(1 + (Pout_1(i)*gamma1)/(Pout_2(i)*gamma1 + 1));
            R2(i) = log2(1 + Pout_2(i)*gamma2);
        end
        R1_repeats(j) = mean(R1);
        R2_repeats(j) = mean(R2);
    end
    R1_avgs(k) = mean(R1_repeats);
    R2_avgs(k) = mean(R2_repeats);
end
        
    
hold on
plot(R1_avgs,R2_avgs,'m-.')
legend('Capacity Region of non-EH broadcast network','Rate Region of EH broadcast network');
title('Figure 4 - Comparison between Rate Region and Capacity region for different N')
xlabel('$\bar{R_{1}}$ (bits/symb)','Interpreter','Latex')
ylabel('$\bar{R_{2}}$ (bits/symb)','Interpreter','Latex')