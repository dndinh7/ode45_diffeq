%here
clear;

D= 7; % intervention period
d= 10; % dilution factor
N_new= 2*10^5; % reasonable value of N_new

init_cond = [100000 0 100 0 100 2*10^5];
time_range = [0 D];
[times,solutions] = ode45(@F,time_range,init_cond);
times_for_plot= times;
Suscep = solutions(:,1);
Infec1 = solutions(:,2); 
Phages1 = solutions(:,3);
Infec2 = solutions(:,4);
Phages2 = solutions(:,5);
Nutrient = solutions(:,6);


fprintf("At end of period 1) Suscep: %d, Infec1: %d, Phages1: %d, Infec2: %d, Phages2: %d, Nutrient: %d\n", Suscep(end), Infec1(end), Phages1(end), Infec2(end), Phages2(end), Nutrient(end))
for i = 2:4 % number of times to do this process
    init_cond(:,1)= Suscep(end)/d;
    init_cond(:,2)= Infec1(end)/d;
    init_cond(:,3)= Phages1(end)/d;
    init_cond(:,4)= Infec2(end)/d;
    init_cond(:,5) = Phages2(end)/d;
    init_cond(:,6)= Nutrient(end)/d + N_new;


    [times, solutions] = ode45(@F, [D*(i-1) D*i], init_cond);
    Suscep = [Suscep; solutions(:,1)]; % creating new arrays with the values appended
    Infec1 = [Infec1; solutions(:,2)]; 
    Phages1 = [Phages1; solutions(:,3)];
    Infec2 = [Infec2; solutions(:,4)];
    Phages2 = [Phages2; solutions(:,5)];
    Nutrient = [Nutrient; solutions(:,4)]; 
    times_for_plot= [times_for_plot; times];
    fprintf("At end of period %d) Suscep: %d, Infec1: %d, Phages1: %d, Infec2: %d, Phages2: %d, Nutrient: %d\n", i, Suscep(end), Infec1(end), Phages1(end), Infec2(end), Phages2(end), Nutrient(end))
end


% this is plotting the graphs
% Suscep and Nutrient
figure(1);
plot(times_for_plot, Suscep)
hold on
plot(times_for_plot, Nutrient)
title("Susceptible Bacteria and Nutrients over Time")
hold off
legend('Suscep', 'Nutrient')
xlabel('Time')
ylabel('Number of Susceptible Bacteria/Amount of Nutrients')

% Infec1 and Phage1
figure(2);
plot(times_for_plot, Infec1)
hold on
plot(times_for_plot, Phages1)
title("Infected Bacteria and Phage #1 over Time")
hold off
legend('Infec1', 'Phages1')
xlabel('Time')
ylabel('Number of Infected Bacteria/Phages')

% Infec2 and Phage2
figure(3);
plot(times_for_plot, Infec2)
hold on
plot(times_for_plot, Phages2)
title("Infected Bacteria and Phage #2 over Time")
hold off
legend('Infec2', 'Phages2')
xlabel('Time')
ylabel('Number of Infected Bacteria/Phages')


function output = F(t,Y) % defines right-hand-side of ODE, in vector formalism
  S = Y(1); I1 = Y(2); P1 = Y(3); I2 = Y(4); P2 = Y(5); N = Y(6);
  A= 1; % the same right now
  K1= A * 10^(-7);
  K2= 1/A * 10^(-7);
  
  if N > 0
      N= 1.34*N/(N+1);
  else      
      N= 0;
  end

  Suscep= N*S*(1 - (S+I1+I2)/2000000) - K1*S*P1 - K2*S*P2;
  Infec1= K1*S*P1 - 3.3*I1;
  Phage1= -K1*S*P1-2*P1+122.1*I1;
  Infec2= K2*S*P2 - 3.3*I2;
  Phage2= -K2*S*P2-2*P2+122.1*I2;
  Nutrient= -0.1*N*S*(1-(S+I1+I2)/2000000) + 0.0033*(I1+I2);
 

  output = [ Suscep;
      Infec1; 
      Phage1;
      Infec2;
      Phage2; 
      Nutrient];
end

