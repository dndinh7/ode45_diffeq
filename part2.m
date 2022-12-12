clear;
D= 7; %intervention period
d= 10; %dilution factor
N_new= 2*10^5; % reasonable value of N_new

init_cond = [100000 0 100 0 100 2*10^5];
time_range = [0 D];
[times,solutions] = ode45(@F,time_range,init_cond);
% times is a vector holding the t values from t=0 to t=50
% solutions is a matrix; column 1 holds the values of S(t); 
%                        column 2 holds the values of H(t)
times_for_plot= times;
Suscep = solutions(:,1);
Infec1 = solutions(:,2); 
Phages1 = solutions(:,3);
Infec2 = solutions(:,4);
Phages2 = solutions(:,5);
Nutrient = solutions(:,6);


fprintf("At end of period 1) Suscep1: %d, Infec1: %d, Phages1: %d, Infec2: %d, Phages2: %d, Nutrient: %d\n", Suscep(end), Infec1(end), Phages1(end), Infec2(end), Phages2(end), Nutrient(end))
for i = 2:15 % number of times to do this process
    init_cond(:,1)= Suscep(end)/d;
    init_cond(:,2)= Infec1(end)/d;
    init_cond(:,3)= Phages1(end)/d;
    init_cond(:,4)= Infec2(end)/d;
    init_cond(:,5) = Phages2(end)/d;
    init_cond(:,6)= Nutrient(end)/d + N_new;


    [times, solutions] = ode45(@F, [7*(i-1) 7*i], init_cond);
    Suscep = [Suscep; solutions(:,1)];
    Infec1 = [Infec1; solutions(:,2)]; 
    Phages1 = [Phages1; solutions(:,3)];
    Infec2 = [Infec2; solutions(:,4)];
    Phages2 = [Phages2; solutions(:,5)];
    Nutrient = [Nutrient; solutions(:,4)]; % creating new arrays with the values appended
    times_for_plot= [times_for_plot; times];
    fprintf("At end of period %d) Suscep: %d, Infec1: %d, Phages1: %d, Infec2: %d, Phages2: %d, Nutrient: %d\n", i, Suscep(end), Infec1(end), Phages1(end), Infec2(end), Phages2(end), Nutrient(end))
end


% this is plotting the graphs
figure(1);
plot(times_for_plot,Suscep)
hold on
plot(times_for_plot, Nutrient)
hold off

legend('Suscep', 'Nutrient')

figure(2);
plot(times_for_plot, Infec1)
hold on
plot(times_for_plot, Phages1)
hold off

legend('Infec1, Phages1')

figure(3);
plot(times_for_plot, Infec2)
hold on
plot(times_for_plot, Phages2)
hold off
legend('Infec2', 'Phages2')




function output = F(t,Y) % defines right-hand-side of ODE, in vector formalism
  S = Y(1); I1 = Y(2); P1 = Y(3); I2 = Y(4); P2 = Y(5); N = Y(6);% input Y is a vector, assigning its slots to scalar variables
  A= 0.98; % the same right now
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
 

  output = [ Suscep; % output defined here as a column vector
      Infec1; 
      Phage1;
      Infec2;
      Phage2; 
      Nutrient];      % (semicolon on line above moves us to the next row)
end

