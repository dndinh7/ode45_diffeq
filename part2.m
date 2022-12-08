% Demo of ode45 used to solve a system of ODEs, specifically 
%   dS/dt = 0.4 S (1-S/2000) - 0.004 S H
%   dH/dt = -0.1 H + 0.001 S H
% with initial condition (S(0),H(0)) = (500, 10).
time_range = [0 50]; init_cond = [500 10];
[times,solutions] = ode45(@F,time_range,init_cond);
% times is a vector holding the t values from t=0 to t=50
% solutions is a matrix; column 1 holds the values of S(t); 
%                        column 2 holds the values of H(t)
Suscep = solutions(:,1); % notation to pick out 1st column
Infec1 = solutions(:,2); % notation to pick out 2nd column
Phages1 = solutions(:,3);
Infec2 = solutions(:,4); % notation to pick out 2nd column
Phages2 = solutions(:,5);
Nutrient = solutions(:,6);
plot(times,Suscep)
hold on
plot(times,Infec1)
plot(times,Phages1)
plot(times,Infec2)
plot(times,Phages2)
plot(times,Nutrient)
hold off
legend('Suscep', 'Infec1', 'Phages1', 'Infec2', 'Phages2', 'Nutrient')
function output = rhs(t,Y) % defines right-hand-side of ODE, in vector formalism
  Suscep = Y(1); Infec1 = Y(2); Phages1 = Y(3); Infec2 = Y(4); Phages2 = Y(5); Nutrient = Y(6);% input Y is a vector, assigning its slots to scalar variables
  output = [ 0.4*S*(1-S/2000)-0.004*S*H; % output defined here as a column vector
      -0.1*H+0.001*S*H];        % (semicolon on line above moves us to the next row)
end
