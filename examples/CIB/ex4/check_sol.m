
clear all;
clc;
A = load('Lambda/U.txt');

t95 = 1.737089476; % time to reach 95% of terminal velocity
V = 1.226083506;   % terminal velocity
t = A(:,1);
v = -A(:,3);

figure(1);
plot(t/t95,v/V,'--','LineWidth', 2); hold on; plot(t/t95,1 - exp(-3*t/t95),'-','LineWidth',2); 
xlabel('t/t_{95}');
ylabel('v/V');
legend('CIB solution','Experimental');
hold off;