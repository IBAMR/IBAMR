x=-1:0.0001:1;
y = -1:0.0001:1;
R = 0.75;
MU = 0.2;
%Ro = 0.75*(1.0+1.0/24);
Ro = 0.75 + 0.0625;
omega = 0.000833;
X1 =0.0;
e = 0.0625 - 0.00390625;%0.5* (Ro - R);
c = Ro - R;
eps = e/(Ro-R);

for i =1:max(size(x))
    if (abs(x(i))<R)
        u(i) = omega*x(i);
    elseif (abs(x(i)-e)>Ro)
         u(i) = 0.0;
    else
         u(i) = omega*x(i)*(1.0 - (sqrt(x(i)*x(i) + X1*X1) - R)/(c + e*x(i)/sqrt(x(i)*x(i) + X1*X1)) - (3*eps*((sqrt(x(i)*x(i) + X1*X1) - R)/(c + e*x(i)/sqrt(x(i)*x(i) + X1*X1)) - ((sqrt(x(i)*x(i) + X1*X1) - R)/(c + e*x(i)/sqrt(x(i)*x(i) + X1*X1)))*((sqrt(x(i)*x(i) + X1*X1) - R)/(c + e*x(i)/sqrt(x(i)*x(i) + X1*X1))))*(2*x(i)/sqrt(x(i)*x(i) + X1*X1)+3*eps+eps*eps*x(i)/sqrt(x(i)*x(i) + X1*X1)))/((2+eps^2)*(1+eps*x(i)/sqrt(x(i)*x(i) + X1*X1))));

    end
end


writematrix(u, 'velocity_0.1.csv');
% hold on
 plot(x,u,'-k');

for i =1:max(size(y))
    if (abs(y(i))<R)
        p(i) = 0.5*omega*omega*(y(i)*y(i));
    elseif (abs(y(i))> sqrt(Ro*Ro-e*e))
         p(i) = 0.0;
    else
        p(i) = 6*MU*eps*omega*R*R*2.0*(y(i)/sqrt(y(i)*y(i)))/((2+eps^2)*c*c) ;
    end
end

writematrix(p, 'pressure.csv');
% figure
%plot(y,p,'-k');