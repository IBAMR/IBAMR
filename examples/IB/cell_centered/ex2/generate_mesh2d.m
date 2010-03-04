NFINEST = 256;       % Cartesian grid spacing parameters
dx = 7.68/NFINEST;

ds = 0.5*dx;         % Lagrangian mesh spacing parameters
length = 3;          % filament length
M = ceil(3/ds);      % number of points on the filament

freq_b = 2.0;        % initial shape of the filament (frequency of sine curve)
amp_a = 0.01;        % initial shape of the filament (amplitude of sine curve)
fix = [4.5 14.75];   % fixed point
L_X = zeros(2,M+1);  % coordinates of the curve

Para_N = 100000;
Para_t = 0:Para_N;
Para_s = Para_t;
inter_S = zeros(1,M+1);

length1 = 0.0;
leng = length/2.0;
while (abs(length1-leng) > eps)
  length1 = leng;
  leng = sum(sqrt(1.0+amp_a^2*freq_b^2*cos(freq_b*leng*Para_t/Para_N).^2));
  leng = length*Para_N/leng;
end

Para_t = linspace(0,leng,Para_N+1);
Para_s(1) = 0.0;
for i = 1:Para_N
  Para_s(i+1) = Para_s(i)+(0.5)*(leng/Para_N)*(sqrt(1+amp_a^2*freq_b^2* ...
                                                    cos(Para_t(i))^2)+ ...
                                               sqrt(1+amp_a^2*freq_b^2* ...
                                                    cos(Para_t(i+1))^2));
end

inter_S = linspace(0,Para_s(Para_N+1),M+1);
L_X(1,M+1) = leng;
j_start = 1;
for k=1:M
  for j = j_start:Para_N;
    if ((Para_s(j) <= inter_S(k)) & (inter_S(k) <= Para_s(j+1)))
      L_X(1,k) = (Para_t(j+1)-Para_t(j))*(inter_S(k)-Para_s(j+1))/ ...
          (Para_s(j+1)-Para_s(j))+Para_t(j+1);
      j_start = j;
      break;
    end
  end
end

inter_S = amp_a*sin(freq_b*L_X(1,:));
L_X(2,:) = -L_X(1,:);
L_X(1,:) = inter_S;
L_X(1,:) = L_X(1,:) + fix(1);
L_X(2,:) = L_X(2,:) + fix(2);
Tan = L_X(:,2:M+1)-L_X(:,1:M);
ds = sum(sqrt(Tan(1,:).^2+Tan(2,:).^2))/M;
plot(L_X(1,:),L_X(2,:),'-o')     %% plot the filament
[Para_s(Para_N+1) ds length/M]   %% check the gridsize : "Para_s(Para_N+1)=length" "ds=length/M"
