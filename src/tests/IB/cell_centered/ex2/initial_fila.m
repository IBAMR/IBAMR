%%%Contributed by Yongsam Kim.

%%%initial

%% Choose the following parameters for a filament to determined;
length=3; % filament length
M=200; % number of points on the filament (length/M < 0.5*9/gridnumber)
freq_b=2.0;  amp_a=0.01; % initial shape of the filament (sin curve);
fix=[4.5 14.75]; %% fixed point
%%% ================================================

Para_t=(0:10000); Para_s=Para_t;
inter_S=zeros(1,M+1); L_X=zeros(2,M+1); Tan=L_X;
leng=length/2.0;
for k=1:1000000
  length1=leng;
  leng=sum(sqrt(1.0+amp_a^2*freq_b^2*cos(freq_b*leng*Para_t/10000).^2));
  leng=length*10000/leng;
  if (abs(length1-leng)<0.00000000000001), break, end
end
Para_t=linspace(0,leng,10001); Para_s(1)=0.0;
for i=1:10000
  Para_s(i+1)=Para_s(i)+(0.5)*(leng/10000) ...
             *(sqrt(1+amp_a^2*freq_b^2*cos(Para_t(i))^2)  ...
             +sqrt(1+amp_a^2*freq_b^2*cos(Para_t(i+1))^2));
end
inter_S=linspace(0,Para_s(10001),M+1);
L_X(1,M+1)=leng;
for j=1:10000; for k=1:M
  if ((Para_s(j)<=inter_S(k)) & (inter_S(k)<Para_s(j+1)))
     L_X(1,k)=(Para_t(j+1)-Para_t(j))*(inter_S(k)-Para_s(j+1))/ ...
              (Para_s(j+1)-Para_s(j))+Para_t(j+1);
end;end;end

%%% Output: L_X which is the coordinates of the curve
%%% 1st row: x, 2nd row: y
inter_S=amp_a*sin(freq_b*L_X(1,:));
L_X(2,:)=-L_X(1,:); L_X(1,:)=inter_S;
L_X(1,:)=L_X(1,:)+fix(1);L_X(2,:)=L_X(2,:)+fix(2);
Tan=L_X(:,2:M+1)-L_X(:,1:M);
ds=sum(sqrt(Tan(1,:).^2+Tan(2,:).^2))/M;
plot(L_X(1,:),L_X(2,:),'-o')  %% check by figure
[Para_s(10001) ds length/M] %% Check the gridsize : "Para_s(10001)=length" "ds=length/M"
