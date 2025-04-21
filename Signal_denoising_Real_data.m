clc
clear all
close all




dataset=1; %1 for Ausgrid 2023, 2 for French current data

if dataset==1
     Data=xlsread("Data_2023.xlsx");
    % reshape into vector
    Data = reshape(Data',[],1);
    data=Data(1:10000);
else
    Data=xlsread('Data_France.xlsx');
    data=Data(1:10000)/5000;
end


x=data;
N=size(x,1);
DCT_basis=dct(eye(N),"Type",2);
IDCT_basis=idct(eye(N));

if dataset==2
w = 0.06 * randn(N,1);      % w : zero-mean Gaussian noise
y = x + w;                % y : noisy speech signal
else
w = 0.2 * randn(N,1);      % w : zero-mean Gaussian noise
y = x + w;                % y : noisy speech signal

end

x_dct=DCT_basis*x;
x0 = zeros(N,1)+0.01;


missing_ratio=0.60;

all_re=[];
all_time=[];
all_iter=[];
all_SNR=[];

for time=1:1



missing = randperm(N, ceil(N*missing_ratio));
Acommon = intersect([1:N],missing);
remaining=setxor([1:N],Acommon);

Phi = eye(N);
Phi(missing,:)=[];

b=Phi*y;   % Measurements
% b=y./amp;
% A=IDCT_basis;
A = Phi*IDCT_basis; % 

x0 = zeros(N,1)+0.01;



x_ref=x_dct;

lambda=0.1;
pm.lambda = lambda;
pm.delta = pm.lambda*100;
pm.xg = x_ref;
pm.maxit=3000;
pm.reltol = 1e-6;




pm_BDR=pm;
pm_BDR.delta=1/20; 
pm_BDR.gamma0=0.3;
pm_BDR.nu=1.4;
pm_BDR.ratio=10;





fprintf('Running BDR ...\n');
[x_BDR,outputBDR] = CS_L1L2_uncon_BDR(A,b,pm_BDR,1);




denoised_BDR=idct(x_BDR);


time=[outputBDR.time];

iter=[size(outputBDR.err,2) ];
err=[outputBDR.err(end)];


SNR_all=[SNR(x,denoised_BDR)];
RE_all=[RE_cal(x,denoised_BDR)];

all_iter=[all_iter;iter];
all_re=[all_re;RE_all];
all_time=[all_time;time];
all_SNR=[all_SNR;SNR_all];


end

ALL_MEAN=[mean(all_time),mean(all_iter),mean(all_re),mean(all_SNR)];


figure;
subplot(4,1,1)
plot(1:N,x,'b','LineWidth',1.2);
grid on
axis tight
title('$\mathbf{u}$','Interpreter', 'latex','FontSize',18)

subplot(4,1,2)
% plot(1:N,y,'b','LineWidth',1.2,'Marker','x','MarkerEdgeColor','red','MarkerIndices',remaining)
plot(1:N,y,'b','LineWidth',1.2)
grid on
axis tight
title('Noisy $\mathbf{u}$','Interpreter', 'latex','FontSize',18)

index_plot=[1:N];
index_plot(missing)=NaN;
value_plot=y;
value_plot(missing)=NaN;

subplot(4,1,3)
plot(index_plot,value_plot,'b','LineWidth',2);
grid on
title('Sampled entries','Interpreter', 'latex','FontSize',18)
% axis tight

subplot(4,1,4)
plot(1:N,denoised_BDR,'b','LineWidth',1.5)
grid on
axis tight
title('$\mathbf{\widehat{u}}$','Interpreter', 'latex','FontSize',18)

hold off


%%

function RE=RE_cal(true,pred)

RE=norm(true-pred,2)/norm(true,2);
end




function R=SNR(true,denoised)

R=20*log10(norm(true,2)/norm(true-denoised,2));
end




