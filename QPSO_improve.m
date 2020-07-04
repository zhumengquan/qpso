%ʵ��qpso��  ��������ϵ�������˷�����ϵ��
%�޸ĸ��¹�ʽ
%����һ���Ż����� Эͬ�������ϵڶ���
clear
clc;
%%��ʼ����Ⱥ
N = 500 ;  %����Ⱥ��ģ
d = 12 ;  % ����ά��,���޲��ϵ������
T = 2000 ; %����������
limit = [0.5 1 ; -0.5 0; 0 0.5; -0.5 0; 0 0.1; -0.1 0; 0 0.1; -0.1 0;...
    0 0.1; -0.1 0; 0 0.1; -0.1 0; 0 0.1; -0.1 0; 0 0.1; -0.1 0];% ����λ�ò�������
range=0.0005; % ������Χ
cutoff_k= [ 0 0.50 0.65 0.72 0.78 0.82]; % ��ֹ����
vlimit=[-0.00005 0.00005];

% for i=1:d
%     for j=1:N
%         x(j,i)=limit(i,1)+rand*(limit(i,2)-limit(i,1));%��ʼ����Ⱥλ��
%     end
% end


qaqa=zeros(N,d);
errorBest=zeros(6,T);%��¼�����̶ȣ����12ά
st=zeros(6,12);%���ϵ�������12ά
gbest=zeros(T,d);%��¼ÿ��gBest������
x=zeros(N,d);%��Ⱥ�����ӵ�λ��
temp=zeros(N,d);
fx=zeros(N,1);


    for i=1:d
        x(:,i)=limit(i,1)+rand(N,1)*(limit(i,2)-limit(i,1));
    end
    
pBest=x;%������ʷ���λ��
 pBestValue=zeros(1,N);
    for i=1:N
        pBestValue(i)=func_fitness(pBest(i),cutoff_k(d/2));
    end

   
[gBestValue,index]  = min(pBestValue);  
 gBest = x(index,:);
% gBest=zeros(1,d);
% gBestValue=-inf;
% for i=1:N
%     pBestValue(i)=func_fitness(pBest(i),cutoff_k(d/2));%������ʷ��Ѻ���ֵ
% end
% [gBestValue,index] = min(pBestValue) ;
% gBest = x(index,:);%Ⱥ����ʷ���λ��
%��ʼ�����

%%
    
       steplength1=0.2; 
       steplength2=0.002;
    
t=1            %��������
while(t<=T)
    qaq=0;
    
    qqq=0;

    
    % ����������ʷ����
    mbest =sum(pBest)/N;        %����Ⱥƽ��ֵ
    b = exp(-(t/T)^2.3);          %��������ϵ�� ��ֵδ��Ĭ��Ϊ1��Ѱ�����ŵ�n  2.3 2.5
    %b = (1-0.5)*(T-t)/T+0.5;
    a = rand;
    u = rand(1,d);
    qaqa(t,:)=log(1./u);
    m = rand;
    
    for i=1:N
        
        % ���¸����λ��
        p = a.*pBest(i,:)+(1-a).*gBest;%���¸���λ��
        
        if (m>=0.5)
            
            temp(i,:) = p - b*(mbest-x(i,:));%.*log(1./u);
        else 
            
            temp(i,:) = p + b*(mbest-x(i,:));%.*log(1./u);
        end
%         x(i,:) = p + b*(mbest-x(i,:)).*(1-2*(u >= 0.5));
        % �߽紦������������Χ��ȡ�÷�Χ��ֵ
        for di=1:d    %���Ʋ���
            if di<5
                if (temp(i,di)-x(i,di)<-steplength1)
                    x(i,di)=x(i,di)-steplength1;
%                     qqq=qqq+1;
                elseif(temp(i,di)-x(i,di)>steplength1)
                    x(i,di)=x(i,di)+steplength1;
%                     qaq=qaq+1;
                else
                    x(i,di)=temp(i,di);
                end
            else
                if (temp(i,di)-x(i,di)<-steplength2)
                    x(i,di)=x(i,di)-steplength2;
%                     qqq=qqq+1;
                elseif(temp(i,di)-x(i,di)>steplength2)
                    x(i,di)=x(i,di)+steplength2;
%                     qaq=qaq+1;
                else
                    x(i,di)=temp(i,di);
                end
            end
        end
        
        for di=1:d %di�����ӵĵ�diά���꣬d�����ӵ���ά��
            if(x(i,di)>limit(di,2))
                x(i,di) = limit(di,2);
            end
            if(x(i,di)<limit(di,1))
                x(i,di) = limit(di,1);
            end
        end
        if(func_fitness(x(i,:),cutoff_k(d/2))< pBestValue(i))
            pBest(i,:) = x(i,:) ;
            pBestValue(i) = func_fitness(x(i,:),cutoff_k(d/2));
        elseif(pBestValue(i) < gBestValue)
            gBest = pBest(i,:) ;
            gBestValue = pBestValue(i) ;
        end
%         if(func_fitness(x(i,:),cutoff_k(d/2))< pBestValue(i))
%             pBest(i,:) = x(i,:) ;
%             pBestValue(i) = func_fitness(x(i,:),cutoff_k(d/2));
%         end
%         if(pBestValue(i) < gBestValue)
%             gBest = pBest(i,:) ;
%             gBestValue = pBestValue(i) ;
%         end
    end
  
    

   
%     for i=1:N   %Эͬ�㷨
%         for j=1:d
%             cooperate=gBest;
%             cooperate(d)=x(i,d);
%         end
%         for k=1:N
%             if(func_fitness(cooperate,cutoff_k(d/2))< gBestValue)
%                 gBest =cooperate ;
%                 gBestValue = func_fitness(cooperate,cutoff_k(d/2));
%             end
%         end
%     end

    % ÿ��ȫ�����ŵ����ֵ
    errorBest(d/2,t) = detect_error(gBest,cutoff_k(d/2),range); 
    gbest(t,:)=gBest; %�洢ÿ�����Ž�����
    t=t+1 %��ʾ����
end
st(d/2,:)=[gBest,zeros(1,12-length(gBest))];%�洢���ϵ��������������ȫ��������������

%% Ƶɢ����
figure(1);
hold on
plot(errorBest(2,:),'g');
plot(errorBest(3,:),'r');
plot(errorBest(4,:),'b');
plot(errorBest(5,:),'m');
plot(errorBest(6,:),'k');
%legend('PSO N=8','PSO N=12','PSO N=16','PSO N=20','PSO N=24','Location','NorthEast');
xlabel('��������') ;
ylabel('��Ӧ��ֵ') ;
title('��Ӧ�Ƚ�������') ;

%Et ����Ⱥ  comEt ���� Ewi��ͳ
Nwin8=9;Nwin12=13;Nwin16=17;Nwin20=21;Nwin24=25;
comst8=[0.854934212837883, -0.262500788026354, 0.072072159204863, -0.012021800384911];
comst12=[0.921761796961476, -0.359745370658747, 0.156644793001838, -0.062448651477732, 0.020265229051348, -0.004410000732067];
comst16=[0.952126075962408, -0.410557251909783, 0.213152450061826, -0.111676426281561, 0.055240459672677, -0.024535628158168,...
    0.009141004208080, -0.002534620380839];
comst24=[0.977204107547906, -0.455896766993109, 0.270678273753916, -0.172431138897995, 0.111576602965736, -0.071426557323346,...
    0.044464050451698, -0.026515441785330, 0.014889980682154, -0.007686282801999, 0.003500427384046, -0.001418797372744];
swi8=[0.800000000000000, -0.200000000000000, 0.038095238095238, -0.003571428571429]; % ����ʽ�������޲��ϵ��
swi12=[0.857142857142857, -0.267857142857143, 0.079365079365079, -0.017857142857143, 0.002597402597403, -1.803751803751804e-04];
swi16=[0.888888888888889, -0.311111111111111, 0.113131313131313, -0.035353535353535, 0.008702408702409, -0.001554001554002,...
    1.776001776001776e-04, -9.712509712509713e-06];
swi24=[0.923076923076923, -0.362637362637363, 0.161172161172161, -0.067994505494505, 0.025597931480284, -0.008295625942685,...
    0.002245432585990, -4.911883781852822e-04, 8.316416985147635e-05, -1.020651175449937e-05, 8.068388738734680e-07, -3.081676254377829e-08];
k=0:0.001:1;
stt8=0;
stt12=0;
stt16=0;
stt20=0;
stt24=0;
comstt8=0;
comstt12=0;
comstt16=0;
comstt24=0;
swit8=0;
swit12=0;
swit16=0;
swit24=0;
for n=1:(Nwin8-1)/2
    stt8=stt8+st(2,n)*sin(n*k*pi);
    comstt8=comstt8+comst8(n)*sin(n*k*pi);
    swit8=swit8+swi8(n)*sin(n*k*pi);
end
for n=1:(Nwin12-1)/2
    stt12=stt12+st(3,n)*sin(n*k*pi);
    comstt12=comstt12+comst12(n)*sin(n*k*pi);
    swit12=swit12+swi12(n)*sin(n*k*pi);
end
for n=1:(Nwin16-1)/2
    stt16=stt16+st(4,n)*sin(n*k*pi);
    comstt16=comstt16+comst16(n)*sin(n*k*pi);
    swit16=swit16+swi16(n)*sin(n*k*pi);
end
for n=1:Nwin20/2
    stt20=stt20+st(5,n)*sin(n*k*pi);
end
for n=1:(Nwin24-1)/2
    stt24=stt24+st(6,n)*sin(n*k*pi);
    comstt24=comstt24+comst24(n)*sin(n*k*pi);
    swit24=swit24+swi24(n)*sin(n*k*pi);
end
Et=(2.*stt8)-(k*pi);
Et12=(2.*stt12)-(k*pi);
Et16=(2.*stt16)-(k*pi);
Et20=(2.*stt20)-(k*pi);
Et24=(2.*stt24)-(k*pi);

comEt8=(2.*comstt8)-(k*pi);
comEt12=(2.*comstt12)-(k*pi);
comEt16=(2.*comstt16)-(k*pi);
comEt24=(2.*comstt24)-(k*pi);

Ewi8=(2.*swit8)-(k*pi);
Ewi12=(2.*swit12)-(k*pi);
Ewi16=(2.*swit16)-(k*pi);
Ewi24=(2.*swit24)-(k*pi);

figure(2)
hold on
plot(k*1000,Et,'b');%����Ⱥ

plot(k*1000,comEt8,'r');%����
plot(k*1000,Ewi8,'k');%��ͳ
h1 = refline(0,0);

set(h1,'color','k','linestyle','--');
legend('QPSO','Combination window','Conventional');

figure(3)
hold on
plot(k*100,Et12,'b');plot(k*100,comEt12,'r');plot(k*100,Ewi12,'k');
h1 = refline(0,0);
set(h1,'color','k','linestyle','--');

figure(4)
hold on
plot(k*100,Et16,'b');plot(k*100,comEt16,'r');plot(k*100,Ewi16,'k');
h1 = refline(0,0);
set(h1,'color','k','linestyle','--');

figure(5)
hold on
plot(k*100,Et20,'b');
%plot(k*100,comEt20,'r');plot(k*100,Ewi20,'k');
h1 = refline(0,0);
set(h1,'color','k','linestyle','--');

figure(6)
hold on
plot(k*100,Et24,'b');plot(k*100,comEt24,'r');plot(k*100,Ewi24,'k');
h1 = refline(0,0);
set(h1,'color','k','linestyle','--');

legend('QPSO','Combination window','Conventional');





