
close all;
%Reading excel data
[num,txt,raw]= xlsread('rawNucleolarData.xlsx','Confocal');
%M= num(:,size(num,2));

% %only take AB nucleus7-10, others are empty
% zeroel=find(~(num(:,2)>6 & num(:,2)<11));
% num(zeroel(:),:)=[];
% 
% %trial
% zeroel=find(~(num(:,3)==5 & num(:,1)==1));
% num(zeroel(:),:)=[];

% Emb=num(:,1);
% Nuc=num(:,2);
% Nucl=num(:,3);

time=num(:,1);
mass=num(:,2);


figure;
plot(time,mass,'.','MarkerSize',30);
title('Nucleoli growth data from Steph');
xlabel('time (min)');
ylabel('Mass (AU)');

alpha = 40; %intensity to number conversion
figure;
 plot(time,mass/alpha,'.','MarkerSize',30);
 title('Nucleoli growth data with alpha = 40');
 xlabel('time (min)');
 ylabel('Number');
 
 alpha = 400;
 figure;
 plot(time,mass/alpha,'.','MarkerSize',30);
 title('Nucleoli growth data with a different alpha');
 xlabel('time (min)');
 ylabel('Number');






