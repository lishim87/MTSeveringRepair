%random walk model for severing competition
%computes squared error between this and vemu
clear all;
close all;
MaxT=1000000;
numTraj=200;
%parameters
kts=.06:.01:.2;
kss=.2:.02:1;
lstsqrsm=zeros(length(kss),length(kts));
lstsqrssd=zeros(length(kss),length(kts));
cts=[0 2 5]; 
ind=2;
for x=1:length(kts)
    kt=kts(x);
    for y=1:length(kss)
        ks=kss(y);
        vals=zeros(length(cts),2);
        for a=ind
            ct=cts(a); %ct is concentration of tubulin
            sevtime=zeros(1,numTraj);
            sevdata=0; %to store number of severing events
            for j=1:numTraj
                t=zeros(1,MaxT);
                dt=zeros(1,MaxT);
                l=zeros(1,MaxT);
                x0=1;
                l(1)=x0;
                L=26;
                for i= 1:MaxT
                    kp=ct*kt; % prob of tubulin adding, based on concentration
                    k0=ks+kp;%rate of something happening is sum of rate of each indepedent event occurring
                    r1=rand;
                    dt(i)=(1/k0)*log(1/r1);
                     if l(i)==0  %if all tubulin spots filled 
                         sever=false;
                         break;
                     end
                    if l(i)==L %if all tubulin gone
                        sever=true;
                         break;
                    end
                    t(i+1)=t(i)+dt(i);
                    p=rand;
                    if p<=(ks/k0)  %a tubulin comes off- go towards severing
                      l(i+1)=l(i)+1;
                    else %repair-fill one unit
                        l(i+1)=l(i)-1;
                    end

                end %for i
                if sever==true
                    sevdata=sevdata+1;
                    sevtime(sevdata)=t(i); %instead add to a vector to make stdev easier
                end
            end %for trajectories
                sevtime=sevtime(1,1:sevdata); %cut off unfilled elements
                vals(a,:)=[mean(sevtime), std(sevtime)];
        end %for concentrations of tubulin
data=[32.3,9.2; %0: mean, std 
      75.4,30.8; %2: mean, std 
      193.8,101.5]; %5: mean, std 
lstsqrm=(data(ind,1)-vals(ind,1))^2;
lstsqrsm(y,x)=lstsqrm;
lstsqrsd=(data(ind,2)-vals(ind,2))^2;
lstsqrssd(y,x)=lstsqrsd;

    end %for ks
end %for kt

ksstore=zeros(1,length(kss)+1);
ksstore(1:length(kss))=kss;
D=[[lstsqrsm;kts],ksstore']; %stores severing probability for each ct
fileID2=fopen('lstsqrsct2tr18m.dat','w');
fprintf(fileID2,'%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\r\n',D');
fclose(fileID2);

E=[[lstsqrssd;kts],ksstore']; %stores severing probability for each ct
fileID3=fopen('lstsqrsct2tr18sd.dat','w');
fprintf(fileID3,'%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f\r\n',E');
fclose(fileID3);
