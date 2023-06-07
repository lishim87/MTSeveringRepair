%least squares but with theory
cts=[0,2,5];
% Ls=20:40;
Ls=24:30;
kts=.05:.01:.2;
kss=.4:.01:.7;
x0s=1;
lstsqrsm=zeros(length(kss),length(kts),length(Ls),length(x0s));
for w=1:length(x0s)
    x0=x0s(w);
    for z=1:length(Ls)
        L=Ls(z);
        for x=1:length(kts)
            kt=kts(x);
            for y=1:length(kss)
                ks=kss(y);
                kT=kt.*cts;
                pt= kT./(kT+ks);
                ps= ks./(kT+ks);
                r=pt./ps;
                one_vec=ones(size(cts));
                k=ones(size(cts)).*x0;
                n=ones(size(cts)).*L;
                tstep=zeros(size(cts));
                %syms j
                for a=1:length(kT)
                    tstep(a)=1/(ks+kT(a));
                end

                C5=tstep.*(((r+one_vec)./(r-one_vec)).*(((r.^n+one_vec)./(r.^n-one_vec)).*n-((r.^k+one_vec)./(r.^k-one_vec)).*k));
                for i=1:length(kT)
                    if kT(i)==ks
                        C5(i)=(L^2-x0^2)/3;
                    end
                end
                data=[32.3,9.2; %0: mean, std 
                      75.4,30.8; %2: mean, std 
                      193.8,101.5]; %5: mean, std 
                lstsqrm=sum(((data(:,1))'-C5).^2);
                lstsqrsm(y,x,z,w)=lstsqrm;
            end
        end
    end
end

for i=5
figure
% lstsqrsm(isnan(lstsqrsm))=max(max(max(max(lstsqrsm)))); %convert nan values to max
lstsqrsm(isnan(lstsqrsm))=1000;
lstsqrsm(lstsqrsm>1000)=1000;
imagesc('XData',Ls,'YData',kss,'CData',squeeze(lstsqrsm(:,i,:,1)));

colorbar
xlabel('N_D')
ylabel('ks')
title("L="+num2str(Lind))
%%
figure
imagesc('XData',kts,'YData',kss,'CData',squeeze(lstsqrsm(:,:,i,1)));
end
%%
[M,I]= min(lstsqrsm,[],'all','linear');
sz=size(lstsqrsm);
[ksmin,ktmin,Lmin,x0min]=ind2sub(sz,I);
best_ks=kss(ksmin);
best_kt=kts(ktmin);
best_L=Ls(Lmin);
best_x0=x0min;
kt_ks=[.15,.7,30];
kt_ind=find(kts==kt_ks(1));
ks_ind=find(kss==kt_ks(2));
L_ind=find(Ls==kt_ks(3));
err=lstsqrsm(ks_ind,kt_ind,L_ind);
Is=find(lstsqrsm<200);
best_inds=zeros(length(Is),4);
best_vals=zeros(length(Is),5);
vec=zeros(1,4);
for i=1:length(Is)
    [I1,I2,I3,I4]=(ind2sub(sz,Is(i)));
    best_inds(i,:)=[I1,I2,I3,I4];
    best_vals(i,1)=kss(best_inds(i,1));
    best_vals(i,2)=kts(best_inds(i,2));
    best_vals(i,3)=Ls(best_inds(i,3));
    best_vals(i,4)=x0s(best_inds(i,4));
    best_vals(i,5)=lstsqrsm(I1,I2,I3,I4);
end
%  D=best_vals; %stores severing probability for each ct
% fileID2=fopen('best_param_sets.dat','w');
% fprintf(fileID2,'%8.4f %3.8f %3.8f %8.4f %8.4 \r\n',D);
% fclose(fileID2);
save('best_param_sets','best_vals')
%store best vals
%do simulation with stdevs
%%
r=best_kt/best_ks; n=best_L; k=1;
C5=tstep.*(((r+one_vec)./(r-one_vec)).*(((r.^n+one_vec)./(r.^n-one_vec)).*n-((r.^k+one_vec)./(r.^k-one_vec)).*k));
plot(cts,C5)
hold on
plot(cts,data(:,1))