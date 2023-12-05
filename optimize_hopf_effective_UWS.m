function optimize_hopf_effective_UWS()

load  empiricalLEiDA.mat;

P1emp=nanmean(P1emp);
P3emp=nanmean(P3emp);

load SC_template.mat;
C=sc_healthy;
C=C/max(max(C))*0.2;

load  fMRI_Liege_2023.mat;

X=tc_gr(:,3)';

TSmax=295;
N=214;
NSUB= 16;
TR=2; % Repetition Time (seconds)
NumClusters= 4; %Number_Clusters;

delt = TR;            % sampling interval
k=2;                  % 2nd order butterworth filter
fnq=1/(2*delt);
flp = .04;           % lowpass frequency of filter
fhi = fnq-0.001;           % highpass
Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt,afilt]=butter(k,Wn);   % construct the filter

flp = .04;           % lowpass frequency of filter
fhi = .07;           % highpass
Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
[bfilt2,afilt2]=butter(k,Wn);   % construct the filter
clear fnq flp fhi Wn k

n_Subjects=16; 


%%%%%%%%%%%%%%
%%% Extracting FC FCD and metastability of data
kk=1;
insub=1;
Isubdiag = find(tril(ones(N),-1));
Tmaxtotal=0;
for nsub=1:n_Subjects
    [N, Tmax0]=size(squeeze(X{1,nsub})); 
    Tmax=min(TSmax,Tmax0);
    Tmaxtotal=Tmaxtotal+Tmax;
    signaldata = squeeze(X{1,nsub});
    signaldata=signaldata(:,1:Tmax);
    Phase_BOLD_data=zeros(N,Tmax);
    timeseriedata=zeros(N,Tmax);
    for seed=1:N
        x=demean(detrend(signaldata(seed,:)));
        x(find(x>3*std(x)))=3*std(x);
        x(find(x<-3*std(x)))=-3*std(x);
        timeseriedata(seed,:) = filtfilt(bfilt2,afilt2,x);    % zero phase filter the data
        Phase_BOLD_data(seed,:) = angle(hilbert(timeseriedata(seed,:)));
    end
    T=10:Tmax-10;
    for t=T
        kudata=sum(complex(cos(Phase_BOLD_data(:,t)),sin(Phase_BOLD_data(:,t))))/N;
        syncdata(t-9)=abs(kudata);
        for i=1:N
            for j=1:i-1
                patt(i,j)=cos(adif(Phase_BOLD_data(i,t),Phase_BOLD_data(j,t)));
            end
        end
        pattern(t-9,:)=patt(Isubdiag);
    end
    metastabilitydata2(nsub)=std(syncdata);
 
    for t=1:Tmax
        for n=1:N
            for p=1:N
                iFC(t,n,p)=cos(Phase_BOLD_data(n,t)-Phase_BOLD_data(p,t));
            end
        end
    end
    FCphasesemp2(nsub,:,:)=squeeze(mean(iFC));
end
FCphasesemp=squeeze(mean(FCphasesemp2));
metastabilitydata=mean(metastabilitydata2);

%%% Extracting peak of data power spectra for determining omega (Hopf)
for nsub=1:n_Subjects
    clear PowSpect PowSpect2;
    [N, Tmax0]=size(squeeze(X{1,nsub})); %size(X{1,nsub});
    Isubdiag = find(tril(ones(N),-1));
    Tmax=min(TSmax,Tmax0);
    TT=Tmax;
    Ts = TT*TR;
    freq = (0:TT/2-1)/Ts;
    signaldata = squeeze(X{1,nsub}); %X{1,nsub};
    signaldata=signaldata(:,1:Tmax);
    FCemp2(nsub,:,:)=corrcoef(signaldata');

    %%%%

    [aux minfreq]=min(abs(freq-0.04));
    [aux maxfreq]=min(abs(freq-0.07));
    nfreqs=length(freq);


    for seed=1:N
        x=detrend(demean(signaldata(seed,:)));
        ts =zscore(filtfilt(bfilt2,afilt2,x));
        pw = abs(fft(ts));
        PowSpect(:,seed,insub) = pw(1:floor(TT/2)).^2/(TT/TR);
        ts2 =zscore(filtfilt(bfilt,afilt,x));
        pw2 = abs(fft(ts2));
        PowSpect2(:,seed,insub) = pw2(1:floor(TT/2)).^2/(TT/TR);
    end
    insub=insub+1;
end

Power_Areas=mean(PowSpect,3);
Power_Areas2=mean(PowSpect2,3);
for seed=1:N
    Power_Areas(:,seed)=gaussfilt(freq,Power_Areas(:,seed)',0.01);
    Power_Areas2(:,seed)=gaussfilt(freq,Power_Areas2(:,seed)',0.01);
    vsig(seed)=sum(Power_Areas2(minfreq:maxfreq,seed))/sum(Power_Areas2(:,seed));
end

vmax=max(vsig);
vmin=min(vsig);

[maxpowdata,index]=max(Power_Areas);
f_diff = freq(index);
FCemp=squeeze(mean(FCemp2));

clear PowSpect PowSpect2 Power_Areas Power_Areas2;

%%%%%%%%%%%%%%%%%%
%% Here we start modelling

omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);

dt=0.1*TR/2;
Tmax=TSmax*n_Subjects;
sig=0.01;
dsig = sqrt(dt)*sig; % to avoid sqrt(dt) at each time step

%%%%%%%%%%%%
%% Optimize
%%
iwe=1;
WE=0:0.01:0.5;  %% G
a=-0.02*ones(N,2);

NWE=length(WE);
PTRsimul=zeros(NWE,NumClusters,NumClusters);
Pstatessimul=zeros(NWE,NumClusters);

for we=WE
    minm=100;
    Cnew=C;
    for iter=1:250 
        wC = we*Cnew;
        sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj
        xs=zeros(Tmax,N);
        %number of iterations, 100 willk�hrlich, weil reicht in diesem Fall
        z = 0.1*ones(N,2); % --> x = z(:,1), y = z(:,2)
        nn=0;
        % discard first 3000 time steps
        for t=0:dt:3000
            suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
        end
        % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
        for t=0:dt:((Tmax-1)*TR)
            suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
            zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
            z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
            if abs(mod(t,TR))<0.01
                nn=nn+1;
                xs(nn,:)=z(:,1)';
            end
        end

        %%%%
        BOLD=xs';
        signal_filt=zeros(N,nn);
        Phase_BOLD=zeros(N,nn);
        for seed=1:N
            BOLD(seed,:)=demean(detrend(BOLD(seed,:)));
            signal_filt(seed,:) =filtfilt(bfilt2,afilt2,BOLD(seed,:));
            Phase_BOLD(seed,:) = angle(hilbert(signal_filt(seed,:)));
        end

        for t=1:nn
            for n=1:N
                for p=1:N
                    iFC(t,n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
                end
            end
        end
        FCphases=squeeze(mean(iFC));

        %% update effective conn matrix Cnew
        for i=1:N
            for j=i+1:N
                if (C(i,j)>0|| j==N/2+i) 
                    Cnew(i,j)=Cnew(i,j)+0.01*(FCphasesemp(i,j)-FCphases(i,j));
                    if Cnew(i,j)<0
                        Cnew(i,j)=0;
                    end
                    Cnew(j,i)=Cnew(i,j);
                end
            end
        end

        Cnew=Cnew/max(max(Cnew))*0.2;

        D = abs(FCphasesemp-FCphases).^2;
        MSE = sum(D(:))/numel(FCphases);
        if MSE<0.001
            break;
        end

        %%%%

    end

    Coptim(iwe,:,:)=Cnew;  %% effective Conn for G (we)

    %%%%%%%%%%%%%%
    %%% Final simul

    xs=zeros(Tmax,N);
    %number of iterations, 100 willk�hrlich, weil reicht in diesem Fall
    z = 0.1*ones(N,2); % --> x = z(:,1), y = z(:,2)
    nn=0;
    % discard first 3000 time steps
    for t=0:dt:3000
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
    end
    % actual modeling (x=BOLD signal (Interpretation), y some other oscillation)
    for t=0:dt:((Tmax-1)*TR)
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(N,2);
        if abs(mod(t,TR))<0.01
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end

    FC_simul=corrcoef(xs(1:nn,:));
    cc=corrcoef(atanh(FCemp(Isubdiag)),atanh(FC_simul(Isubdiag)),'rows','complete');
    fitt(iwe)=cc(2);

    %%%%% Meta & FCD
    BOLD=xs';
    Phase_BOLD=zeros(N,nn);
    signal_filt=zeros(N,nn);
    for seed=1:N
        BOLD(seed,:)=demean(detrend(BOLD(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt2,afilt2,BOLD(seed,:));
        Phase_BOLD(seed,:) = angle(hilbert(signal_filt(seed,:)));
    end
    T=10:Tmax-10;
    for t=T
        ku=sum(complex(cos(Phase_BOLD(:,t)),sin(Phase_BOLD(:,t))))/N;
        sync(t-9)=abs(ku);
        for i=1:N
            for j=1:i-1
                patt(i,j)=cos(adif(Phase_BOLD(i,t),Phase_BOLD(j,t)));
            end
        end
        pattern(t-9,:)=patt(Isubdiag);
    end
    metastability(iwe)=abs(metastabilitydata-std(sync));


    %%%% KL dist between PTR2emp

    [PTRsim,Pstates]=LEiDA_fix_cluster(xs',NumClusters,Vemp,TR);


    %% PMS fitting
    klpstatesawake(iwe)=0.5*(sum(Pstates.*log(Pstates./P1emp))+sum(P1emp.*log(P1emp./Pstates)));
    klpstatesUWS(iwe)=0.5*(sum(Pstates.*log(Pstates./P3emp))+sum(P3emp.*log(P3emp./Pstates)));

    %% extra fitting
    kldistawake(iwe)=KLdist(PTR1emp,PTRsim);
    kldistUWS(iwe)=KLdist(PTR3emp,PTRsim);
    entropydistawake(iwe)=EntropyMarkov(PTR1emp,PTRsim);
    entropydistUWS(iwe)=EntropyMarkov(PTR3emp,PTRsim);

    PTRsimul(iwe,:,:)=PTRsim;

    Pstatessimul(iwe,:)=Pstates;

    iwe=iwe+1

    save Modelling_Workspace_UWS.mat


end
save Modelling_Workspace_UWS_all.mat metastability metastabilitydata WE PTRsimul Pstatessimul  klpstatesawake klpstatesUWS kldistUWS kldistawake entropydistawake entropydistUWS fitt Coptim n_Subjects f_diff;


