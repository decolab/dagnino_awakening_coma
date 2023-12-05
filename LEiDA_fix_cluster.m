function [PTRANSITION,Pstates] = LEiDA_fix_cluster(BOLD,NumClusters,Center,TR)

[N_areas, Tmax]=size(BOLD);

% Preallocate variables to save FC patterns and associated information
%Leading_Eig=zeros(Tmax,2*N_areas); % All leading eigenvectors %PD
Leading_Eig=zeros(Tmax,N_areas); % All leading eigenvectors %PD

% Bandpass filter settings
fnq=1/(2*TR);                 % Nyquist frequency
flp = 0.04;                    % lowpass frequency of filter (Hz)
fhi = 0.07;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
clear fnq flp fhi Wn k

t_all=0; % Index of time (starts at 0 and will be updated until n_Sub*Tmax)

Phase_BOLD=zeros(N_areas,Tmax);

% Get the BOLD phase using the Hilbert transform
for seed=1:N_areas
    BOLD(seed,:)=BOLD(seed,:)-mean(BOLD(seed,:));
    signal_filt =filtfilt(bfilt,afilt,BOLD(seed,:));
    Phase_BOLD(seed,:) = angle(hilbert(signal_filt));
end

for t=1:Tmax
    
    %Calculate the Instantaneous FC (BOLD Phase Synchrony)
    iFC=zeros(N_areas);
    for n=1:N_areas
        for p=1:N_areas
            iFC(n,p)=cos(Phase_BOLD(n,t)-Phase_BOLD(p,t));
        end
    end
    
    % Get the leading eigenvector
    [V1 ~]=eigs(iFC,1);

    
    % Make sure the largest component is negative
    if mean(V1>0)>.5
        V1=-V1;
    elseif mean(V1>0)==.5 && sum(V1(V1>0))>-sum(V1(V1<0))
        V1=-V1;
    end
    
    
    % Save V1 from all frames in all fMRI sessions in Leading eig
    t_all=t_all+1; % Update time
    Leading_Eig(t_all,:)=V1;
end

clear signal_filt iFC V1 Phase_BOLD

%% 2 - Cluster the Leading Eigenvectors
IDX=zeros(t_all,1);

for t=1:t_all
    for j=1:NumClusters
     di(j)=sqrt(sum((Leading_Eig(t,:)-Center(j,:)).^2));
    end
    [aux indmin]=min(di);
    IDX(t)=indmin;
end

Pstates=zeros(1,NumClusters);
for c=1:NumClusters
    Pstates(c)=mean(IDX==c);
end

Pstates=Pstates/sum(Pstates);

PTRANSITION=zeros(NumClusters,NumClusters);
i=1;
for c1=1:NumClusters
    j=1;
    for c2=1:NumClusters
        sumatr=0;
        for t=1:length(IDX)-1
            if IDX(t)==c1 && IDX(t+1)==c2
                sumatr=sumatr+1;
            end
        end
        if length(find(IDX(1:length(IDX)-1)==c1)) ~= 0
            PTRANSITION(i,j)=sumatr/(length(find(IDX(1:length(IDX)-1)==c1)));
        end
        j=j+1;
    end
    i=i+1;
end

 