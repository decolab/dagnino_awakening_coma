function perturb_hopf_syncbynode_MCS(node)

load Modelling_Workspace_MCS.mat

PERTURB=0:0.01:0.2;  %% This is synchronisation protocol

we= 0.07;  %% Insert here G value from optimal MCS model fitting
ITER=10;
N=214;
a=-0.02*ones(N,2);
C=squeeze(Coptim(find(abs(WE-we)<0.0001),:,:));

TSmax=295;
NSUB=29; %n_Subjects;
TR=2;  % Repetition Time (seconds)
NumClusters=Number_Clusters;

%%%%%%%%%%%%%%%%%%

omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);

dt=0.1*TR/2;
Tmax=NSUB*TSmax;
sig=0.01;
dsig = sqrt(dt)*sig; % to avoid sqrt(dt) at each time step

%%%%%%%%%%%%
%% Optimize
%%
wC = we*C;
sumC = repmat(sum(wC,2),1,2); % for sum Cij*xj

node
iwe=1;
for perturb=PERTURB
    a=-0.02*ones(N,2);
    a(node,:)=a(node,:)+perturb;
    for iter=1:ITER
        iter
        xs=zeros(Tmax,N);
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

        [PTRsim,Pstates]=LEiDA_fix_cluster(xs',NumClusters,Vemp,TR);
        KLps_p(iter)=0.5*(sum(Pstates.*log(Pstates./P1emp))+sum(P1emp.*log(P1emp./Pstates)));  %% Target condition fitting
        iPTRsim(iter,:,:)=PTRsim; 
        iPstates(iter,:)=Pstates;
    end
    KLpstatesMCS_perturbed(node,iwe)=nanmean(KLps_p);
    Pstates_perturbed(node,iwe,:)=nanmean(iPstates);
    PTRsim_perturbed(node,iwe,:,:)=squeeze(nanmean(iPTRsim,1));
    iwe=iwe+1;
end
save(sprintf('MCS_KLpstates_perturbed_syncProt_bynode_%s',num2str(node)),'PTRsim_perturbed','KLpstatesMCS_perturbed','PTRsim','Pstates_perturbed','Pstates')

end
