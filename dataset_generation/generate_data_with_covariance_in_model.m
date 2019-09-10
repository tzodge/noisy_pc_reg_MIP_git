clear;
clc;

close all
warning('off');
shapename='bunny';%name of the stl file
d2r = pi/180;
noise= 1;% If you have noise, turn this to 1
Ns=1000;
misalignment = 20;
shapename_save = strcat(shapename,string(Ns),'noisy','_misalignment_',string(misalignment));
Ns_sampled = 20;


%% 
%%% to generate partially overlapping datas, set both sensor and model
%%% switch as 1. run code for view = 1 then without changing the
%%% shapename_save run the same code for view = 2 .
partial_sensor_switch = 0;
partial_model_switch = 0;
%%% view == 1 --> cam_position= [1,0,0]; 
%%% view == 2 --> cam_position= [1,1,0]; 
view = 1; %% can take value 1 and 2. only relevant for partial data.  
theta = 30; %deg
if view == 1
    cam_position= [0,1,0]; 
elseif view == 2
    cam_position= [1*cos(theta*d2r),1*theta*d2r,0]; 
end
%%% keep view = 1 if not partial
%%


% I sampled 10 points from each face to form Msampled
pointPerFace=1;
plot_switch = 0;
[tF,tP]=stlread(strcat("STL_files/",shapename,'.stl'));

% next two lines to scale and set to origin
tP=tP/(max(max(tP)-min(tP)))*1;
tP=tP-mean(tP);

M=unique(tP,'rows');

gg=M;
 
% vertices are repeated in stl file. So I need to select only unique ones.
tF1=tF;
tF_partial=tF;
for ii=1:length(tF)
    for jj=1:3
        tF1(ii,jj)=find(gg(:,1)==tP(tF(ii,jj),1) & gg(:,2)==tP(tF(ii,jj),2) & gg(:,3)==tP(tF(ii,jj),3));
        tF_partial(ii,jj)=find(gg(:,1)==tP(tF(ii,jj),1) & gg(:,2)==tP(tF(ii,jj),2) & gg(:,3)==tP(tF(ii,jj),3));
    end
end

 


Nf=length(tF1);
Nm=length(M);

%%

F=zeros(Nf,Nm);

for ii=1:Nf
    F(ii,tF1(ii,:))=1;
end

%%
noiseSigmaM_z= 1e-3 ;% also try 1e-1 for more noise and 1e-2 for less noise
noiseSigmaM_xy = 1; %%%%%%%%%%%% last values tested in log_1 --> z = 1e-3, xy=1

noiseSigmaS=2e-2;  %%% default value is 1e-2


Nmsampled=Nf*pointPerFace;
Msampled=zeros(Nmsampled,3);
tt=1;


        
%%
% d2r=pi/180;
% Xgt=[0.5*(2*rand(1,3)-1),180*d2r*(2*rand(1,3)-1)];
% Xgt(4:6)=rotm2eul(eul2rotm(Xgt(4:6)));
% 
% gt=eye(4);
% Rgt=eul2rotm(Xgt(4:6));
% tgt=Xgt(1:3)';
% gt(1:3,1:3)=Rgt;
% gt(1:3,4)=tgt;

%%% temporary change to generate gt of varying misalignemnt
gt=eye(4);
Rgt_axang = [1,0,0,misalignment*d2r];
Rgt = axang2rotm(Rgt_axang);
tgt = 2*rand(3,1)-1;
gt(1:3,1:3) = Rgt;
gt(1:3,4)=tgt;
%%% temporary change to generate gt of varying misalignemnt

 
invTgt=inv(gt);


%%

outliers=0;%IF you have outliers, turn this to 1
outlierPercentage=60;

partial=0;
partialPercentage=60;% also try 80 and 30

if partial==1
    Ns=floor(Ns*100/40);
end

o=zeros(Ns,1);
S=zeros(Ns,3);
S_noise_free = zeros(Ns,3);
C=zeros(Ns,Nm);
Cb=zeros(Ns,Nmsampled);
f1=zeros(Ns,Nf);
SigmaS=zeros(Ns,6);
B=zeros(Ns,6);
t=1;
% wb=waitbar(0);
for ii=1:ceil(Ns/Nf)*Nf
    
%     waitbar(ii/(ceil(Ns/Nf)*Nf),wb);
    idx=randi(Nf);
    w=rand(1,3);
    w=w/sum(w);
    C(t,tF1(idx,:))=w;
    f1(t,idx)=1;
    
    S(t,:)= w(1)*M(tF1(idx,1),:)+w(2)*M(tF1(idx,2),:)+w(3)*M(tF1(idx,3),:);
        idx1=knnsearch(Msampled((idx-1)*pointPerFace+1:idx*pointPerFace,:),S(t,:));
    Cb(t,(idx-1)*pointPerFace+idx1)=1;
    S_noise_free(t,:) = S(t,:);
    
    CovS=0*eye(3);
    SigmaS(t,:)=[CovS(1,1:3),CovS(2,2:3),CovS(3,3)];
    B(t,:)=SigmaS(t,:);
    % if the variable outliers in 1 then I add outliers. 1% of points
    if outliers==1 && rand>1-(outlierPercentage/100)
        S(t,:)=S(t,:)+(2*rand(1,3)-1)*0.5;
        C(t,:)=0;
        o(t,:)=1;
        f1(t,:)=0;
        
        Cb(t,:)=0;
       
    end
    Rn=eye(3);
    if noise==1
        nvec=cross((tP(tF(idx,2),:)-tP(tF(idx,1),:)),(tP(tF(idx,3),:)-tP(tF(idx,1),:)));
        nvec=nvec/norm(nvec);
        tmp2=cross(rand(1,3),nvec);tmp2=tmp2/norm(tmp2);
        tmp1=cross(tmp2,nvec);
        Rn(:,1)=tmp1';
        Rn(:,2)=tmp2';
        Rn(:,3)=nvec';
        S(t,:)=S(t,:)+(Rn*[normrnd(0,noiseSigmaS),normrnd(0,noiseSigmaS),normrnd(0,noiseSigmaS)]')';
        CovS=invTgt(1:3,1:3)*Rn*diag([noiseSigmaS^2,noiseSigmaS^2,noiseSigmaS^2])*Rn'*invTgt(1:3,1:3)';
        cholInvCov=chol(inv(CovS));
        B(t,:)=[cholInvCov(1,1:3),cholInvCov(2,2:3),cholInvCov(3,3)];
        SigmaS(t,:)=[CovS(1,1:3),CovS(2,2:3),CovS(3,3)];
        
    end
    
   
    
    S(t,:)=transpose(invTgt(1:3,1:3)*S(t,:)'+invTgt(1:3,4));
    S_noise_free(t,:) = transpose(invTgt(1:3,1:3)*S_noise_free(t,:)'+invTgt(1:3,4));

    
    
    
    
    if t==Ns
        break;
    end
    
    t=t+1;
end

%%
for jj=1:Nf
        nvec=cross((tP(tF(jj,2),:)-tP(tF(jj,1),:)),(tP(tF(jj,3),:)-tP(tF(jj,1),:)));
        nvec=nvec/norm(nvec);
        tmp2=cross(rand(1,3),nvec);tmp2=tmp2/norm(tmp2);
        tmp1=cross(tmp2,nvec);
        Rn(:,1)=tmp1';
        Rn(:,2)=tmp2';
        Rn(:,3)=nvec';
        
        if noise==0
            CovM=Rn*diag([noiseSigmaM_xy,noiseSigmaM_xy,noiseSigmaM_z])*Rn';
        else
            CovM=Rn*diag([noiseSigmaM_xy,noiseSigmaM_xy,noiseSigmaM_z])*Rn'*((3*noiseSigmaS^2)/(2*noiseSigmaM_xy+noiseSigmaM_z));
        end
        cholInvCov=chol(inv(CovM));
        cholInvCov_M_and_S = chol(inv(CovM+CovS));
    
 
        if pointPerFace == 1
            for ii=1:pointPerFace
                Msampled(tt,:)= (M(tF1(jj,1),:) + M(tF1(jj,2),:)+ M(tF1(jj,3),:))/3;
                BMsampled(tt,:)=[cholInvCov(1,1:3),cholInvCov(2,2:3),cholInvCov(3,3)];
                B_combined(tt,:) = [cholInvCov_M_and_S(1,1:3),cholInvCov_M_and_S(2,2:3),cholInvCov_M_and_S(3,3)];
                SigmaMsampled(tt,:)=[CovM(1,1:3),CovM(2,2:3),CovM(3,3)];
                tt=tt+1;
            end
        else
            for ii=1:pointPerFace
                w=rand(1,3);
                w=w/sum(w);
                Msampled(tt,:)=w(1)*M(tF1(jj,1),:)+w(2)*M(tF1(jj,2),:)+w(3)*M(tF1(jj,3),:);
                BMsampled(tt,:)=[cholInvCov(1,1:3),cholInvCov(2,2:3),cholInvCov(3,3)];
                B_combined(tt,:) = [cholInvCov_M_and_S(1,1:3),cholInvCov_M_and_S(2,2:3),cholInvCov_M_and_S(3,3)];
                SigmaMsampled(tt,:)=[CovM(1,1:3),CovM(2,2:3),CovM(3,3)];
                tt=tt+1;
            end
        end
end
       

%%
if partial==1
    
    idx1=randi(length(S));
    idx2=knnsearch(S,S(idx1,:),'k',floor(partialPercentage/100*Ns));
    C=C(idx2,:);
    f1=f1(idx2,:);
    
    Cb=Cb(idx2,:);
    o=o(idx2,:);

    S=S(idx2,:);
        B=B(idx2,:);
    SigmaS=SigmaS(idx2,:);

    idx3=randperm(length(S));
    S=S(idx3,:);
    B=B(idx3,:);
    SigmaS=SigmaS(idx3,:);
    Ns=length(S);
end

for ii=1:Ns
St=transpose(gt(1:3,1:3)*S'+gt(1:3,4));
St_noise_free = transpose(gt(1:3,1:3)*S_noise_free'+gt(1:3,4));
phi=0;
    phi=phi+sum(abs([B(ii,1),B(ii,2),B(ii,3);B(ii,2),B(ii,4),B(ii,5);B(ii,3),B(ii,5),B(ii,6)]*(gt(1:3,1:3)'*M'*C(ii,:)'-gt(1:3,1:3)'*gt(1:3,4)-S(ii,:)')));
    
end
phi=phi/Ns;


kdOBJ=KDTreeSearcher(Msampled);
GTcorrespondenceIdx=knnsearch(kdOBJ,St); 
GTcorrespondenceIdx = GTcorrespondenceIdx-1;

%%
distTotal_gicp=0;
GICPcorrespondenceIdx=zeros(Ns_sampled,1);
for ii=1:Ns_sampled
     disp(ii);
     tmp_min=1e100;
     SigmaSi=[SigmaS(ii,1),SigmaS(ii,2),SigmaS(ii,3);...
            SigmaS(ii,2),SigmaS(ii,4),SigmaS(ii,5);...
            SigmaS(ii,3),SigmaS(ii,5),SigmaS(ii,6)];

         Q=(Rgt*SigmaSi*Rgt');
%          Q=eye(3);
         
     for jj=1:size(St_noise_free,1)
         tmp_dist=(Rgt*S(ii,:)'+tgt-St_noise_free(jj,:)')'*(Q\(Rgt*S(ii,:)'+tgt-St_noise_free(jj,:)'));
         if tmp_dist<tmp_min
             GICPcorrespondenceIdx(ii,:)=jj;
             tmp_min=tmp_dist;
         end
     end
     distTotal_gicp=distTotal_gicp+tmp_min;
     
end


%% generating partial data using HPR 

if partial_sensor_switch == 1
    indices_partial_St = HPR(St, cam_position, exp(0.3));
    indices_partial_M = HPR(M, cam_position, exp(0.3));
    indices_partial_Msampled = HPR(Msampled, cam_position, exp(0.3));
     
end




 %% This part of the coed is to generate correspondence with model points, not S. 

% kdOBJ=KDTreeSearcher(Msampled);
%  idx1=knnsearch(kdOBJ,St_noise_free(1:100,:)); 
%  Mcgt=Msampled(idx1,:);
% distTotal_gt=0;
%  for ii=1:100
%      
%      SigmaSi=[SigmaS(ii,1),SigmaS(ii,2),SigmaS(ii,3);...
%             SigmaS(ii,2),SigmaS(ii,4),SigmaS(ii,5);...
%             SigmaS(ii,3),SigmaS(ii,5),SigmaS(ii,6)];
% 
%          Q=Rgt*SigmaSi*Rgt';
%      distTotal_gt=distTotal_gt+(Rgt*S(ii,:)'+tgt-Mcgt(ii,:)')'*(Q\(Rgt*S(ii,:)'+tgt-Mcgt(ii,:)'));
%  end
% distTotal_gicp=0;
% Mcgicp=0*Mcgt;
% for ii=1:100
%      disp(ii);
%      tmp_min=1e100;
%      SigmaSi=[SigmaS(ii,1),SigmaS(ii,2),SigmaS(ii,3);...
%             SigmaS(ii,2),SigmaS(ii,4),SigmaS(ii,5);...
%             SigmaS(ii,3),SigmaS(ii,5),SigmaS(ii,6)];
% 
%          Q=(Rgt*SigmaSi*Rgt');
% %          Q=eye(3);
%          
%      for jj=1:size(Msampled,1)
%          tmp_dist=(Rgt*S(ii,:)'+tgt-Msampled(jj,:)')'*(Q\(Rgt*S(ii,:)'+tgt-Msampled(jj,:)'));
%          if tmp_dist<tmp_min
%              Mcgicp(ii,:)=Msampled(jj,:);
%              tmp_min=tmp_dist;
%          end
%      end
%      distTotal_gicp=distTotal_gicp+tmp_min;
%      
%  end
 

%  figure();
%     trisurf(tF1,M(:,1),M(:,2),M(:,3),'FaceAlpha',0.5); hold on; axis equal
%    shading interp
%  axis equal; 
%  hold on;
%  scatter3(Mcgt(1:10,1),Mcgt(1:10,2),Mcgt(1:10,3),'g','fill');
%  hold on
%  scatter3(Mcgicp(1:10,1),Mcgicp(1:10,2),Mcgicp(1:10,3)',"*r");
%  scatter3(St(1:10,1),St(1:10,2),St(1:10,3),'b','fill');
%     for i = 1:10
%           
%         line([St(i,1),Mcgicp(i,1)],[St(i,2),Mcgicp(i,2)],[St(i,3),Mcgicp(i,3)],'Color','red');
%     end
%   
%     for i = 1:10
%         line([St(i,1),St_noise_free(i,1)],[St(i,2),St_noise_free(i,2)],[St(i,3),St_noise_free(i,3)],'Color','black');
%     end
% 
0
%%
if partial_sensor_switch ==1 
   
    S = S(indices_partial_St,:); 
    S_noise_free= S_noise_free(indices_partial_St,:);
    St = St(indices_partial_St,:);
    St_noise_free = St_noise_free(indices_partial_St,:);
    SigmaS = SigmaS(indices_partial_St,:);
    B = B(indices_partial_St,:);
    GTcorrespondenceIdx = GTcorrespondenceIdx(indices_partial_St,:);

end

if partial_model_switch == 1
    M=M(indices_partial_M,:);
    Msampled = Msampled(indices_partial_Msampled,:);
    % %% added by arun for covariance in model.
    SigmaMsampled = SigmaMsampled(indices_partial_Msampled,:);
    BMsampled = BMsampled(indices_partial_Msampled,:);
    B_combined = B_combined(indices_partial_Msampled,:);

end

 %%
if plot_switch == 1
%  figure();
    % trisurf(tF,tP(:,1),tP(:,2),tP(:,3),'FaceAlpha',0.8); hold on; axis equal
    % scatter3(Msampled(:,1),Msampled(:,2),Msampled(:,3),'b','fill');
    trisurf(tF1,M(:,1),M(:,2),M(:,3),'FaceAlpha',0.8); hold on; axis equal
%     scatter3(St(:,1),St(:,2),St(:,3),'.m');
    hold on
    scatter3(S(:,1),S(:,2),S(:,3),'.b');
%     scatter3(M(indices_partial_M,1),M(indices_partial_M,2),M(indices_partial_M,3),'.r');    
    scatter3(Msampled(:,1),Msampled(:,2),Msampled(:,3),'.r');
    scatter3(St(1:20,1),St(1:20,2),St(1:20,3),'k',"fill");
    scatter3(St(:,1),St(:,2),St(:,3),'*g');
    
%     scatter3(St_noise_free(:,1),St_noise_free(:,2),St_noise_free(:,3)',"*g");
%     scatter3(S(:,1),S(:,2),S(:,3),'.b');
    scatter3(cam_position(1),cam_position(2),cam_position(3), '.r')
    
    %%%%%%%%%%% joining corresponding noisy and noisefree sensor points
%     
%     a = size(St);
%     for i = 1:a(1)
%         line([St(i,1),St_noise_free(i,1)],[St(i,2),St_noise_free(i,2)],[St(i,3),St_noise_free(i,3)]);
%     end

end

%% 
% % % filename="/home/biorobotics/Desktop/tejas/gurobiCodesPython/datasets";
filename = "/home/biorobotics/Desktop/tejas/gurobiCodesPython/datasets";

mkdir(strcat(filename,'/',shapename_save))

%%%%%%%%%% changes with sensor pts %%%%%%%%
dlmwrite(strcat(filename,'/',shapename_save,'/S.txt'),S);
dlmwrite(strcat(filename,'/',shapename_save,'/S_noise_free.txt'),S_noise_free);
dlmwrite(strcat(filename,'/',shapename_save,'/St.txt'),St);
dlmwrite(strcat(filename,'/',shapename_save,'/St_noise_free.txt'),St_noise_free);
dlmwrite(strcat(filename,'/',shapename_save,'/gt.txt'),gt);
dlmwrite(strcat(filename,'/',shapename_save,'/SigmaS.txt'),SigmaS);
dlmwrite(strcat(filename,'/',shapename_save,'/B.txt'),B);
dlmwrite(strcat(filename,'/',shapename_save,'/GICPcorrespondenceIdx.txt'),GICPcorrespondenceIdx);
dlmwrite(strcat(filename,'/',shapename_save,'/o.txt'),o);
dlmwrite(strcat(filename,'/',shapename_save,'/GTcorrespondenceIdx.txt'),GTcorrespondenceIdx);% added by arun for covariance in model.
dlmwrite(strcat(filename,'/',shapename_save,'/Msampled.txt'),Msampled);
dlmwrite(strcat(filename,'/',shapename_save,"/F.txt"),F);
dlmwrite(strcat(filename,'/',shapename_save,'/M.txt'),M);
dlmwrite(strcat(filename,'/',shapename_save,'/SigmaMsampled.txt'),SigmaMsampled); % added by arun for covariance in model.
dlmwrite(strcat(filename,'/',shapename_save,'/BMsampled.txt'),BMsampled);% added by arun for covariance in model
dlmwrite(strcat(filename,'/',shapename_save,'/B_combined.txt'),B_combined);% added by tejas for corrected cov values
%%%%%%%%%% changes with model pts %%%%%%%%   for partial data in two views,
%%%%%%%%%% comment out this section for second view.
if view== 1
    dlmwrite(strcat(filename,'/',shapename_save,'/Msampled.txt'),Msampled);
    dlmwrite(strcat(filename,'/',shapename_save,"/F.txt"),F);
    dlmwrite(strcat(filename,'/',shapename_save,'/M.txt'),M);
    dlmwrite(strcat(filename,'/',shapename_save,'/SigmaMsampled.txt'),SigmaMsampled); % added by arun for covariance in model.
    dlmwrite(strcat(filename,'/',shapename_save,'/BMsampled.txt'),BMsampled);% added by arun for covariance in model
    dlmwrite(strcat(filename,'/',shapename_save,'/B_combined.txt'),B_combined);% added by tejas for corrected cov values
end
    % 

f = fopen(strcat(filename,'/',shapename_save,'/parameters.txt'),'wt');
fprintf(f,"noiseSigmaM_z = %f \n",noiseSigmaM_z);
fprintf(f,"noiseSigmaM_xy = %f \n",noiseSigmaM_xy);
fprintf(f,"noiseSigmaS = %f \n",noiseSigmaS);
fclose(f);
a = rotm2axang(gt(1:3,1:3));
'misalignment'
a(4)/d2r
'done,saved in'
disp( filename)