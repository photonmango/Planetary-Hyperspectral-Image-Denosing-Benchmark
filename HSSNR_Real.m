
clear,clc
currentFolder = pwd;
addpath(genpath(currentFolder));

load img.mat
f=img;
[M,N,B] =   size(f);


tic;
iba=[1,1:B-1];
df=f-f(:,:,iba);
windowsize  = 7;
windowfilt = ones(1,windowsize)/windowsize;
Scale=4;
wavelet_type='db1';
W=cell(Scale+1,1);
for i=1:Scale
    W{i}=cell(3,1);
end
u=zeros(M,N,B);
ufinal=zeros(M,N,B);
for i=2:B
    W{Scale+1}=df(:,:,i);
    for j=1:Scale%wavelet decomposition
        [W{Scale+1},W{j}{1},W{j}{2},W{j}{3}]=dwt2(W{Scale+1},wavelet_type);
    end
    Nsig=median(abs(W{1}{3}(:)))/0.6745;
    %     Nsig=45;
    for j=Scale:-1:1
        for dir=1:3
            Wsig = conv2(windowfilt,windowfilt,(W{j}{dir}).^2,'same');
            %                Wsig=sum(sum(W{j}{dir}.^2))/numel(W{j}{dir});
            Ssig = sqrt(max(Wsig-Nsig.^2,eps));
            T=Nsig^2./Ssig;
            W{j}{dir}=soft(W{j}{dir},T);
        end
        W{Scale+1}=idwt2(W{Scale+1},W{j}{1},W{j}{2},W{j}{3},wavelet_type);
        if j>1
            W{Scale+1}=W{Scale+1}(1:size(W{j-1}{1},1),1:size(W{j-1}{1},2));
        end
    end
    
    u(:,:,i)=W{Scale+1}(1:M,1:N);
end
for i=1:M
    for j=1:N
        W{Scale+1}=u(i,j,:);
        for k=1:Scale
            [W{Scale+1},W{k}{1}]=dwt(W{Scale+1},wavelet_type);
        end
        Nsig=median(abs(W{1}{1}(:)))/0.6745;
        %         Nsig=45;
        for k=Scale:-1:1
            
            Wsig = conv((W{k}{1}).^2,windowfilt,'same');
            %              Wsig=sum(sum(W{k}{1}.^2))/numel(W{k}{dir});
            Ssig = sqrt(max(Wsig-Nsig.^2,eps));
            T=Nsig^2./Ssig;
            W{k}{1}=soft(W{k}{1},T);
            
            W{Scale+1}=idwt(W{Scale+1},W{k}{1},wavelet_type);
            if k>1
                W{Scale+1}=W{Scale+1}(1:size(W{k-1}{1},1),1:size(W{k-1}{1},2));
            end
        end
        u(i,j,:)=W{Scale+1}(1:B);
    end
end

W{Scale+1}=f(:,:,1);
for j=1:Scale%wavelet decomposition
    [W{Scale+1},W{j}{1},W{j}{2},W{j}{3}]=dwt2(W{Scale+1},wavelet_type);
end
Nsig=median(abs(W{1}{3}(:)))/0.6745;
% Nsig=45;
for j=Scale:-1:1
    for dir=1:3
        Wsig = conv2(windowfilt,windowfilt,(W{j}{dir}).^2,'same');
        %                Wsig=sum(sum(W{j}{dir}.^2))/numel(W{j}{dir});
        Ssig = sqrt(max(Wsig-Nsig.^2,eps));
        T=Nsig^2./Ssig;
        W{j}{dir}=soft(W{j}{dir},T);
    end
    W{Scale+1}=idwt2(W{Scale+1},W{j}{1},W{j}{2},W{j}{3},wavelet_type);
    if j>1
        W{Scale+1}=W{Scale+1}(1:size(W{j-1}{1},1),1:size(W{j-1}{1},2));
    end
end
ufinal(:,:,1)=W{Scale+1}(1:M,1:N);

% ufinal(:,:,1)=f(:,:,1);
% reshape
for i=2:B
    ufinal(:,:,i)=ufinal(:,:,i-1)+u(:,:,i);
end

windowsize  = 7;
windowfilt = ones(1,windowsize)/windowsize;
% ufinal2=zeros(M*N,B);
f=reshape(f,M*N,B);
ufinal=reshape(ufinal,M*N,B);
% for i=1:M*N
%             ufinal2(i,:)=conv(f(i,:),windowfilt,'same')-conv((ufinal(i,:)),windowfilt,'same');
% end
ufinal2=imfilter(f,windowfilt,'symmetric')-imfilter(ufinal,windowfilt,'symmetric');
ufinal2=ufinal+ufinal2;
ufinal2=reshape(ufinal2,M,N,B);


HSSNR=ufinal2;
save('HSSNR.mat','HSSNR');


toc;



