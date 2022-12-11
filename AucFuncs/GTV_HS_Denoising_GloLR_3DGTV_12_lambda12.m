function [MH, MS,BH] = GTV_HS_Denoising_GloLR_3DGTV_12_lambda12(NH,Direc,para)

[mm,nn,bb] = size(NH);
block_sz = para.block_sz;

beta = para.beta;
lambda1 = para.lambda;
lambda2 = para.lambda2;

R=Direc.R;
A=Direc.A;
sha=Direc.sha;
clear Direc;

[ind_edge,ind_smooth]=find_block(A,sha);
len_edge=length(ind_edge);
len_smooth=length(ind_smooth);

X = reshape(NH,[mm*nn,bb]);
S = zeros(mm*nn,bb);

Y_f = Extract_blocks(NH,para);

Z = X;
norm_two = lansvd(Z, 1, 'L');
norm_inf = norm( Z(:), inf) / lambda1;
dual_norm = max(norm_two, norm_inf);
Z = Z / dual_norm;
Z_G = reshape(Z,[mm,nn,bb]);
Z_f = Extract_blocks(Z_G,para);

Z2 = Z;

P1_f = Z_f(:,:,:,ind_edge);
P2_f = Z_f(:,:,:,ind_edge);
P3_f = Z_f(:,:,:,ind_edge);
[Hx,Hy,Hz]=gradient(NH);
Mx_f = Extract_blocks(Hx,para);
My_f = Extract_blocks(Hy,para);
Mz_f = Extract_blocks(Hz,para);
M1_f = Mx_f(:,:,:,ind_edge);
M2_f = My_f(:,:,:,ind_edge);
M3_f = Mz_f(:,:,:,ind_edge);

mu1 = 6/norm_two;
mu2 = 6/norm_two;
mu3 = 6/norm_two;
rho=1.5;
mu_bar = 10e15;
max_Iter = 50;
tol = 10e-10;
converged = false;
iter = 0;
sv = 10;

while ~converged
    iter = iter + 1;
    
    % update Y
    A2 = X - S + Z2/mu2;
    A1_f = Y_f + Z_f/mu1;
    A1_f_G = Aggre_blocks(A1_f,para,size(NH));
    A1 = reshape(A1_f_G,[mm*nn,bb]);
    Y_aux = (A1*mu1+mu2*A2)/(mu2+mu1);
    [Y, sv] = SVShrinkage(Y_aux, 1/(mu1+mu2), sv);
    
    % update S
    S_aux = X - Y + Z2/mu2;
    S = Shrinkage(S_aux, beta/mu2);
    
    Y_G = reshape(Y,[mm,nn,bb]);
    Y_G_all = Extract_blocks(Y_G,para);
    % update Yi_f
    % update Yi_f: smooth blocks
    for ii = 1:len_smooth
        i_smooth = ind_smooth(ii);
        Y_G_all_smooth = Y_G_all(:,:,:,i_smooth);
        Z_f_smooth = Z_f(:,:,:,i_smooth);
        
        r = R(:,:,i_smooth);
        a = A(:,:,i_smooth);
        a = diag(a);
        
        [f1,f2,f3] = filter_direc3DTV(r,a);
        
        f1_FFT=fftn(f1,block_sz);
        f2_FFT=fftn(f2,block_sz);
        f3_FFT=fftn(f3,block_sz);
        
        A_smooth = mu1 + lambda1*conj(f1_FFT).*f1_FFT + lambda1*conj(f2_FFT).*f2_FFT...
            + lambda1*conj(f3_FFT).*f3_FFT;
        T_smooth_FFT = fftn(Y_G_all_smooth - Z_f_smooth/mu1);
        Y_f(:,:,:,i_smooth) = ifftn(mu1*T_smooth_FFT./A_smooth,'symmetric');
    end
    
    % update Yi_f: edge blocks
    for jj=1:len_edge
        i_edge = ind_edge(jj);
        Y_G_all_edge = Y_G_all(:,:,:,i_edge);
        Z_f_edge = Z_f(:,:,:,i_edge);
        
        r = R(:,:,i_edge);
        a = A(:,:,i_edge);
        a = diag(a);
        
        [f1,f2,f3] = filter_direc3DTV(r,a);
        
        f1_FFT=fftn(f1,block_sz);
        f2_FFT=fftn(f2,block_sz);
        f3_FFT=fftn(f3,block_sz);
        
        % update Yi_f
        R1_edge = M1_f(:,:,:,jj) + P1_f(:,:,:,jj)/mu3;
        R2_edge = M2_f(:,:,:,jj) + P2_f(:,:,:,jj)/mu3;
        R3_edge = M3_f(:,:,:,jj) + P3_f(:,:,:,jj)/mu3;
        DTR = mu3*(conj(f1_FFT).*fftn(R1_edge) + conj(f2_FFT).*fftn(R2_edge) + conj(f3_FFT).*fftn(R3_edge));
        T_edge_FFT = fftn(Y_G_all_edge - Z_f_edge/mu1);
        A_edge = mu1 + mu3*conj(f1_FFT).*f1_FFT + mu3*conj(f2_FFT).*f2_FFT...
            + mu3*conj(f3_FFT).*f3_FFT;
        Y_f(:,:,:,i_edge) = ifftn((mu1*T_edge_FFT+DTR)./A_edge,'symmetric');
        
        % update Mi_f
        RconvY1= ifftn(f1_FFT.*fftn(Y_f(:,:,:,i_edge)),'symmetric');
        RconvY2= ifftn(f2_FFT.*fftn(Y_f(:,:,:,i_edge)),'symmetric');
        RconvY3= ifftn(f3_FFT.*fftn(Y_f(:,:,:,i_edge)),'symmetric');
        MT1_edge = RconvY1 - P1_f(:,:,:,jj)/mu3;
        MT2_edge = RconvY2 - P2_f(:,:,:,jj)/mu3;
        MT3_edge = RconvY3 - P3_f(:,:,:,jj)/mu3;
        
        M1_f(:,:,:,jj) = Shrinkage(MT1_edge, lambda2/mu3);
        M2_f(:,:,:,jj) = Shrinkage(MT2_edge, lambda2/mu3);
        M3_f(:,:,:,jj) = Shrinkage(MT3_edge, lambda2/mu3);
        
        % update multiplier for edge blocks
        P1_f(:,:,:,jj) = P1_f(:,:,:,jj) + mu3*(M1_f(:,:,:,jj) - RconvY1);
        P2_f(:,:,:,jj) = P2_f(:,:,:,jj) + mu3*(M2_f(:,:,:,jj) - RconvY2);
        P3_f(:,:,:,jj) = P3_f(:,:,:,jj) + mu3*(M3_f(:,:,:,jj) - RconvY3);
    end
       
    % update multipliers
    Z2 = Z2 + mu2*(X - Y - S);
    Z_f = Z_f + mu1*(Y_f - Y_G_all);
    
    % update mu
    mu1 = min(rho*mu1, mu_bar);
    mu2 = min(rho*mu2, mu_bar);
    mu3 = min(rho*mu3, mu_bar);
    
    stopCriterion = norm(X - Y - S, 'fro') / norm(X, 'fro');
    if stopCriterion < tol
        converged = true;
    end
    
    if ~converged && iter >= max_Iter
        converged = 1 ;
    end
%     disp(['********Iter:' num2str(iter) '********Error:' num2str(stopCriterion) ...
%         '********Converged:' num2str(converged)]);
end

BH = Aggre_blocks(Y_f,para,size(NH));

MH = reshape(Y,[mm,nn,bb]);
MS = reshape(S,[mm,nn,bb]);





