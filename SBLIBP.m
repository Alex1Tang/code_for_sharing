function res = SBLIBP(paras)
%% IBPSBL
% Code of this paper:
% T. Tang, C. Yang, S. Yan, L. Xu and D. Chen, "A Sparse Bayesian Learning-Based Approach With Indian Buffet 
% Process Prior for Joint Wideband DOA and Frequency Band Estimation", IEEE Trans. Aerosp. Electron. Syst. .
% 
% Copyright 2025 by Tao Tang
% The code is sorted out at the University of Chinese Academy of Sciences
% Date: 01-Sept-2025
% Matlab R2022b
%%

Y = paras.Y;
A = paras.A;
[M, T, fft_num] = size(Y);
N = size(A, 2);

alpha0 = 1 / paras.sigma2;
rho = paras.rho / T;
alpha = paras.alpha;

maxiter = paras.maxiter;
tol = paras.tolerance;               % stopping criterion
converged = false;

% hyper-parameter of Gamma distribution
a = 1e-4; b = 1e-4;
c = 1e-4; d = 1e-4;
e = 1e-4; h = 1e-4;
     
featureNum = 30;    % featureNum of IBP

alpha_prev = alpha.';
alpha0_update = alpha0;

for ff = 1:fft_num
    sigmaMat_update(:,:,ff) = diag(ones(N,1));
    inv_sigmaMat_update(:,:,ff) = inv(sigmaMat_update(:,:,ff));
end
x_update = zeros(N,T,fft_num);

alpha1 = 0.1; 

tau = rand(featureNum,2);
Phi_update = ones(N,N,featureNum);

Zema_update = rand(fft_num, featureNum);    
Zema_update = Zema_update./repmat(sum(Zema_update,2),1,featureNum);    % IBP过程的Z的伯努利分布的参数矩阵

phi_update = rand(featureNum,N);

% OGSBI initialization
for ff = 1:fft_num
    iter = 0;
    while iter < 100
        iter = iter + 1; Phi = A(:,:,ff);
        C = 1 / alpha0 * eye(M) + Phi * diag(alpha_prev(ff,:)) * Phi'; Cinv = inv(C);
        Sigma(:,:,ff) = diag(alpha_prev(ff,:)) - diag(alpha_prev(ff,:)) * Phi' * Cinv * Phi * diag(alpha_prev(ff,:));
        mu(:,:,ff) = alpha0 * Sigma(:,:,ff) * Phi' * Y(:,:,ff);
        % update alpha
        musq = mean(abs(mu(:,:,ff)).^2, 2); alpha_prev(ff,:) = musq + real(diag(Sigma(:,:,ff)));
        if rho ~= 0
            alpha_prev(ff,:) = -.5 / rho + sqrt(.25 / rho^2 + alpha_prev(ff,:) / rho);
        end
        resid = Y(:,:,ff) - Phi * mu(:,:,ff); gamma1 = 1 - real(diag(Sigma(:,:,ff))) ./ (alpha_prev(ff,:).' + 1e-10);
        alpha0 = (T * M + a - 1) / (norm(resid, 'fro')^2 + T / alpha0 * sum(gamma1) + b);
    end
end
%% start IBP-SBL
iter = 0;
while ~converged
    iter = iter + 1;
    phi_last = phi_update;
    % update of xf
    for ff = 1:fft_num
        sigmaMat_update(:,:,ff) = alpha0_update*(A(:,:,ff)')*A(:,:,ff) + diag(1./alpha_prev(ff,:));
        inv_sigmaMat_update(:,:,ff) = inv(sigmaMat_update(:,:,ff));
        x_update(:,:,ff) = alpha0_update*inv_sigmaMat_update(:,:,ff)*A(:,:,ff)'*Y(:,:,ff);
    end

    % undate of alpha0
    a_update = (M*T*fft_num+a);
    temp = 0;
    for ff=1:fft_num
        temp = temp + norm((Y(:,:,ff)-A(:,:,ff)*x_update(:,:,ff)),'fro').^2 + T * trace(A(:,:,ff)'*A(:,:,ff)*inv_sigmaMat_update(:,:,ff));
    end
    b_update = real(temp+b);
    alpha0_update = a_update/b_update;

    % update of alphaG
    c_update = 0.5 + c;
    for k = 1:featureNum
        d_update(k,:) = d + 0.5*phi_update(k,:).^2 + 0.5 * diag(Phi_update(:,:,k)).';
    end
    alphaG_update = c_update./d_update;

    % update of alphaA
    e_update = e + N*fft_num/2;

    temp = 0;
    for ff = 1:fft_num
        temp1 = 0;
        for k = 1:featureNum
            temp1 = temp1 + Zema_update(ff,k)*(phi_update(k,:)*phi_update(k,:).' + trace(Phi_update(:,:,k)) ...
                + phi_update(k,:)*(Zema_update(ff,:)*phi_update-Zema_update(ff,k).*phi_update(k,:)).');
        end
        temp = temp + alpha_prev(ff,:)*alpha_prev(ff,:).' - 2*alpha_prev(ff,:)*(Zema_update(ff,:)*phi_update).' + temp1;
    end

    h_update = h + 0.5*temp;
    alphaA_update = e_update/h_update; 

    % update of Z
    for k = 1:featureNum
        theta = (psi(tau(k,1))-psi(tau(k,2))- 0.5*alphaA_update*(trace(Phi_update(:,:,k))+phi_update(k,:)*phi_update(k,:).') +...
            alphaA_update*phi_update(k,:)*(alpha_prev-Zema_update*phi_update+Zema_update(:,k)*phi_update(k,:)).');
        Zema_update(:,k) = 1./(1+exp(-theta));
    end

    % update of Gk
    for k = 1:featureNum
        Phi_update(:,:,k) = diag(1./(alphaG_update(k,:) + alphaA_update * sum(Zema_update(:,k))));
        phi_update(k,:) = alphaA_update * Zema_update(:,k)'*(alpha_prev - Zema_update * phi_update + Zema_update(:,k) * phi_update(k,:)) * Phi_update(:,:,k);
        phi_update(k, find(phi_update(k,:) < 0)) = 1e-5;
        if iter<=20
        [~,idx] = max(phi_update(k,:));
        phi_update(k,[1:idx-1, idx+1:end]) = 0;
        end
    end
    
    % update of pi
    Zema_sum = sum(Zema_update)';
    tau(:,1) = alpha1/featureNum + Zema_sum;
    tau(:,2) = fft_num + 1 - Zema_sum;

    % update of Lambda 
    for ff = 1:fft_num
        for n = 1:N
            % solve 
            Stemp=zeros(1,4);

            Stemp(1) = alphaA_update; Stemp(2) = -alphaA_update.*Zema_update(ff,:) * phi_update(:,n); 
            Stemp(3) = T; Stemp(4) = - (x_update(n,:,ff)*(x_update(n,:,ff)') + T*inv_sigmaMat_update(n,n,ff));
            
            Stemp = real(Stemp);
            S_solve = roots(Stemp);
            % choose the positive real root
            S_solveReal = S_solve(~logical(imag(S_solve)));
            idxTemp = S_solveReal > 0;
            S_solveReal = S_solveReal(idxTemp);
            % check
%             check = Stemp(1)*S_solveReal.^3 + Stemp(2)*S_solveReal.^2 + Stemp(3)*S_solveReal + Stemp(4);
            if length(S_solveReal) == 1
                alpha_current(ff,n) = S_solveReal;
            else
                % choose the root maximum variational posterior
                logPostPro = -T*log(S_solveReal) + Stemp(4)./S_solveReal - 0.5*alphaA_update.*(S_solveReal).^2 - Stemp(2).*S_solveReal;
                [~,idxTemp] = max(logPostPro);
                alpha_current(ff,n) = S_solveReal(idxTemp);
            end
        end
    end

    % stopping criteria
    if norm(phi_update - phi_last,'fro')/norm(phi_last,'fro') < tol || iter >= maxiter
        converged = true;
    end 
    alpha_prev = real(alpha_current);
end
res.alpha = alpha_current;
res.iter = iter;
res.Zema = Zema_update;
res.sigma2 = 1/alpha0_update;
res.phi = phi_update;
res.x = x_update;
res.featureNum = featureNum;

end
