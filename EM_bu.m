function [Q,cr,bet1,alp1,rho1,itr] = EM_bu(A)
% [Q,cr,bet1,alp1,rho1,itr] = EM_bu(A)
% Find the *true* network using Mark Newman's ER model 
% assumes the input is binary and directed networks.
% The input should be a 3D matrix, subject's adj matrix stacked at the
% top of each others. 
%
%%%REFERENCE
% Newman, M. E. J. "Network structure from rich but noisy data." 
% Nature Physics 14.6 (2018): 542.
%_________________________________________________________________________
% Soroosh Afyouni, University of Oxford, 2018
% srafyouni@gmail.com
fnnf=mfilename; if ~nargin; help(fnnf); return; end; clear fnnf;
%_________________________________________________________________________

    N = size(A,1);
    M = size(A,3);
    
    crfnc = @(E1,E0) abs((E1-E0)./E0);
    tu    = @(G) triu(G,1);

    E = sum(A,3);

    %initialisation
    alp1 = 0.8;
    bet1 = 0.01;
    rho1 = 0.5;

    itr = 0;
    cr_tmp = ones(1,3).*0.1;
    epsilon = 10e-6;

    while true

        itr = itr + 1; 

        alp0 = alp1; bet0 = bet1; rho0 = rho1; 

        for i = 1:N
            for j = 1:N
                nom    = rho0.*alp0.^E(i,j).*(1-alp0).^(M-E(i,j));
                dnom   = nom + ( (1-rho0).*bet0^E(i,j).*(1-bet0).^(M-E(i,j)));
                Q(i,j) = nom./dnom;
            end
        end

        rho1 = (1./nchoosek(N,2)) .* sum(sum(tu(Q)));
        bet1 = sum(sum(tu(E.*(1-Q))))./sum(sum(tu(M.*(1-Q))));
        alp1 = sum(sum(tu(E.*Q)))./sum(sum(tu(M.*Q)));

        cr_tmp = [crfnc(alp0,alp1) crfnc(bet0,bet1) crfnc(rho0,rho1)]; 
        cr(:,itr) = cr_tmp;

        if sum(cr_tmp<epsilon)==3 || sum(isnan(cr_tmp<epsilon))==3

            disp(['================================='])
            disp(['Converged on: ' num2str(itr) ' iteration.'])
            disp(['FPR: ' num2str(bet1)])
            disp(['TPR: ' num2str(alp1)])
            disp(['RHO: ' num2str(rho1)])
            disp(['================================='])
            break; 
        end;
    end

end


