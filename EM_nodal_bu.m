function [Q,bet1,alp1,rho1,FDR,itr] = EM_nodal_bu(A,alp1,bet1,rho1)
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
    idx = find(triu(ones(N),1));
    
    crfnc = @(E1,E0) abs((E1-E0)./E0);
    tu    = @(G) triu(G,1);

    E = sum(A,3);

if nargin==1    
    disp(['No input for initialisation values... So they are set by script!'])
    %initialisation
    alp1 = ones(1,N).*0.4;
    bet1 = ones(1,N).*0.3;
    rho1 = 0.5;
end

    itr = 0;
    cr_tmp = ones(1,3).*0.1;
    epsilon = 10e-6;

    while true

        itr = itr + 1; 

        alp0 = alp1; bet0 = bet1; rho0 = rho1; 

        for i = 1:N
            alp0i = alp0(i); 
            bet0i = bet0(i);
            for j = 1:N
                alp0j = alp0(j); 
                bet0j = bet0(j);
                nom    = rho0.*alp0i.^E(i,j).*(1-alp0i).^(M-E(i,j)) .* alp0j.^E(i,j).*(1-alp0j).^(M-E(i,j));                
                dnom   = nom + ( (1-rho0).*bet0i^E(i,j).*(1-bet0i).^(M-E(i,j)) .* bet0j^E(i,j).*(1-bet0j).^(M-E(i,j)) );
                Q(i,j) = nom./dnom;
            end
        end

        rho1 = (1./nchoosek(N,2)) .* sum(Q(idx));
        bet1 = sum(E.*(1-Q))      ./ sum(M.*(1-Q));
        alp1 = sum(E.*Q)          ./ sum(M.*Q);

        %cr_tmp = [crfnc(alp0,alp1); crfnc(bet0,bet1); crfnc(rho0,rho1)]; 
        %cr(:,:,itr) = cr_tmp;

        Alp_Changes = crfnc(alp0,alp1); Alp_Changes(isnan(Alp_Changes)) = [];
        Bet_Changes = crfnc(bet0,bet1); Bet_Changes(isnan(Bet_Changes)) = [];
        
        if sum(Alp_Changes<epsilon)==numel(Alp_Changes) && sum(Bet_Changes<epsilon)==numel(Bet_Changes) && crfnc(rho0,rho1)<epsilon %|| sum(isnan(cr_tmp<epsilon))==3

            disp(['================================='])
            disp(['Converged on: ' num2str(itr) ' iteration.'])
            disp(['FPR: ' num2str(bet1)])
            disp(['TPR: ' num2str(alp1)])
            disp(['RHO: ' num2str(rho1)])
            disp(['================================='])
            
            FDR = ((1-rho1).*bet1)./(rho0.*alp1+(1-rho1).*bet1);
            
            break; 
        end;
    end

end