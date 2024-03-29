
function [Phi,dPhi,Ind] = adaex_sk2_only_FLIMFRET_fixed_tauD_pp_69_70(alpha)

    global t; % picoseconds

    global DT;
        
    global IRF;  
    
    global Tp;
    
    global y;
    
    global E; % efficiency
    
    global E_avr;
    global tau_FRET_avr;
    
    global fixed_tauD;
        
    N = numel(t);

    Ndecays = length(y)/N;

    tauD = fixed_tauD;
    RDA_div_R0 = alpha(1);

    Nb = numel(t);
    Np = numel(E);

    tau = tauD.*(1-E);

    pE = efficiency_distribution(RDA_div_R0,E);

    FRETSUM = zeros(1,Nb);
    pE_deficit = 0;
    for p = 1:Np        
        conv_irf = conv_irf_pp_69_70(tau(p));
        if 0~=sum(isnan(conv_irf))
            %disp([p tau(p)]);
            pE_deficit = pE_deficit + pE(p);
            continue;
        end
        FRETSUM = FRETSUM + pE(p)/DT*conv_irf;
    end  
    %
%     if 0~=pE_deficit
%         disp(pE_deficit); % don't happen at all? 
%     end
    %
    A1 = FRETSUM/(sum(pE)-pE_deficit); % ??
    
% standard convolution method            
%     for p = 1:Np
%         FRETSUM = FRETSUM + pE(p)*(1-exp(-DT/tau(p)))*exp(-t/tau(p))/(1-exp(-Tp/tau(p)));
%         %FRETSUM = FRETSUM + pE(p)/tau(p)*exp(-t/tau(p))/(1-exp(-Tp/tau(p)));
%     end
%     
%     A1 = FRETSUM/sum(pE);                
%     
%     A1 = conv([ A1'; A1' ],IRF)';
%     A1 = A1(N+1:2*N);

    Phi = zeros(N*Ndecays,Ndecays);
        
    for k=1:Ndecays
        Phi((k-1)*N+1:k*N,k) = A1; 
    end

    dPhi = [];
    Ind = [];
    %
    E_avr = sum(pE.*E)/sum(pE);
    tau_FRET_avr = sum(tau.*pE)/sum(pE);
    
    Phi = sparse(Phi);
    
end
