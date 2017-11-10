
function [Phi,dPhi,Ind] = adaex_sk2_sk2_FLIMFRET_pp_69_70(alpha)

    global t; % picoseconds

    global DT;
        
    global IRF;  
    
    global Tp;
    
    global y;
    
    global E; % efficiency
    
    global E_avr_1;
    global tau_FRET_avr_1;
    global E_avr_2;
    global tau_FRET_avr_2;
    
    N = numel(t);

    Ndecays = length(y)/N;

    tauD = alpha(1);
    RDA_div_R0_1 = alpha(2);
    RDA_div_R0_2 = alpha(3);

    Nb = numel(t);
    Np = numel(E);

    tau = tauD.*(1-E);

    pE_1 = efficiency_distribution(RDA_div_R0_1,E);
    pE_2 = efficiency_distribution(RDA_div_R0_2,E);

    FRETSUM_1 = zeros(1,Nb);
    FRETSUM_2 = zeros(1,Nb);
    
    pE_deficit_1 = 0;
    pE_deficit_2 = 0;
    for p = 1:Np        
        conv_irf = conv_irf_pp_69_70(tau(p));
        if 0~=sum(isnan(conv_irf))
            %disp([p tau(p)]);
            pE_deficit_1 = pE_deficit_1 + pE_1(p);
            pE_deficit_2 = pE_deficit_2 + pE_2(p);
            continue;
        end
        FRETSUM_1 = FRETSUM_1 + pE_1(p)/DT*conv_irf;
        FRETSUM_2 = FRETSUM_2 + pE_2(p)/DT*conv_irf;
    end      
    %
    A1 = FRETSUM_1/(sum(pE_1)-pE_deficit_1); % ??
    A2 = FRETSUM_2/(sum(pE_2)-pE_deficit_2); % ??
    
% standard convolution method    
% %     for p = 1:Np
% %         FRETSUM_1 = FRETSUM_1 + pE_1(p)*(1-exp(-DT/tau(p)))*exp(-t/tau(p))/(1-exp(-Tp/tau(p)));
% %         FRETSUM_2 = FRETSUM_2 + pE_2(p)*(1-exp(-DT/tau(p)))*exp(-t/tau(p))/(1-exp(-Tp/tau(p)));
% %     end
% %     
% %     A1 = FRETSUM_1/sum(pE_1);
% %     A2 = FRETSUM_2/sum(pE_2);    
% %     
% %     A1 = conv([ A1'; A1' ],IRF)';
% %     A1 = A1(N+1:2*N);
% %     A2 = conv([ A2'; A2' ],IRF)';
% %     A2 = A2(N+1:2*N);

    Phi = zeros(N*Ndecays,2*Ndecays);
        
    for k=1:Ndecays
        index = k*2-1;
        Phi((k-1)*N+1:k*N,index) = A1;
        Phi((k-1)*N+1:k*N,index+1) = A2;        
    end

    dPhi = [];
    Ind = [];
    %
    E_avr_1 = sum(pE_1.*E)/sum(pE_1);
    tau_FRET_avr_1 = sum(tau.*pE_1)/sum(pE_1);
    E_avr_2 = sum(pE_2.*E)/sum(pE_2);
    tau_FRET_avr_2 = sum(tau.*pE_2)/sum(pE_2);
    
    Phi = sparse(Phi);
    
end
