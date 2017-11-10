
function [Phi,dPhi,Ind] = adaex_sk2_FLIMFRET_fixed_tauD(alpha)

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
    
    for p = 1:Np
        FRETSUM = FRETSUM + pE(p)*(1-exp(-DT/tau(p)))*exp(-t/tau(p))/(1-exp(-Tp/tau(p)));
        %FRETSUM = FRETSUM + pE(p)/tau(p)*exp(-t/tau(p))/(1-exp(-Tp/tau(p)));
    end
    
    A1 = FRETSUM/sum(pE);                
    A2 = (1-exp(-DT/tauD))/(1-exp(-Tp/tauD)).*exp(-t/tauD);
    %A2 = 1/tauD/(1-exp(-Tp/tauD)).*exp(-t/tauD);
    
    A1 = conv([ A1'; A1' ],IRF)';
    A1 = A1(N+1:2*N);
    A2 = conv([ A2'; A2' ],IRF)';
    A2 = A2(N+1:2*N);

    Phi = sparse(zeros(N*Ndecays,2*Ndecays));
        
    for k=1:Ndecays
        index = k*2-1;
        Phi((k-1)*N+1:k*N,index) = A1;
        Phi((k-1)*N+1:k*N,index+1) = A2;        
    end

    dPhi = [];
    Ind = [];
    %
    E_avr = sum(pE.*E)/sum(pE);
    tau_FRET_avr = sum(tau.*pE)/sum(pE);
    
end
