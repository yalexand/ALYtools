
function [Phi,dPhi,Ind] = adaex_sk2_only_FLIMFRET(alpha)

    global t; % picoseconds

    global DT;
        
    global IRF;  
    
    global Tp;
    
    global y;
    
    global E; % efficiency
    
    global E_avr;
    global tau_FRET_avr;
        
    N = numel(t);

    Ndecays = length(y)/N;

    tauD = alpha(1);
    RDA_div_R0 = alpha(2);

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
    
    A1 = conv([ A1'; A1' ],IRF)';
    A1 = A1(N+1:2*N);

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
