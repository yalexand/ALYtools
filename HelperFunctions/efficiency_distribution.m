function pE = efficiency_distribution(RDA_div_R0,E)

    pE = zeros(1,numel(E));

    F =  3/2*(RDA_div_R0)^(-6);

    G = sqrt(E/F./(1-E));

    H = (2.*(1-E).*sqrt(3.*E.*F.*(1-E))).^(-1);

    for k=1:numel(E) 
        E_k = E(k);
        %
        if E_k >= 0 && E_k <= F/(1+F)
            pE(k) = H(k).*log(2+sqrt(3));
        elseif E_k > F/(1+F) && E_k < 4*F/(1+4*F)
            pE(k) = H(k).*log( (2+sqrt(3))./(G(k) + sqrt(G(k).^2-1)) );
        else
            pE(k) = 0;        
        end
        %
    end
        
end
