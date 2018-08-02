function AM = adjacency_matrix(x,y,mode)

N = numel(x);
AM = zeros(N);

if 1==strcmp(mode,'Gabriel')
    
    for k=1:N
        for j=1:N
            if k<j
                AM(k,j) = 1;
                AM(j,k) = 1;
                v_k = [x(k) y(k)];
                v_j = [x(j) y(j)];
                v_c = (v_k+v_j)/2;
                r = norm(v_k-v_c);
                %
                for m=1:N
                    if m~=k && m~=j
                        v_m = [x(m) y(m)];
                        if norm(v_m-v_c)<r %not neighbours
                            AM(k,j) = 0;
                            AM(j,k) = 0;                           
                            break;
                        end
                    end
                end
            end
        end
    end    
end    

if 1==strcmp(mode,'SOI')
   
    % calculate min distances    
    d_min = zeros(1,N);    
    for k=1:N
        v_k = [x(k) y(k)];                        
        D = 1e27;
            for j=1:N
                if k~=j
                    v_j = [x(j) y(j)];    
                    d_kj = norm(v_k-v_j);
                    %
                    if d_kj < D
                        D = d_kj;
                    end
                end
            end
        d_min(1,k) = D;
    end

    for k=1:N
        v_k = [x(k) y(k)];                        
        for j=1:N
            if k~=j
                v_j = [x(j) y(j)];    
                d_kj = norm(v_k-v_j);
                if d_kj < ( d_min(k) + d_min(j) ) %neighbours
                    AM(k,j)=1;
                    AM(j,k)=1;
                end
            end
        end
    end        
end

if 1==strcmp(mode,'Delaunay')    
        
    if 2==length(x)
        AM = [ 0 1 ; 1 0];
    end
                
    if length(x)>2
        
    TRI = delaunay(x,y);    
    [L_TRI, tri] = size(TRI);

    for k=1:L_TRI
        i1 = TRI(k,1);
        i2 = TRI(k,2);
        i3 = TRI(k,3);
        AM(i1,i2)=1;
        AM(i2,i1)=1;    
        AM(i1,i3)=1;
        AM(i3,i1)=1;    
        AM(i2,i3)=1;    
        AM(i3,i2)=1;        
    end
    
    end
    
end

% % figure();
% % plot(x,-y,'r.');
% %     for k=1:N
% %         for j=1:N
% %             if k<j && 1==AM(k,j)
% %                 v1 = [x(k) x(j)];
% %                 v2 = [y(k) y(j)];
% %                 line(v1,-v2);
% %             end
% %         end
% %     end
% % end
