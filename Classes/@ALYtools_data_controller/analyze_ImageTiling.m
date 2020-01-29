 function [datas, captions, table_names, fig] = analyze_ImageTiling(obj,~,~) 

tic

     datas = [];
     captions = [];
     table_names = 'default';
     fig = [];
     
     sgm = obj.do_ImageTiling_Segmentation(false);
  
     [sX,sY,sC,sZ,sT] = size(sgm);
     
     rws = obj.ImageTiling_Nrows;
     cls = obj.ImageTiling_Ncols;
     prct_x = obj.ImageTiling_Ovlp_X;
     prct_y = obj.ImageTiling_Ovlp_Y;
     QT = obj.ImageTiling_QT;

     N = rws;
     M = cls;
     
     if 0==N || 0==M
         disp('error - number of row/columns can not be zero')
         return;
     end
     
     sT = N*M;
     
     A = zeros(N*M);
     
     for k=1:N
        for m=1:M
            I0 = sub2ind([N M],k,m);
            %
            I1 = 0;
            I2 = 0;
            I3 = 0;
            I4 = 0;
            %
            try I1 = sub2ind([N M],k-1,m); catch; end
                if 0~=I1, A(I0,I1) = 1; A(I1,I0) = 1; end
            try I2 = sub2ind([N M],k+1,m); catch; end
                if 0~=I2, A(I0,I2) = 1; A(I2,I0) = 1; end
            try I3 = sub2ind([N M],k,m-1); catch; end
                if 0~=I3, A(I0,I3) = 1; A(I3,I0) = 1; end
            try I4 = sub2ind([N M],k,m+1); catch; end
                if 0~=I4, A(I0,I4) = 1; A(I4,I0) = 1; end
            %                                
        end
     end
  
    G = graph(A);
%    T = minspantree(G);
%     figure;
%     a1=subplot(1,2,1);
%     plot(a1,G);
%     title(a1,num2str(numedges(G)));
%     a2=subplot(1,2,2);
%     plot(a2,T);
%     title(a2,num2str(numedges(T)));

    edges = table2cell(G.Edges);

    Nedges = size(edges,1);

    edges_shifts_xy = zeros(Nedges,2);
    edges_q = zeros(Nedges,1);

    % is this needed?
    nodes = [];
    for k=1:size(edges,1)
        edge = cell2mat(edges(k));
        nodes = [nodes edge(:)];
    end
    nodes = unique(nodes);
    % is this needed?    
   
    Dx = fix(sX*prct_x);
    Dy = fix(sY*prct_y);
    rx1 = 1:Dx;
    rx2 = sX-Dx+1:sX;
    ry1 = 1:Dy;
    ry2 = sY-Dy+1:sY;
    rx_full=1:sX;
    ry_full=1:sY;    

    along_row_acc_x1x2 = [];
    along_row_acc_x2x1 = [];
    along_row_acc_y1y2 = [];
    along_row_acc_y2y1 = [];
    %
    along_col_acc_x1x2 = [];
    along_col_acc_x2x1 = [];
    along_col_acc_y1y2 = [];
    along_col_acc_y2y1 = [];        
    %
    s1 = 3/2;
    s2 = 7/2;
    for k=1:Nedges
        edge = edges{k,1};
        in_node = edge(1);
        out_node = edge(2);
        disp([k in_node out_node]);
        %
        u1 = squeeze(sgm(:,:,1,1,in_node));
        u2 = squeeze(sgm(:,:,1,1,out_node));
        %
        %
        [r_in,col_in] = ind2sub([N M],in_node);
        [r_out,col_out] = ind2sub([N M],out_node);
        dr = abs(r_in-r_out);
        mode = 'along row';
        if 0~=dr
             mode = 'along col';
        end
        disp([num2str(in_node) ' ' num2str(out_node) ' ' mode]);
        if strcmp(mode,'along row')
            % x1x2
                    z = xcorr2_fft(imresize(u1(rx1,ry_full),0.5),imresize(u2(rx2,ry_full),0.5));
                    g1 = gsderiv(z,s1,0);
                    g2 = gsderiv(z,s2,0);
                    z = (g1-g2);                    
                    z(z<0)=0;                    
                    q = max(z(:))/mean(z(:));
                    along_row_acc_x1x2 = [along_row_acc_x1x2; q];
                    %icy_imshow(uint16(z),[num2str(k) ' : ' num2str(in_node) ' -> ' num2str(out_node) ', q = ' num2str(q) ' x1x2']);
            % x2x1
                    z = xcorr2_fft(imresize(u1(rx2,ry_full),0.5),imresize(u2(rx1,ry_full),0.5));
                    g1 = gsderiv(z,s1,0);
                    g2 = gsderiv(z,s2,0);
                    z = (g1-g2);                    
                    z(z<0)=0;                    
                    q = max(z(:))/mean(z(:));
                    along_row_acc_x2x1 = [along_row_acc_x2x1; q];            
                    %icy_imshow(uint16(z),[num2str(k) ' : ' num2str(in_node) ' -> ' num2str(out_node) ', q = ' num2str(q) ' x2x1']);                    
            % y1y2
                    z = xcorr2_fft(imresize(u1(rx_full,ry1),0.5),imresize(u2(rx_full,ry2),0.5));
                    g1 = gsderiv(z,s1,0);
                    g2 = gsderiv(z,s2,0);
                    z = (g1-g2);                    
                    z(z<0)=0;                    
                    q = max(z(:))/mean(z(:));
                    along_row_acc_y1y2 = [along_row_acc_y1y2; q];            
                    %icy_imshow(uint16(z),[num2str(k) ' : ' num2str(in_node) ' -> ' num2str(out_node) ', q = ' num2str(q) ' y1y2']);                    
            % y2y1
                    z = xcorr2_fft(imresize(u1(rx_full,ry2),0.5),imresize(u2(rx_full,ry1),0.5));
                    g1 = gsderiv(z,s1,0);
                    g2 = gsderiv(z,s2,0);
                    z = (g1-g2);                    
                    z(z<0)=0;                    
                    q = max(z(:))/mean(z(:));
                    along_row_acc_y2y1 = [along_row_acc_y2y1; q];                        
                    %icy_imshow(uint16(z),[num2str(k) ' : ' num2str(in_node) ' -> ' num2str(out_node) ', q = ' num2str(q) ' y2y1']);                    
        elseif strcmp(mode,'along col')
            % x1x2
                    z = xcorr2_fft(imresize(u1(rx1,ry_full),0.5),imresize(u2(rx2,ry_full),0.5));
                    g1 = gsderiv(z,s1,0);
                    g2 = gsderiv(z,s2,0);
                    z = (g1-g2);                    
                    z(z<0)=0;                    
                    q = max(z(:))/mean(z(:));
                    along_col_acc_x1x2 = [along_col_acc_x1x2; q];
                    %icy_imshow(uint16(z),[num2str(k) ' : ' num2str(in_node) ' -> ' num2str(out_node) ', q = ' num2str(q) ' x1x2']);
            % x2x1
                    z = xcorr2_fft(imresize(u1(rx2,ry_full),0.5),imresize(u2(rx1,ry_full),0.5));
                    g1 = gsderiv(z,s1,0);
                    g2 = gsderiv(z,s2,0);
                    z = (g1-g2);                    
                    z(z<0)=0;                    
                    q = max(z(:))/mean(z(:));
                    along_col_acc_x2x1 = [along_col_acc_x2x1; q];            
                    %icy_imshow(uint16(z),[num2str(k) ' : ' num2str(in_node) ' -> ' num2str(out_node) ', q = ' num2str(q) ' x2x1']);                    
            % y1y2
                    z = xcorr2_fft(imresize(u1(rx_full,ry1),0.5),imresize(u2(rx_full,ry2),0.5));
                    g1 = gsderiv(z,s1,0);
                    g2 = gsderiv(z,s2,0);
                    z = (g1-g2);                    
                    z(z<0)=0;                    
                    q = max(z(:))/mean(z(:));
                    along_col_acc_y1y2 = [along_col_acc_y1y2; q];            
                    %icy_imshow(uint16(z),[num2str(k) ' : ' num2str(in_node) ' -> ' num2str(out_node) ', q = ' num2str(q) ' y1y2']);                    
            % y2y1
                    z = xcorr2_fft(imresize(u1(rx_full,ry2),0.5),imresize(u2(rx_full,ry1),0.5));
                    g1 = gsderiv(z,s1,0);
                    g2 = gsderiv(z,s2,0);
                    z = (g1-g2);                    
                    z(z<0)=0;                    
                    q = max(z(:))/mean(z(:));
                    along_col_acc_y2y1 = [along_col_acc_y2y1; q];                        
                    %icy_imshow(uint16(z),[num2str(k) ' : ' num2str(in_node) ' -> ' num2str(out_node) ', q = ' num2str(q) ' y2y1']);                                
        end        
        %
    end
        
%     figure; 
%     subplot(2,2,1)
%     hist(along_row_acc_x1x2); title('row x1x2');
%     subplot(2,2,2)    
%     hist(along_row_acc_x2x1); title('row x2x1');
%     subplot(2,2,3)    
%     hist(along_row_acc_y1y2); title('row y1y2');
%     subplot(2,2,4)    
%     hist(along_row_acc_y2y1); title('row y2y1');
%     %
%     figure; 
%     subplot(2,2,1)
%     hist(along_col_acc_x1x2); title('col x1x2');
%     subplot(2,2,2)    
%     hist(along_col_acc_x2x1); title('col x2x1');
%     subplot(2,2,3)    
%     hist(along_col_acc_y1y2); title('col y1y2');
%     subplot(2,2,4)    
%     hist(along_col_acc_y2y1); title('col y2y1');
    
    row_tokens = {'row x1x2','row x2x1','row y1y2','row y2y1'}
    row_q = [median(along_row_acc_x1x2) median(along_row_acc_x2x1) median(along_row_acc_y1y2) median(along_row_acc_y2y1)]
    col_tokens = {'col x1x2','col x2x1','col y1y2','col y2y1'}
    col_q = [median(along_col_acc_x1x2) median(along_col_acc_x2x1) median(along_col_acc_y1y2) median(along_col_acc_y2y1)]   
    try
                row_token = row_tokens{row_q==max(row_q(:))};    
                switch row_token
                    case 'row x1x2' % z = xcorr2_fft(u1(rx1,ry_full),u2(rx2,ry_full));
                        rRX1 = rx1;
                        rRY1 = ry_full;
                        rRX2 = rx2;
                        rRY2 = ry_full;
                    case 'row x2x1' % z = xcorr2_fft(u1(rx2,ry_full),u2(rx1,ry_full));
                        rRX1 = rx2;
                        rRY1 = ry_full;
                        rRX2 = rx1;
                        rRY2 = ry_full;
                    case 'row y1y2' % z = xcorr2_fft(u1(rx_full,ry1),u2(rx_full,ry2));
                        rRX1 = rx_full;
                        rRY1 = ry1;
                        rRX2 = rx_full;
                        rRY2 = ry2;
                    case 'row y2y1' % z = xcorr2_fft(u1(rx_full,ry2),u2(rx_full,ry1));                        
                        rRX1 = rx_full;
                        rRY1 = ry2;                        
                        rRX2 = rx_full;
                        rRY2 = ry1;                                                
                end
    catch
        disp('no rows');
    end
    
    try
                col_token = col_tokens{col_q==max(col_q(:))};
                switch col_token
                    case 'col x1x2' % z = xcorr2_fft(u1(rx1,ry_full),u2(rx2,ry_full));
                        cRX1 = rx1;
                        cRY1 = ry_full;
                        cRX2 = rx2;
                        cRY2 = ry_full;
                    case 'col x2x1'
                        cRX1 = rx2;
                        cRY1 = ry_full;
                        cRX2 = rx1;
                        cRY2 = ry_full;
                    case 'col y1y2'
                        cRX1 = rx_full;
                        cRY1 = ry1;
                        cRX2 = rx_full;
                        cRY2 = ry2;
                    case 'col y2y1'                        
                        cRX1 = rx_full;
                        cRY1 = ry2;                        
                        cRX2 = rx_full;
                        cRY2 = ry1;                                                
                end
    catch
        disp('no columns');
    end
                    
    s1 = 3;
    s2 = 7;
    
    % define edges characteristics
    for k=1:Nedges
        edge = edges{k,1};
        in_node = edge(1);
        out_node = edge(2);
        disp([k in_node out_node]);
                u1 = squeeze(sgm(:,:,1,1,in_node));
                u2 = squeeze(sgm(:,:,1,1,out_node));                                                   
                    %
        [r_in,col_in] = ind2sub([N M],in_node);
        [r_out,col_out] = ind2sub([N M],out_node);
        dr = abs(r_in-r_out);
        mode = 'along row';
        if 0~=dr
             mode = 'along col';
        end
        if strcmp(mode,'along col')
            RX1 = cRX1;
            RY1 = cRY1;
            RX2 = cRX2;
            RY2 = cRY2;
        else
            RX1 = rRX1;
            RY1 = rRY1;
            RX2 = rRX2;
            RY2 = rRY2;            
        end                                
        %
        z = xcorr2_fft(u1(RX1,RY1),u2(RX2,RY2)); 
        %
        g1 = gsderiv(z,s1,0);
        g2 = gsderiv(z,s2,0);
        z = (g1-g2);
        %
        z(z<0)=0;
        %
        [wc,hc] = size(z);
        wc=fix(wc/2);
        hc=fix(hc/2);
        rx = wc-s2:wc+s2;
        ry = hc-s2:hc+s2;
        z(rx,ry)=0;    
        %
        [x,y] = find(z==max(z(:)));
        %
        edges_q(k) = z(x,y)/mean(z(:)); % quality                                        
        d_shift = [x-wc y-hc];        
        edges_shifts_xy(k,:) = d_shift + [max(RX1)-sX max(RY1)-sY];
%                     if edges_q(k)<100
%                     icy_imshow(uint16(z),[num2str(k) ' : ' num2str(in_node) ' -> ' num2str(out_node) ', q = ' num2str(edges_q(k))]);
%                     end
    end
                          
%chosen = N*M/2+M/2;
%chosen = 64;
chosen = 1;

linked_nodes = chosen;
coord = nan(numel(nodes),2); % X,Y
coord(chosen,:) = [0 0];

% compile the list of bad (unreachable) nodes
bad_nodes = [];
for n=1:numel(nodes)
    node = nodes(n);
    n_is_bad = true;
    for e = 1:numel(edges)
        if ismember(node,edges{e}) && edges_q(e) > QT
            n_is_bad = false;
            break;
        end
    end
    if n_is_bad 
        bad_nodes = [bad_nodes node];
    end
end

while numel(linked_nodes) < numel(nodes) - numel(bad_nodes)
    %
    for k=1:Nedges
        %
        edge = edges{k};
        %
        in_node = edge(1);
        out_node = edge(2);
        %
        if ismember(in_node,bad_nodes) || ismember(out_node,bad_nodes), continue, end        
        %
        if edges_q(k) < QT, continue, end
        %
        in_node_linked = ~isempty(intersect(in_node,linked_nodes));
        out_node_linked = ~isempty(intersect(out_node,linked_nodes));
        %
        if in_node_linked && ~out_node_linked % use shift directly
            %
            coord(out_node,:) = coord(in_node,:) + edges_shifts_xy(k,:);
            %
            linked_nodes = [linked_nodes out_node];
            break;
            %
        elseif out_node_linked && ~in_node_linked % invert the shift
            %
            coord(in_node,:) = coord(out_node,:) - edges_shifts_xy(k,:);
            %
            linked_nodes = [linked_nodes in_node];
            break;
        else
            continue, 
        end           
    end    
end

coord = coord - min(coord);

% debug
% figure;plot(coord(:,1),coord(:,2),'bo')        

SX = max(coord(:,1))+sX;
SY = max(coord(:,2))+sY;
scene = zeros(SX,SY);
CX = 1;
CY = 1;
%
if strcmp(obj.ImageTiling_mode,'bleached_fluor')

    % %%%%%%%%%%%%%%%%%%%%%%%%% beautify
    if sT<10
        scene = zeros(SX,SY,sT,1,1);
        for k=1:sT
            if ~isnan(coord(k,1))
                rx = coord(k,1) + CX : coord(k,1) + CX+sX-1;
                ry = coord(k,2) + CY : coord(k,2) + CY+sY-1;
                scene(rx,ry,k,1,1) = sgm(:,:,2,1,k);
            end
        end
        try
            icy_imshow(uint16(scene));
        catch
            disp('cannot show image in icy - not running?')
        end
    end
    % %%%%%%%%%%%%%%%%%%%%%%%%%    
 
% primitive "blending"    
scene = zeros(SX,SY);
        for k=1:sT
            if ~isnan(coord(k,1))            
                rx = coord(k,1) + CX : coord(k,1) + CX+sX-1;
                ry = coord(k,2) + CY : coord(k,2) + CY+sY-1;
                scene(rx,ry) = pixelwise_max(scene(rx,ry),sgm(:,:,2,1,k));
            end
        end
%    
elseif strcmp(obj.ImageTiling_mode,'brightfield')
    % median along all acquired images for brightfield - to compensate camera dirt - thanks Sunil
    z = squeeze(sgm(:,:,2,1,:));
    medsgm = median(z,3);
    for k=1:sT
        if ~isnan(coord(k,1))
            rx = coord(k,1) + CX : coord(k,1) + CX+sX-1;
            ry = coord(k,2) + CY : coord(k,2) + CY+sY-1;    
            %scene(rx,ry) = sgm(:,:,2,1,k);
            scene(rx,ry) = sgm(:,:,2,1,k) - medsgm;
        end
    end
    scene = scene - min(scene(:));
end

% debug
%icy_imshow(scene);

fig = uint16(scene);    

for k=1:sT
    record = {char(obj.M_filenames{k}), coord(k,1), coord(k,2)};
    datas = [datas; record];
end
captions = {'filename','X','Y'};

disp(['execution time ' num2str(toc/60) ' min']);

end

function r3 = rect_x_rect(x1,y1,w1,h1,x2,y2,w2,h2)
    %r1 = [x1 y1 w1 h1];
    %r2 = [x2 y2 w2 h2];
    l = max(x1,x2);
    r = min(x1+w1,x2+w2);
    b = max(y1,y2);
    t = min(y1+h1,y2+h2);
    r3 =  [l b r-l t-b];
end



