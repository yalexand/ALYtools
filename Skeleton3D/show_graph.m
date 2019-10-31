function show_graph(node2,link2,w,l,h,TITLE) % :)
figure();
title(TITLE);
hold on;
for i=1:length(node2)
    x1 = node2(i).comx;
    y1 = node2(i).comy;
    z1 = node2(i).comz;
    
    if(node2(i).ep==1)
        ncol = 'c';
    else
        ncol = 'y';
    end;
    
    for j=1:length(node2(i).links)    % draw all connections of each node
        if(node2(node2(i).conn(j)).ep==1)
            col='k'; % branches are black
        else
            col='r'; % links are red
        end;
        if(node2(i).ep==1)
            col='k';
        end;
        
        % draw edges as lines using voxel positions
        for k=1:length(link2(node2(i).links(j)).point)-1            
            [x3,y3,z3]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k));
            [x2,y2,z2]=ind2sub([w,l,h],link2(node2(i).links(j)).point(k+1));
            line([y3 y2],[x3 x2],[z3 z2],'Color',col,'LineWidth',2);
        end;
    end;
    
    % draw all nodes as yellow circles
    plot3(y1,x1,z1,'o','Markersize',9,...
        'MarkerFaceColor',ncol,...
        'Color','k');
end;
axis image;axis off;
set(gcf,'Color','white');
drawnow;
view(-17,46);

end