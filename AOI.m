function [ListofEdges, ListofNodes] = AOI(edges,nodes,center,radius)
ListofEdges=[];
ListofNodes=[];
for i=1:length(edges) % loop through all the edges
        a= edges(i,1); % find node a of the edge 
        b= edges(i,2); % find node b of the edge
        x1=nodes(a,1); % getting the (x,y) coord of node a
        y1=nodes(a,2);
        x2=nodes(b,1); % getting the (x,y) coord of node b
        y2=nodes(b,2);
        X=[x1,y1 ; center(1), center(2)];
        Y=[x2,y2 ; center(1), center(2)];
        D1 = pdist(X); % find the euclidean distance from the center to node a of edge i
        D2 = pdist(Y); % find the euclidean distance from the center to node b of edge i
    
            if (D1 <= radius) || (D2 <= radius) % if either distance is less than the radius add the edge to the list
                ListofEdges(end+1)= i;
            end
end
for i=1:length(ListofEdges) % loop through the edges within the radius
        node1=edges(ListofEdges(i),1);
        ListofNodes(end+1)=node1;
        node2=edges(ListofEdges(i),2);
        ListofNodes(end+1)=node2; % add the nodes of each edge to ListofNodes 
end
ListofNodes=unique(ListofNodes);
end 
    
