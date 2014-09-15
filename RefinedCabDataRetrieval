function [ListofRoutes]= RefinedCabDataRetrieval(edges,nodes,center,radius,numberofcabs)
constant = 2; 
threshold = radius * constant; 
counter = 0; %% counter for how many cabs have been selected
CabIDs=[];
while(counter<numberofcabs)
    Data = QueryDiscreteCabData('OSMSF.mat',1,10,5,2);
    if ismember(Data.id(1),CabIDs)==1 %% check to see that the Cab you pick hasn't already been selected
        continue;
    end
    F0 = find(Data.mflag == 0); %% returns list of all unmetered edges of the Cab
    [ListofEdges,ListofNodes] = AOI(edges,nodes,center,threshold); %% returns the List of Edges within the threshold area
    Endingindex=0;
    listofRoutesindex=1;
    NumofRoutes=0; %% number of routes generated 
    for i=1:length(F0)-1
        if(i<=Endingindex)
            continue;
        end
            
        if(ismember(Data.edge(F0(i)),ListofEdges)==1) %% check to see if the unmetered edge is within the threshold region
            NumofRoutes=NumofRoutes+1;
            Startingedge=F0(i);
            while((F0(i+1))==F0(i)+1) %% generate the path of the taxi while its unmetered
                if (i >= length(F0)-1)
                    i=i+1;
                    break;
                else
                i=i+1;
                end
            end
            Endingindex=i;
            Endingedge=F0(i);
            ListofRoutes(listofRoutesindex,1)= Data.id(1);
            ListofRoutes(listofRoutesindex,2)=Startingedge;
            ListofRoutes(listofRoutesindex,3)=Endingedge;
            listofRoutesindex=listofRoutesindex+1;
        end
    end
    if(NumofRoutes >0)
        counter=counter+1;
        %%CabIDs(end+1)=Data.id(1);
    end
end 

end

