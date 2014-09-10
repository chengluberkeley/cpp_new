%% Date: June 4, 2014
%% Work for the nuclear detection project.
%% Data: San Francisco map

clc;
clear;

%% Load basic info of the problem
%load SFStreetGraph.mat; % INPUT: edges and nodes
% We may load the protection probabilities later!
load('./Santiago/edgeList.mat');
load('./Santiago/XY_coord.mat');
% Rename edges
edges = edgeList;
% Rename node coordinates
%xyCoord = nodes;
xyCoord = XY_coord;
nodes = xyCoord;
%load edgeProbs15.mat;
%load attackSet15.mat; % INPUT
% (INPUT) Utility vectors of defender and attacker. They each is a vector
%   corresponding to the nodes in the map
%load UtDefCov.mat;
%load UtDefUnc.mat;
%load UtAtkCov.mat;
%load UtAtkUnc.mat;

%% Constants for initialization
BIGNUM = inf;
<<<<<<< HEAD
filename = '08-30-2014.mat'; % Output filename

% Switch to determine whether it is a route or a tour
IS_BACK = 1;

% Maximum number of configurations to output
MAX_CONFIG_NUM = 20;

%% Obtain AOI (Areas of Interests)
AOI_CENTER = [0.1,0.5]; % X-Y Coordinates
AOI_RADIUS = 1000;
[edgesAOI, nodesAOI] = AOI(edgesSF, nodesSF, AOI_CENTER, AOI_RADIUS); % Set of target edges
MDA_LEVEL = 5; % Required uniform MDA level for all target edges.
edgeNumAOI = size(edgesAOI,1);
% Assume uniform speed limit over all the target edges
edgeVelocityAOI = zeros(edgeNumAOI,1) + 13.5; % 30 miles/hour ~= 13.5m/s
% Compute edge distances
edgeDistAOI = zeros(edgeNumAOI,1);
for i = 1:edgeNumAOI
    edgeIndex = edgesAOI(i);
    edgeDistAOI(i) = distMat(nodesSF(edgesSF(edgeIndex,1)), nodesSF(edgesSF(edgeIndex,2)));
end



%% Time interval constants
START_TIME_HOUR = 10;  % 24 hour format
START_TIME_MIN = 30;
END_TIME_HOUR = 11;
END_TIME_MIN = 30;
% Transform uniform format
START_TIME = START_TIME_HOUR*60 + START_TIME_MIN; % Unit: minute
END_TIME = END_TIME_HOUR*60 + END_TIME_MIN;
TIME_INTERVAL = END_TIME - START_TIME;  % Used to compute capacity

%% Types of vehicles are hard-coded now! Two types are considered currently:
% 21: Incentive taxis;
% 33: Fully controled trucks;

%% Initialization of datastructures for shortest path computing. This is only for trucks!
=======
filename = '07-17-2014_1_basic.mat'; % Output filename
PM_SWITCH = 1; % Switch to determine which minimum cost perfect matching subroutine to use.
               % 0: the downloaded exact algorithm (slow)
               % 1: the heuristic algortihm (fast)
% Switches to determine which part is run
%SIMPLE_TEST_SWITCH = false;
%SK_SWITCH = true;
%UK_SWITCH = true;
%SN_SWITCH = true;
%UN_SWITCH = true;

% Switch to do initial checking for the nearest neighbor algorithm
%NN_SWITCH = false;

%% Initialization
m = size(edges, 1); % Number of edges
n = size(nodes, 1); % Number of nodes.
%prob = edgeProbs; % (INPUT) Probabilities over edges. This is computed by solving the Stackelberg game (Need to prepare the data specifically.)
% Only for m = 115:
% prob = prob/sum(prob);
%xyCoord = XY_coord; % (INPUT) x-y coordinates of each nodes. Use to compute the Euclidean distances between nodes.

%[edgeTargetList, ~, probs] = find(prob); % Obtain the list of target edges. Assume no parallel edges. The target edges are defined to be those edges with positive probabilities.
%edgeTargetList = edgeList(edgeTargetList,:);
%numEdgeTarget = size(edgeTargetList, 1);
% FOR COMPARISON: Generate uniform distribution
%   Notice: this uniform distribution should range over all the existing
%   edges, not only the edgeTargetList from the Stackelberg game solution.
unifEdgeTargetList = edges;
unifNumEdgeTarget = m;
unifProbs = zeros(unifNumEdgeTarget, 1) + 1/unifNumEdgeTarget;
% Fernando
%MM = sum(probs);

% Compute the geometric distances
distMat = zeros(n,n) + BIGNUM;
for i = 1:n
    distMat(i,i) = 0;
end

for i = 1:m  % Symmetric case
    node1 = edges(i,1);
    node2 = edges(i,2);
    distMat(node1, node2) = sqrt((xyCoord(node1,1)-xyCoord(node2,1))^2+(xyCoord(node1,2)-xyCoord(node2,2))^2);
    distMat(node2, node1) = distMat(node1,node2);
end
 
>>>>>>> parent of 3b3563f... Code clearance. Make main.m as the entry file
% Prepare for the computation of the all-pair shortest paths between any two nodes (including target edges)
% We use Bellman-Ford Algorithm
% Initialization
spMat = zeros(n,n) + BIGNUM; % Initialize;
index = sub2ind(size(spMat), 1:n, 1:n);
spMat(index) = 0;
spRteMat = zeros(n,n) + diag(1:n); % Matrix to record the shortest path route info
hasDoneSp = zeros(n,1); % Record whether we have computed the shortest paths from the source j. 1: has computed.
% We only compute the shortest paths when needed.


%% Compare two possible potential routing strategies
% 1. Partition a single route into k routes (PK)
% 2. Add one target egde at a time, to the nearest truck that still has
% capacity to cover it. (NN)

% Scenario constants, subject to change. Note: now we only consider trucks!
DAYS = 100; % Total days to consider. Assume one attack per day
V_AVAILABLE = [3,7,11,15,30]; % This is the number of available vehicles each day;
%V_CAPACITY = [4000, 4250, 4500];  % Every day's capacity of each truck. Note that now every truck does not need to return to its respective starting point
% DEPOT_ID = 28; % Depot ID % No depot. Instead, let's assume that every
% day the k trucks will start at an edge uniformly!
N_RUN = 10; % Number of runs for the same configuration.
N_TARGET_EDGE = 20; % The number of target edges to protect.

%% Initialize the global variables to record the final success rates
% Here we don't have the "capacity". We only compute the shortest total distance (and the maximum single truck distance) to cover a specific set of target edges 
% Total distances
aveTotalDistPK = zeros(length(V_AVAILABLE),1); 
aveTotalDistNN = zeros(length(V_AVAILABLE),1);

% For standard deviation computing
sdTotalDistPK = zeros(length(V_AVAILABLE), N_RUN);
sdTotalDistNN = zeros(length(V_AVAILABLE), N_RUN);

% Maximum single truck distance
aveMaxDistPK = zeros(length(V_AVAILABLE),1);
aveMaxDistNN = zeros(length(V_AVAILABLE),1);

% For standard deviation computing
sdMaxDistPK = zeros(length(V_AVAILABLE), N_RUN);
sdMaxDistNN = zeros(length(V_AVAILABLE), N_RUN);


%% Run N_RUN number of different realization of the targets to defense/attack and compute the average MAXIMUM TOTAL DISTANCE/SINGLE TRUCK DISTANCE BY EACH METHOD!
for runtime = 1:N_RUN    
<<<<<<< HEAD
    % Base edge coverage matrix
    baseEdgeCoverCount = zeros(edgeNumAOI, 2);  % dimension 1: taxis; dimension 2; trucks.
    % Special case: no taxis
    taxi_num = 0;
    modeConfig(runtime, 1, 1) = taxi_num;
    % Compute the number of trucks needed in this case
    % Initialization
    targetEdgeIndex = (1:edgeNumAOI)';
    edgeCoverCount = baseEdgeCoverCount;  % This is only for the extreme case of having only trucks! 
    while (~isempty(targetEdgeIndex))        
        % Compute the number of trucks needed for covering all the
        % targetEdges
        
        % NN method
        % IS_NEXT == 1, IS_RETURN_COUNT == 0
        % We use a binary search algorithm to determine the
        % minimum number of trucks needed.
        high = 1;
        low = 0;
        % First find the high end of the search space
        shareLoc = zeros(high, 1) + truckLoc; % Need to recompute truckLoc.

        [leftFlag, spMat, spRteMat, hasDoneSp] = capNearestNeighbors(edgesSF, nodesSF, edgesAOI(targetEdgeIndex,:), spMat, spRteMat, hasDoneSp, high, TRUCK_CAPACITY, shareLoc, IS_BACK, 1, 0);
        while (leftFlag)
            low = high;
            high = high*2;
            shareLoc = zeros(high, 1)+truckLoc; % Need to recompute
            [leftFlag, spMat, spRteMat, hasDoneSp] = capNearestNeighbors(edgesSF, nodesSF, edgesAOI(targetEdgeIndex,:), spMat, spRteMat, hasDoneSp, high, TRUCK_CAPACITY, shareLoc, IS_BACK, 1, 0);
        end
=======
    %% Generate the potential list of edges to protect at every day!
    % Current simple case: Uniform distribution
    unifSumProb = zeros(unifNumEdgeTarget, 1);
    unifSofar = 0;
    for i = 1:unifNumEdgeTarget
        unifSumProb(i,1) = unifProbs(i) + unifSofar;
        unifSofar = unifSofar + unifProbs(i);
    end
>>>>>>> parent of 3b3563f... Code clearance. Make main.m as the entry file

    %% Generate the defenders' candidate target edges to defend for all the DAYS
    % Current simple case: Uniform distribution
    defenceVectorU = zeros(unifNumEdgeTarget,DAYS);
    for i = 1:DAYS
        varProb = unifProbs;
        varSumProb = unifSumProb;
        varEdgeTargetIndexList = (1:unifNumEdgeTarget)';
        for j = 1:unifNumEdgeTarget
            randDefProb = rand(1);
            index = f_whichisit(randDefProb, varSumProb);
            defenceVectorU(j,i) = varEdgeTargetIndexList(index);
            % Update the probabilities and the sum of probabilities:
            % conditional probabilities!
            varProb(index) = [];
            varEdgeTargetIndexList(index) = [];
            varSumProb(index) = [];
            varProb = varProb/(sum(varProb)); % Re-normalization
            sofar = 0;
            for k = 1:length(varSumProb)
                varSumProb(k,1) = varProb(k) + sofar;
                sofar = sofar + varProb(k);
            end
        end
<<<<<<< HEAD
        
        % Update the number of trucks needed
        modeConfig(runtime, 1, 2) = modeConfig(runtime, 1, 2) + high;
        % Update the edge coverages
        edgeCoverCount(targetEdgeIndex, 2) = edgeCoverCount(targetEdgeIndex, 2) + 1;
        
        % Check MDA and compute the new set of target edges to further
        % cover
        % Step 1: prepare the input to call the MDA computation subroutine
        totalCoverCount = sum(sum(edgeCoverCount));
        edgeIndex = zeros(totalCoverCount,1);
        edgeLength_m = zeros(totalCoverCount,1);
        edgeVelocity_m_per_s = zeros(totalCoverCount,1);
        modeId = zeros(totalCoverCount,1);
        mark = 1;
        for i = 1:edgeNumAOI
            for j = 1:edgeCoverCount(i,1)
                edgeIndex(mark) = edgesAOI(i);
                edgeLength_m(mark) = edgeDistAOI(i);
                edgeVelocity_m_per_s(mark) = edgeVelocityAOI(i);
                modeId(mark) = 21;
                mark = mark + 1;
            end
            
            for j = 1:edgeCoverCount(i,2)
                edgeIndex(mark) = edgesAOI(i);
                edgeLength_m(mark) = edgeDistAOI(i);
                edgeVelocity_m_per_s(mark) = edgeVelocityAOI(i);
                modeId(mark) = 33;
                mark = mark + 1;
            end
        end
        
        inStruct = struct('edgeIndex', edgeIndex, 'edgeLength_m', edgeLength_m, 'edgeVelocity_m_per_s', edgeVelocity_m_per_s, 'modeId', modeId);
        outStruct = mda(inStruct);
        
        
=======
>>>>>>> parent of 3b3563f... Code clearance. Make main.m as the entry file
    end

    %% Prepare the set of edges to protect every day
    defenceVector = defenceVectorU(1:N_TARGET_EDGE, :);
    
    %% Case 1: Uniform distribution + PK; Case 2: Uniform distribution + NN
    % Note that here we will put the two scenarios into the same for loops,
    % because at every day they need to share the same initial positions
    % for the k trucks
    for numV = 1:length(V_AVAILABLE)
        totalDistPK = zeros(1, DAYS);
        maxDistPK = zeros(1, DAYS);
        totalDistNN = zeros(1, DAYS);
        maxDistNN = zeros(1, DAYS);

        for day = 1:DAYS % Check day-by-day
            % Generate the initial points of the k trucks. Note that
            % here we assume that the trucks starting at nodes (intersection)
            truckLoc = randi(n, V_AVAILABLE(numV), 1);
            % PK method
            [totalDistPK(day), maxDistPK(day), spMat, spRteMat, hasDoneSp] = PK(distMat, unifEdgeTargetList(defenceVector(:, day), :), spMat, spRteMat, hasDoneSp, V_AVAILABLE(numV), truckLoc, edges);
            % NN method
            [totalDistNN(day), maxDistNN(day), spMat, spRteMat, hasDoneSp] = nearestNeighbors(distMat, unifEdgeTargetList(defenceVector(:, day), :), spMat, spRteMat, hasDoneSp, V_AVAILABLE(numV), truckLoc, edges);
        end

        % PK
        aveTotalDistPK(numV) = aveTotalDistPK(numV) + sum(totalDistPK)/DAYS;
        sdTotalDistPK(numV, runtime) = sum(totalDistPK)/DAYS;
        aveMaxDistPK(numV) = aveMaxDistPK(numV) + sum(maxDistPK)/DAYS;
        sdMaxDistPK(numV, runtime) = sum(maxDistPK)/DAYS;

        % NN
        aveTotalDistNN(numV) = aveTotalDistNN(numV) + sum(totalDistNN)/DAYS;
        sdTotalDistNN(numV, runtime) = sum(totalDistNN)/DAYS;
        aveMaxDistNN(numV) = aveMaxDistNN(numV) + sum(maxDistNN)/DAYS;
        sdMaxDistNN(numV, runtime) = sum(maxDistNN)/DAYS;    
    end
end


%% Compute the final result
% Compute averages
aveTotalDistPK = aveTotalDistPK/N_RUN;
aveMaxDistPK = aveMaxDistPK/N_RUN;
aveTotalDistNN = aveTotalDistNN/N_RUN;
aveMaxDistNN = aveMaxDistNN/N_RUN;

% Compute standard deviations
sdTPK = zeros(length(V_AVAILABLE));
sdMPK = zeros(length(V_AVAILABLE));
sdTNN = zeros(length(V_AVAILABLE));
sdMNN = zeros(length(V_AVAILABLE));

for numV = 1:length(V_AVAILABLE)
    for runtime = 1:N_RUN
        sdTPK(numV) = sdTPK(numV) + (sdTotalDistPK(numV,runtime) - aveTotalDistPK(numV))*(sdTotalDistPK(numV,runtime) - aveTotalDistPK(numV));
        sdMPK(numV) = sdMPK(numV) + (sdMaxDistPK(numV,runtime) - aveMaxDistPK(numV))*(sdMaxDistPK(numV,runtime) - aveMaxDistPK(numV));
        sdTNN(numV) = sdTNN(numV) + (sdTotalDistNN(numV,runtime) - aveTotalDistNN(numV))*(sdTotalDistNN(numV,runtime) - aveTotalDistNN(numV));
        sdMNN(numV) = sdMNN(numV) + (sdMaxDistNN(numV,runtime) - aveMaxDistNN(numV))*(sdMaxDistNN(numV,runtime) - aveMaxDistNN(numV));        
    end
end

for numV = 1:length(V_AVAILABLE)
    sdTPK(numV) = sqrt(1/(N_RUN-1)*sdTPK(numV));
    sdMPK(numV) = sqrt(1/(N_RUN-1)*sdMPK(numV));
    sdTNN(numV) = sqrt(1/(N_RUN-1)*sdTNN(numV));
    sdMNN(numV) = sqrt(1/(N_RUN-1)*sdMNN(numV));

end

% Write to file
save(filename, 'V_AVAILABLE', 'aveTotalDistPK', 'aveTotalDistNN', 'aveMaxDistPK', 'aveMaxDistNN', 'sdTPK', 'sdTNN', 'sdMPK', 'sdMNN');
