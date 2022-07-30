% Yiwen Mei (yiwenm@umich.edu)
% SEAS, University of Michigan
% Last update: 7/2/2022

%% Functionality
% This function implements the concept of minimal depth and interactive depth
%  to quantify predictor importance for Random Forest (RF) models (Ishwaran et al. 2010;
%  Konapala & Mishra 2020).

%% Input
% Mdl: a compact treebagger object for the RF model.

%% Output
% MD : minimal depth of predictor to the root node;
% NMD: normalized minimal depth;
% ID : minimal depth of predictor to the root node of maximal subtree of another
%       predictor (Interactive Depth);
% NID: normazlized Interactive Depth.

function [MD,NMD,ID,NID]=RF_MinDepth(Mdl)
%% Check the inputs
narginchk(1,1);
ips=inputParser;
ips.FunctionName=mfilename;

addRequired(ips,'Mdl',@(x) validateattributes(x,{'CompactTreeBagger'},{'nonempty'},mfilename,'Mdl'));

parse(ips,Mdl);
clear ips

%% Calculation
MD=nan(length(Mdl.Trees),length(Mdl.PredictorNames));
NMD=nan(length(Mdl.Trees),length(Mdl.PredictorNames));
NID=nan(length(Mdl.PredictorNames),length(Mdl.PredictorNames),length(Mdl.Trees));
ID=nan(length(Mdl.PredictorNames),length(Mdl.PredictorNames),length(Mdl.Trees));
for j=1:length(Mdl.Trees)
  tr=Mdl.Trees{j};
  TrD=TreeDepth(tr); % Depth of tree
  LeafNode_ls=find(~tr.IsBranchNode);

  Np_nm=unique(tr.CutPredictor); % List of variable used in the tree
  Np_nm(cellfun(@isempty,Np_nm))=[];

  for n=1:length(Np_nm) % For every predictor
%% Minimal Depth
    np_id=strcmp(tr.PredictorNames,Np_nm{n}); % Predictor ID

    k=strcmp(tr.CutPredictor,Np_nm{n}); % Cut point Booleans of the predictor
    PrnNode_id=tr.Parent(k); % Parent node ID of the cut point
    SubTrD=zeros(size(PrnNode_id)); % Initialize the sub-tree depth
    while all(PrnNode_id~=0) % Find the root node upward
      SubTrD=SubTrD+ones(size(PrnNode_id));
      PrnNode_id=tr.Parent(PrnNode_id);
    end

    [SubTrD_min,SubTrD_i]=min(SubTrD);
    MD(j,np_id)=SubTrD_min; % minimal depth in distance
    NMD(j,np_id)=SubTrD_min/TrD; % Normalize to [0 1], larger more important

%% Interactive depth
% Create node structure for the connectivity
    Pid=find(k); % Cut point IDs of a predictor
    Pid=Pid(SubTrD_i); % Cut point ID of the maximal subtree defined by the predictor
    Nd_stru=[]; % Node structure (connectivity)
    while any(Pid~=0)
      ChlNode=[];
      for i=1:length(Pid)
        if Pid(i)~=0
          ChlNode_id=tr.Children(Pid(i),:)'; % The 2 child nodes of every cut point
          Nd_stru=[Nd_stru;[repelem(Pid(i),2,1) ChlNode_id]];
          ChlNode=[ChlNode;ChlNode_id]; % Collection of child nodes for the cut points
        end
      end
      Pid=ChlNode;
    end
    k=ismember(Nd_stru(:,2),[LeafNode_ls;0]);
    Nd_stru(k,:)=[]; % Remove the leaf nodes

% Calculate the shortest distance from root node of subtree to other nodes
    if ~isempty(Nd_stru)
      Gobj=digraph(Nd_stru(:,1),Nd_stru(:,2)); % Create a directed graph object
      [~,IntD]=shortestpathtree(Gobj,Nd_stru(1,1),Nd_stru(:,2));
      MTD=max(IntD)+1; % Maximal depth of subtree

% Convert node structure from cut point ID to predictor ID
      Nd_stru_pid=cellfun(@(X) find(strcmp(tr.PredictorNames,X)),tr.CutPredictor(Nd_stru(:,2)));
      IntD=accumarray(Nd_stru_pid,IntD',fliplr(size(Mdl.PredictorNames)),@min,NaN); % If a combination
      ID(np_id,:,j)=IntD;       % does not exist, that means no interaction between the two variables.
      NID(np_id,:,j)=IntD/MTD;  % In this case, set ID=TrD and NID=1.
    end
  end
end

%% Average for the trees
MD=mean(MD,1,'omitnan');
NMD=mean(NMD,1,'omitnan');
ID=mean(ID,3,'omitnan');
ID(isnan(ID))=Inf;
NID=mean(NID,3,'omitnan');
NID(isnan(NID))=1;
end

function TrD=TreeDepth(tree)
if isa(tree,'classreg.learning.regr.CompactRegressionTree')
  parent=tree.Parent;
else
  parent=tree;
end

TrD=0;
node=parent(end);
while node~=0
  TrD=TrD+1;
  node=parent(node);
end
end
