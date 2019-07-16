
%% Last Edit 07/05/2017
clear all
close all


%% Initialization
%delete(gcp('nocreate'))
global NI
global AI
global IC
global ControlCode
NI = [];
AI = [];
% Initialize Wc.m (Used in "CombineColumn", call at line 122)
FP = fopen('Wc.m','w');
fprintf(FP,'Wci = [];');
fclose(FP);
%
%% Create the robot dynamic model in symbolic form
% Computes the model, linear in the parameters, using data previously obtained in Mathematica
disp('Creating the symbolic dynamic model...');
tic
CreateModel
disp('Done.');
toc

n = size(Wi,1); % Number of DOFs
N = 100;         % Number of Data Points in the Trajectory

%% Optimization Variables
Bq = [-3,-3;3,3];       % Pos. Bounds
Bdq = [-5,-5;5,5];      % Vel. Bounds
Bddq = [-10,-10;10,10]; % Acc. Bounds
%
LB = zeros(N,n*3)+[Bq(1,:),Bdq(1,:),Bddq(1,:)];
LB = reshape(LB,[1,3*n*N]);
UB = zeros(N,n*3)+[Bq(2,:),Bdq(2,:),Bddq(2,:)];
UB = reshape(UB,[1,3*n*N]);
% Configure Optimization Algorithm
UsePar = 1;
%% Start the Hybrid optimization of the 1st Trajectory (PSW,GA or Fmincon)
% Set up Parallel Computation
if UsePar
    disp('Enstablishing workers for parallel computation...');
    % generate workers for parallel computation 
    %parpool();
    disp('Done.');
else
    disp('No workers for parallel computation have been enstablished');
end
% Run the Optimization 
disp('Running the Hybrid Optimization for the 1st trajectory (needed to apply model reduction meth.)...');
%{
% GENETIC ALGORITHM
%
% Custom otptions for the hybrid optimization
hybridopts = optimoptions('patternsearch');
hybridopts = optimoptions(hybridopts,'Display','iter',...
                                     'MaxIterations', 100000, ...
                                     'MaxFunctionEvaluations', 1e7, ...
                                     'AccelerateMesh', true, ...
                                     'PlotFcn', @psplotbestf);
% Start with the default options of GA
options = optimoptions('ga');
% Modify options setting
options = optimoptions(options,'Display', 'iter', ... %to visualize the temporary result at each iteration
                               'UseParallel', true, ...
                               'MaxGenerations',1000, ...
                               'PlotFcn', @gaplotbestf, ...
                               'PopulationSize', 100);
                     %%'HybridFcn', {@patternsearch, hybridopts}, ...
tic
[x0,fval0,exitflag0,output0,population0,score0] = ...
ga(@Trj,3*n*N,[],[],[],[],LB,UB,@Trj_con,options);
disp('Done.');
toc
disp(['The best cost fcn value is ',num2str(fval0)]);
%}
% FMINCON
%{
% Start with the default options
options = optimoptions('fmincon');
% Modify options setting
options = optimoptions(options,'Display', 'iter', ... %to visualize the temporary result at each iteration
                               'MaxFunEvals', 2e6, ...
                               'GradObj', 'off', ...
                               'UseParallel', true, ...
                               'TolFun', 1e-99, ...
                               'TolX', 1e-99); 
tic
[x0,fval0,exitflag0] = fmincon(@(Xi) Trj(Xi),pi/2*rand(1,3*n*N),[],[],[],[],LB,UB,[],options);
toc
%}
%
% HYBRID PARTICLE SWARM
hybridopts = optimoptions('patternsearch');
hybridopts = optimoptions(hybridopts,'Display','iter',...
                                     'MaxIterations', 100000, ...
                                     'MaxFunctionEvaluations', 1e7, ...
                                     'AccelerateMesh', true, ...
                                     'PlotFcn', @psplotbestf);
% Start with the default options
options = optimoptions('particleswarm');
% Modify options setting
%options = saoptimset('PlotFcns',{@saplotbestf,@saplotf});                    
options = optimoptions(options,'Display', 'iter', ... %to visualize the temporary result at each iteration
'HybridFcn', {@patternsearch, hybridopts}, ...                            
'PlotFcn', @pswplotbestf);%,...
                           %'SwarmSize', 200, ...    
                           %'UseParallel', true);
                           %'HybridFcn', {@patternsearch, hybridopts}, ...

tic
[x0,fval0] = ...
particleswarm(@(Xi) CostFcn(Xi,n,N),3*n*N,LB,UB,options);
disp('Done.');
toc
disp(['The best cost fcn value is ',num2str(fval0)]);

%}

%% Convert the result of the optimization process in a "readable" form for the other functions
 % (the matrices N*n: q, dq, ddq containing the optimal measurement pt.s for pos. , vel. and acc. respectively)
%T = v2m(x0,3*n);
T = reshape(x0,[N,3*n]);
T = [zeros(3*n,1),T',zeros(3*n,1)]; % add a column of zero as initial point of q, dq, ddq
Q = T';  
q0 = Q(:,1:n);
dq0 = Q(:,n+1:2*n);
ddq0 = Q(:,2*n+1:3*n);
t = 0:1:N;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% NUMERICAL SVD method for the model reduction
disp('Performing model Reduction...');
[gammaR,IC1,IC2,M] = ModelReduction(q0,dq0,ddq0,gamma,n,N+2);
CombineColumns;
disp('Done.');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% (Re)Calculate the optimal trajectory with the reduced model
disp('(Re)optimizing the trajectories for the reduced model...')
%delete(gcp('nocreate'));
%pool = parpool(c);
% HYBRID GENETIC ALGORITHM
%{
% Custom otptions for the hybrid optimization
hybridopts = optimoptions('patternsearch');
hybridopts = optimoptions(hybridopts,'Display','iter',...
                                     'UseParallel', always, ...
                                     'UseVectorized', false, ...
                                     'MaxIterations', 100000, ...
                                     'MaxFunctionEvaluations', 1e7, ...
                                     'AccelerateMesh', true, ...
                                     'PlotFcn', @psplotbestf);
% Start with the default options
options = optimoptions('ga');
% Modify options setting
options = optimoptions(options,'Display', 'iter', ... %to visualize the temporary result at each iteration
                               'HybridFcn', {@patternsearch, hybridopts}, ...
                               'UseParallel', true, ...
                               'PlotFcn', @gaplotbestf, ...
                               'PopulationSize', 200);
                     
tic
[x,fval,exitflag,output,population,score] = ...
ga(@Trj,3*n*N,[],[],[],[],LB,UB,[],[],options);
disp('Done.');
toc
%}
% HYBRID PARTICLE SWARM
hybridopts = optimoptions('patternsearch');
hybridopts = optimoptions(hybridopts,'Display','iter',...
                                     'MaxIterations', 100000, ...
                                     'MaxFunctionEvaluations', 1e7, ...
                                     'AccelerateMesh', true, ...
                                     'PlotFcn', @psplotbestf);
                                     %'UseParallel',true,...
                                     
% Start with the default options
options = optimoptions('particleswarm');
%options = saoptimset('PlotFcns',{@saplotbestf,@saplotf});                 
% Modify options setting
options = optimoptions(options,'Display', 'iter', ... %to visualize the temporary result at each iteration
'HybridFcn', {@patternsearch, hybridopts}, ...
'SwarmSize', 100, ...
'PlotFcn', @pswplotbestf);%,...                               
%'HybridFcn', {@patternsearch, hybridopts});%, ...
                            %'UseParallel', true);    
                           %'SwarmSize', 300, ...
                               %                            'PlotFcn', @pswplotbestf,...
                             
nq = randn(N,n);
nqd = randn(N,n);
nqdd = randn(N,n);
tic

[x,fval] = particleswarm(@(Xi) CostFcn2(Xi,n,N,nq,nqd,nqdd),3*n*N,LB,UB,options);
%[x,fval] = ...
%simulannealbnd(@(Xi) CostFcn(Xi,n,N),randn(1,3*n*N),LB,UB,options);
disp('Done.');
toc
disp(['The new best cost fcn value is ',num2str(fval)]);
T = reshape(x,[N,3*n]);
T = T';%[zeros(3*n,1) T' zeros(3*n,1)]; % add a column of zero as initial point of q, dq, ddq
Q = T';
q = Q(:,1:n);
qd = Q(:,n+1:2*n);
qdd = Q(:,2*n+1:3*n);
disp('Generation of the exciting trajectories completed.');
%
save('traj_opti.mat')
%