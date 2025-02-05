% Copyright (c) 2025 Matteo Giordano
%
% Codes accompanying the article "Bayesian inference for initial heat states 
% with Gaussian series priors" by Matteo Giordano
%
%%
% Discrete noisy observations of heat equation solution
%
% Let O be a smooth domain in R^2. Consider the 2D heat equation with 
% homogeneous Dirichlet boundary conditions:
%
%   d_tu(t,x) - div(c grad u)(t,x) = 0, in (0,T]xO
%   u(t,x)=0, on (0,T]xboundary of O
%   u(0,x) = f(x), in O
%   
%
% where f:O -> R is the unknown (sufficiently smooth) initial condition and 
% c:O -> (0,+infty) is a given smooth diffusivity.
%
% There exists a unique classical solution G(f)=u_f, that depends linearly on
% f. We observe data
%
%   Y_i = u_f(T,x_i) + sigma W_i,   i = 1,...,n
%
% where x_i are design points in O, sigma>0 and W_i are i.i.d. 
% N(0,1) random variables.
%
% The following code generates n observations Y_1,...,Y_n.

%%
% Create rotate elliptically-shaped domain and triangular mesh

% Display more digits
format long

ax_h = 1; 
    % length of horizontal semiaxis
ax_v = .75; 
    % length of vertical semiaxis
rot = pi/6;
    % angle of rotation
t = linspace(0,2*pi,1000);
pgon = polyshape({ax_h*cos(t)*cos(rot) - ax_v*sin(t)*sin(rot)},...
    {ax_v*sin(t)*cos(rot) + ax_h*cos(t)*sin(rot)});
vol = pi*ax_h*ax_v;

% Create a triangulation representation of pgon
tr = triangulation(pgon);

% Create a PDE model
model = createpde;

% With the triangulation data as a mesh, use the geometryFromMesh function
% to create a geometry
tnodes = tr.Points';
telements = tr.ConnectivityList';
geometryFromMesh(model,tnodes,telements);
%pdegplot(model)

% Generate and plot triangular mesh
generateMesh(model,'Hmax',.05);
%figure()
%axes('FontSize', 20, 'NextPlot','add')
%pdemesh(model)
%xticks([-1,-.5,0,.5,1])
%yticks([-1,-.5,0,.5,1])
%xlabel('x', 'FontSize', 20);
%ylabel('y', 'FontSize', 20);
mesh_nodes = model.Mesh.Nodes; 
    % 2 x mesh_size matrix whose columns contain the (x,y) coordinates 
    % of the nodes in the mesh
mesh_nodes_num = size(mesh_nodes); 
mesh_nodes_num=mesh_nodes_num(2); 
    % number of nodes in the mesh
mesh_elements = model.Mesh.Elements; 
    % 6 x mesh_elements_num whose columns contain the 6 node indices 
    % identifying each triangle. The first 3 elements of each column contain 
    % the indices of the 3 vertices of the triangle 
mesh_elements_num = size(mesh_elements); 
mesh_elements_num = mesh_elements_num(2); 
    % number of triangles in the mesh
[~,mesh_elements_area] = area(model.Mesh); 

% Compute barycenters of triangular mesh elements
barycenters = zeros(2,mesh_elements_num);
for i=1:mesh_elements_num
    barycenters(:,i) = mean(mesh_nodes(:,mesh_elements(1:3,i)),2);
end


%%
% Specify diffusivity c and true initial condition f_0

% Specify the diffusivity c as a function of (x,y)
c_min = 2.5;
c = @(x,y) c_min - exp(-(5*x-2).^2-(2.5*y-.5).^2);
c_fun=@(location,state) c(location.x,location.y);
    % function of (location,state) to pass to PDE solver
c_mesh = c(mesh_nodes(1,:),mesh_nodes(2,:));
    % evaluates c at mesh
%c_bary = c(barycenters(1,:),barycenters(2,:));
    % evaluates c at barycenters

% Plot c
figure()
axes('FontSize', 20, 'NextPlot','add')
pdeplot(model,'XYData',c_mesh)
colorbar('Fontsize',20)
% pdeplot(model,'XYData',f0_num,'ZData',f0_num,'ColorMap',jet) 
    % 3D plot
title('Known diffusivity c','FontSize',20);
xticks([-1,-.5,0,.5,1])
yticks([-1,-.5,0,.5,1])
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
crameri vik

% Specify the unknown source f_0 as a function of (x,y)
f0 = @(x,y) 5*exp(-(5*x-2).^2-(5*y-1).^2);
f0_mesh = f0(mesh_nodes(1,:),mesh_nodes(2,:));
    % evaluates f_0 at mesh
f0_bary=f0(barycenters(1,:),barycenters(2,:));
    % evaluates f_0 at barycenters
f0_norm=sqrt(sum((f0_bary).^2.*mesh_elements_area));
    % approximates the L^2 norm of f_0

% Plot f0
figure()
axes('FontSize', 20, 'NextPlot','add')
pdeplot(model,'XYData',f0_mesh,'ColorMap',jet)
colorbar('Fontsize',20)
title('True initial condition f_0','FontSize',20);
xticks([-1,-.5,0,.5,1])
yticks([-1,-.5,0,.5,1])
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
crameri -lajolla

%%
% Find eigenbasis of the associated infinitesimal generator

% Specify zero Dirichlet boundary conditions on all edges
%applyBoundaryCondition(model,'dirichlet','Edge',1:model.Geometry.NumEdges,'u',0);

% Specify zero Neumann boundary conditions on all edges
applyBoundaryCondition(model,'neumann','Edge',1:model.Geometry.NumEdges,...
    'g',0,'q',0); 

% Specify coefficients for eigenvalue equation associated to infinitesimal
% generator
specifyCoefficients(model,'m',0,'d',1,'c',c_fun,'a',0,'f',0);

% Specify range of searc of eigenvalues
range = [-1,750]; 

%c_fun=@(location,state) 10*c(location.x,location.y);

% Solve eigenvalue problem
results = solvepdeeig(model,range); 
lambdas_SVD = results.Eigenvalues; 
J_SVD = length(lambdas_SVD); 
e_SVD_mesh = results.Eigenvectors;
e_SVD_mesh(:,1:J_SVD) = e_SVD_mesh(:,1:J_SVD)*sqrt(mesh_nodes_num/vol);
    % normalisation of eigenfunctions to have unit area
%disp(['lambda_SVD_1= ', num2str(lambdas_SVD(1))])

% Plot the eigenfunctions
%figure()
%subplot(1,3,1)
%pdeplot(model,'XYData',e_SVD_mesh(:,1));
%title('e_0','FontSize',20)
%subplot(1,3,2)
%pdeplot(model,'XYData',e_SVD_mesh(:,2));
%title('e_2','FontSize',20)
%subplot(1,3,3)
%pdeplot(model,'XYData',e_SVD_mesh(:,J_SVD));
%title('e_J','FontSize',20)

% Plot the eigenvalues
%figure()
%axes('FontSize', 20, 'NextPlot','add')
%plot(lambdas_SVD,'.','Linewidth',3)
%xlabel('j', 'FontSize', 20);
%ylabel('\lambda_j', 'FontSize', 20);

% Evaluate eigenfunctions at barycenters
e_SVD_bary = zeros(mesh_elements_num,J_SVD);
for j=1:J_SVD
    ej_SVD_interp=scatteredInterpolant(mesh_nodes(1,:)', ...
        mesh_nodes(2,:)',e_SVD_mesh(:,j));
    e_SVD_bary(:,j)=ej_SVD_interp(barycenters(1,:),barycenters(2,:));
end

%%
% Projection of true initial condition f_0 onto eigenbasis for computation
% of heat equation solution via spectral characterisation

% Compute Fourier coefficients of f_0 in the eigenbasis
f0_SVD_coeff=zeros(J_SVD,1); 
for j=1:J_SVD
    f0_SVD_coeff(j)=sum(mesh_elements_area.*f0_bary.*e_SVD_bary(:,j)');
end

% Compute projection
%f0_proj_mesh=zeros(1,mesh_nodes_num);
%f0_proj_bary=zeros(1,mesh_elements_num);
%for j=1:J_SVD
%    f0_proj_mesh = f0_proj_mesh+f0_SVD_coeff(j)*e_SVD_mesh(:,j)';
%    f0_proj_bary = f0_proj_bary+f0_SVD_coeff(j)*e_SVD_bary(:,j)';
%end

% Compare f_0 to its projection and the approximation error
%figure()
%subplot(1,3,1)
%pdeplot(model,'XYData',f0_mesh,'ColorMap',jet)
%title('True f_0','FontSize',20)
%subplot(1,3,2)
%pdeplot(model,'XYData',f0_proj_mesh,'ColorMap',jet)
%title('Projection of f_0','FontSize',20)
%subplot(1,3,3)
%pdeplot(model,'XYData',f0_mesh-f0_proj_mesh,'ColorMap',jet)
%title('Approximation error','FontSize',20)

% Approximate L^2 distance between f_0 and posterior mean
%approx_error = sqrt(sum((f0_bary-f0_proj_bary).^2.*mesh_elements_area));
%disp(['L^2 approximation error via projection = ', num2str(approx_error)])
%disp(['L^2 norm of f_0 = ', num2str(f0_norm)])
%disp(['Relative error = ', num2str(round(approx_error/f0_norm,4))])

%%
% Heat equation solution corresponding to diffusivity c and true initial
% condition f_0

% Specify time horizon
T = .01;
%times=linspace(0,T,5);
times=[T];

u_f0 = zeros(mesh_nodes_num,length(times));

for t=1:length(times)
    for j=1:J_SVD
        u_f0(:,t) = u_f0(:,t) + exp(-lambdas_SVD(j)*times(t))*...
            f0_SVD_coeff(j)*e_SVD_mesh(:,j);
    end
%    figure()
%    axes('FontSize', 20, 'NextPlot','add')
%    pdeplot(model,'XYData',u_f0(:,t),'ColorMap',jet)
%    xticks([-1,-.5,0,.5,1])
%    yticks([-1,-.5,0,.5,1])
%    xlabel('x', 'FontSize', 20);
%    ylabel('y', 'FontSize', 20);
%    clim([min(f0_mesh) max(f0_mesh)]);
%    crameri -lajolla
end

% Solution at time T
u0=u_f0(:,length(times));

% Norm of solution
u0_interp=scatteredInterpolant(mesh_nodes(1,:)',mesh_nodes(2,:)',u0);
u0_bary=u0_interp(barycenters(1,:),barycenters(2,:));
u0_norm = sqrt(sum((u0_bary).^2.*mesh_elements_area));
%disp(['Norm of u0 = ', num2str(u0_norm)])

%%
% Noisy observations of PDE solution

rng(1)

% Sample design points and noise variance
sample_size=1000;
if sample_size>mesh_nodes_num
    sample_size=mesh_nodes_num;
end
sigma=0.05;
disp(['Signal to noise ratio = ', num2str(u0_norm/sigma)])
    % noise standard deviation
rand_index=sort(randsample(mesh_nodes_num,sample_size)); 
    % random indices in the mesh
rand_mesh=mesh_nodes(:,rand_index); 
    % random sample of mesh points drawn uniformly at random
%figure()
%axes('FontSize', 20, 'NextPlot','add')
%scatter(rand_mesh(1,:),rand_mesh(2,:),'filled') 
    % plot of sampled locations
%title('Design points X_1,...,X_n','FontSize', 20)
%xticks([-1,-.5,0,.5,1])
%yticks([-1,-.5,0,.5,1])
%xlabel('x','FontSize', 20)
%ylabel('y', 'FontSize', 20)

% Sample observations
observations=u0(rand_index)+(mvnrnd(zeros(sample_size,1),...
    sigma^2*eye(sample_size)))'; 
    % add i.i.d N(0,sigma^2) noise to the observation
u0_noisy=u0;
u0_noisy(rand_index)=observations;

% Plot corrupted PDE solution
figure()
axes('FontSize', 20, 'NextPlot','add')
pdeplot(model,'XYData',u0_noisy)
%title('Observations Y_i=u_{f_0}(T,X_i)+\sigma W_i','FontSize', 20);
xlabel('x','FontSize', 20)
ylabel('y', 'FontSize', 20)
xticks([-1,-.5,0,.5,1])
yticks([-1,-.5,0,.5,1])
clim([min(f0_mesh) max(f0_mesh)]);
crameri -lajolla
