% Copyright (c) 2025 Matteo Giordano
%
% Codes accompanying the article "Bayesian inference for initial heat states 
% with Gaussian series priors" by Matteo Giordano
%
%%
% Bayesian nonparametric inference for the initial condition f:O -> R 
% with the Gaussian series priors on the Dirichlet Laplacian eigebasis
% via conjugate formulae

% Requires output of GenerateObservations.m (including f_0, observations 
% and geometry) and K_mat.m (Matérn covariance kernel).

%%
% Mesh for computation of the Dirichlet-Laplacian eigenpairs

model_prior = createpde(); 
geometryFromMesh(model_prior,tnodes,telements);
generateMesh(model_prior,'Hmax',0.075);
mesh_nodes_prior=model_prior.Mesh.Nodes;
mesh_nodes_num_prior=size(mesh_nodes_prior); 
mesh_nodes_num_prior=mesh_nodes_num_prior(2);
mesh_elements_prior = model_prior.Mesh.Elements;
mesh_elements_num_prior = size(mesh_elements_prior); 
mesh_elements_num_prior = mesh_elements_num_prior(2); 
[~,mesh_elements_area_prior] = area(model_prior.Mesh); 

% Compute barycenters of triangular mesh elements
barycenters_prior = zeros(2,mesh_elements_num_prior);
for i=1:mesh_elements_num_prior
    barycenters_prior(:,i) = mean(mesh_nodes_prior(:,mesh_elements_prior(1:3,i)),2);
end

%%
% Solve elliptic eigenvalue problem for the Dirichlet-Laplacian

tic

% Specity homogeneous Dirichlet boundary conditions
applyBoundaryCondition(model_prior,'dirichlet','Edge', ...
    1:model.Geometry.NumEdges,'u',0); 
% Specify coefficients for eigenvalue equation
specifyCoefficients(model_prior,'m',0,'d',1,'c',1,'a',0,'f',0);
range = [-1,500]; 
    % range of search for eigenvalues
results = solvepdeeig(model_prior,range); 
    % solve eigenvalue equation
lambdas_basis = results.Eigenvalues; 
    % extract eigenvalues
J_basis = length(lambdas_basis); 
    % number of eigenvalues (dimension of discretised parameter space)
e_basis = results.Eigenvectors; 
    % extract eigenfunctions
e_basis(:,1:J_basis) = e_basis(:,1:J_basis)*sqrt(mesh_nodes_num_prior/vol);
    % normalise eigenfunctions

toc

figure() 
subplot(1,3,1)
pdeplot(model_prior,'XYData',e_basis(:,1)); 
    % plot first eigenfunction
title('e_0','FontSize',20)
subplot(1,3,2)
pdeplot(model_prior,'XYData',e_basis(:,2)); 
    % plot second eigenfunction
title('e_2','FontSize',20)
subplot(1,3,3)
pdeplot(model_prior,'XYData',e_basis(:,J_basis)); 
% plot eigenfunction corresponding to the largest found eigenvalue
title('e_J','FontSize',20)

% Plot the eigenvalues
figure()
axes('FontSize', 20, 'NextPlot','add')
plot(lambdas_basis,'*','Linewidth',1)
xlabel('j', 'FontSize', 20);
ylabel('\lambda_j', 'FontSize', 20);
%legend('\lambda_j=O(j)','FontSize',25)

%%
% Projection of f_0 onto the Dirichlet-Laplacian eigenbasis

f0_mesh_prior = f0(mesh_nodes_prior(1,:),mesh_nodes_prior(2,:));
f0_interp=scatteredInterpolant(mesh_nodes_prior(1,:)',mesh_nodes_prior(2,:)',f0_mesh_prior');
f0_bary_prior=f0_interp(barycenters_prior(1,:),barycenters_prior(2,:));

f0_coeff=zeros(J_basis,1); 
    % initialises vector to store the Fourier coefficients of F0 in the 
    % Dirichlet-Laplacian eigenbasis

for j=1:J_basis
    ej_basis_interp=scatteredInterpolant(mesh_nodes_prior(1,:)',...
        mesh_nodes_prior(2,:)',e_basis(:,j));
    ej_basis_bary=ej_basis_interp(barycenters_prior(1,:),...
        barycenters_prior(2,:));
    f0_coeff(j)=sum(mesh_elements_area_prior.*f0_bary_prior.*ej_basis_bary);
end

f0_proj=zeros(1,mesh_nodes_num_prior);
for j=1:J_basis
    f0_proj = f0_proj+f0_coeff(j)*e_basis(:,j)';
end

figure()
axes('FontSize', 20, 'NextPlot','add')
subplot(1,3,1)
pdeplot(model_prior,'XYData',f0_mesh_prior,'ColorMap',jet)
title('True f_0','FontSize',20)
%colorbar('FontSize',20)
%xlabel('x','FontSize', 20)
%ylabel('y', 'FontSize', 20)
xticks([-1,-.5,0,.5,1])
yticks([-1,-.5,0,.5,1])

subplot(1,3,2)
pdeplot(model_prior,'XYData',f0_proj,'ColorMap',jet)
title('Projection of f_0','FontSize',20)
%colorbar('FontSize',20)
%xlabel('x','FontSize', 20)
%ylabel('y', 'FontSize', 20)
xticks([-1,-.5,0,.5,1])
yticks([-1,-.5,0,.5,1])

subplot(1,3,3)
pdeplot(model_prior,'XYData',f0_mesh_prior-f0_proj,'ColorMap',jet)
title('Approximation error','FontSize',20)
%colorbar('FontSize',20)
%xlabel('x','FontSize', 20)
%ylabel('y', 'FontSize', 20)
xticks([-1,-.5,0,.5,1])
yticks([-1,-.5,0,.5,1])

% Compute piecewise constant approximations of the projection of f_0 at
% the triangle baricenters
f0_proj_bary_prior = griddata(mesh_nodes_prior(1,:),mesh_nodes_prior(2,:),f0_proj,...
    barycenters_prior(1,:),barycenters_prior(2,:));

% Approximate L^2 distance between f_0 and posterior mean
approx_error = sqrt(sum((f0_bary_prior-f0_proj_bary_prior).^2.*mesh_elements_area_prior));
disp(['L^2 approximation error via projection = ', num2str(approx_error)])
% L^2 norm of f0
f0_norm=norm(f0_coeff);
%disp(['L^2 norm of f_0 = ', num2str(f0_norm)])
disp(['Relative error = ', num2str(round(approx_error/f0_norm,4))])

%%
% Prior covariance matrix for Gaussian series prior draws

prior_regularity=.5; 
prior_cov=diag(lambdas_basis.^(-2*prior_regularity)); 
    % diagonal prior covariance matrix
prior_cov_inv=diag(lambdas_basis.^(2*prior_regularity)); 
    % inverse of prior covariance matrix

%%
% Sample and plot a prior draw

rng(1)

theta_rand=mvnrnd(zeros(J_basis,1),prior_cov,1)'; 
% sample Fourier coefficients from prior

f_rand=zeros(1,mesh_nodes_num_prior);

for j=1:J_basis
    f_rand = f_rand+theta_rand(j)*e_basis(:,j)';
end

figure()
axes('FontSize', 20, 'NextPlot','add')
pdeplot(model_prior,'XYData',f_rand,'ColorMap',jet);
%title('f\sim\Pi(\cdot)','FontSize',20)
colorbar('FontSize',20)
xlabel('x','FontSize', 20)
ylabel('y', 'FontSize', 20)
xticks([-1,-.5,0,.5,1])
yticks([-1,-.5,0,.5,1])
%clim([min(f0_mesh_prior),max(f0_mesh_prior)])
crameri -lajolla

%%
% Discretisation of solution operator

% Discretisation of operator for the elements of basis over full mesh
fwd_operator=zeros(sample_size,J_basis); 
    % mesh_size x coeff_num matrix, whose (i,j)-th element is 
    % L^{-1}e_j((x_i,y_i)) where L^{-1} is the solution operator, e_j is
    % the j-th basis function and (x_i,y_i) is the i-th design point.

tic

for j=1:J_basis
    ej_bary= griddata(mesh_nodes_prior(1,:),mesh_nodes_prior(2,:),e_basis(:,j),...
        barycenters(1,:),barycenters(2,:));
    ej_SVD_coeff=zeros(J_SVD,1); 
    for k=1:J_SVD
        ej_SVD_coeff(k)=sum(mesh_elements_area.*ej_bary.*e_SVD_bary(:,k)');
    end
    u_ej = zeros(mesh_nodes_num,1);
    for k=1:J_SVD
        u_ej = u_ej + exp(-lambdas_SVD(k)*T)*ej_SVD_coeff(k)*e_SVD_mesh(:,k);
    end
    %figure()
    %axes('FontSize', 20, 'NextPlot','add')
    %pdeplot(model,'XYData',u_phim,'ColorMap',jet)
    %xticks([-1,-.5,0,.5,1])
    %yticks([-1,-.5,0,.5,1])
    %xlabel('x', 'FontSize', 20);
    %ylabel('y', 'FontSize', 20);
    %crameri -lajolla
    fwd_operator(:,j) = u_ej(rand_index);
end

toc

%%
% Conjugate formulae for posterior variance and mean

% Posterior covariance matrix
post_cov=inv(fwd_operator'*fwd_operator/(sigma^2)+prior_cov_inv);
% Posterior mean
post_mean = post_cov*fwd_operator'*observations/(sigma^2);

f_mean_mesh_prior=zeros(mesh_nodes_num_prior,1);
for j=1:J_basis
    f_mean_mesh_prior = f_mean_mesh_prior+post_mean(j)*e_basis(:,j);
end

% Plot posterior mean and estimation erorr
figure()
%subplot(1,2,1)
axes('FontSize', 20, 'NextPlot','add')
pdeplot(model_prior,'XYData',f0_mesh_prior,'ColorMap',jet)
title('True f_0','FontSize',20)
clim([min(f0_mesh_prior),max(f0_mesh_prior)])
colorbar('Fontsize',20)
xticks([-1,-.5,0,.5,1])
yticks([-1,-.5,0,.5,1])
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
crameri -lajolla

figure()
%subplot(1,2,2)
axes('FontSize', 20, 'NextPlot','add')
pdeplot(model_prior,'XYData',f_mean_mesh_prior,'ColorMap',jet)
colorbar('Fontsize',20)
title('Posterior mean estimate','FontSize',20)
clim([min(f0_mesh_prior),max(f0_mesh_prior)])
clim([min(f0_mesh_prior),max(f0_mesh_prior)])
xticks([-1,-.5,0,.5,1])
yticks([-1,-.5,0,.5,1])
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
crameri -lajolla

figure()
axes('FontSize', 20, 'NextPlot','add')
pdeplot(model_prior,'XYData',f0_mesh_prior-f_mean_mesh_prior','ColorMap',hot)
title('Estimation error','FontSize',20)
colorbar('Fontsize',20)
xticks([-1,-.5,0,.5,1])
yticks([-1,-.5,0,.5,1])
xlabel('x', 'FontSize', 20);
ylabel('y', 'FontSize', 20);
crameri -lajolla

% Approximate L^2 distance between f_0 and posterior mean
estim_error = norm(post_mean-f0_coeff);
disp(['L^2 estimation error = ', num2str(estim_error)])
%disp(['L^2 norm of f_0 = ', num2str(f0_norm)])
disp(['Relative error = ', num2str(round(estim_error/f0_norm,4))])