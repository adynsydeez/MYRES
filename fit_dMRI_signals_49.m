%% loading trace data
clear;
PATH   = 'C:\Users\yangq\OneDrive - Queensland University of Technology\MGH Connectome';
SUBJECT    = 'sub_001';
datatype = '_real'; % '_real' or empty string ''
FOLDER = 'C:\Users\adyna\Downloads\dwi_real\dwi_real\'; % [PATH '\' SUBJECT '\dwi' datatype '\'];

delta = 8; % in units of ms
gyroRatio = 267.522E6; % in units of rad/s/T
G = [0 31, 68, 105, 142, 179, 216, 253, 290]' * 1E-3; % in units of T/m
qvec = gyroRatio*G*(delta*1e-3)*1e-3; % in units of 1/mm

% brain mask 
brainmask = double(niftiread([FOLDER SUBJECT '_dwi' datatype '_brainmask.nii.gz']));

S0 = double(niftiread([FOLDER SUBJECT '_dwi' datatype '_b0_delta0_image.nii'])).*brainmask;
%figure, imagesc(squeeze(S0(:,:,73))); axis image; view(-90,90); axis off;

num_nz = nnz(S0)

Delta = 49;
diff_time = (Delta - delta/3)*1E-3; % in units of s

% dwi with Delta = 49 ms
bvec49 = [0 200 950 2300 4250 6750 9850 13500 17800]';

S9  = double(niftiread([FOLDER SUBJECT '_dwi' datatype '_b200_delta49_image.nii'])).*brainmask;
S10 = double(niftiread([FOLDER SUBJECT '_dwi' datatype '_b950_delta49_image.nii'])).*brainmask;
S11 = double(niftiread([FOLDER SUBJECT '_dwi' datatype '_b2300_delta49_image.nii'])).*brainmask;
S12 = double(niftiread([FOLDER SUBJECT '_dwi' datatype '_b4250_delta49_image.nii'])).*brainmask;
S13 = double(niftiread([FOLDER SUBJECT '_dwi' datatype '_b6750_delta49_image.nii'])).*brainmask;
S14 = double(niftiread([FOLDER SUBJECT '_dwi' datatype '_b9850_delta49_image.nii'])).*brainmask;
S15 = double(niftiread([FOLDER SUBJECT '_dwi' datatype '_b13500_delta49_image.nii'])).*brainmask;
S16 = double(niftiread([FOLDER SUBJECT '_dwi' datatype '_b17800_delta49_image.nii'])).*brainmask;

figure
tiledlayout(1,9,'TileSpacing','tight','Padding','compact')
nexttile
imagesc(squeeze(S0(:,5:95,73))); axis image; view(-90,90); axis off; title('b = 0 s/mm^2', FontSize=14)
nexttile
imagesc(squeeze(S9(:,5:95,73))); axis image; view(-90,90); axis off; title('b = 200 s/mm^2', FontSize=14)
nexttile
imagesc(squeeze(S10(:,5:95,73))); axis image; view(-90,90); axis off; title('b = 950 s/mm^2', FontSize=14)
nexttile
imagesc(squeeze(S11(:,5:95,73))); axis image; view(-90,90); axis off; title('b = 2300 s/mm^2', FontSize=14)
nexttile
imagesc(squeeze(S12(:,5:95,73))); axis image; view(-90,90); axis off; title('b = 4250 s/mm^2', FontSize=14)
nexttile
imagesc(squeeze(S13(:,5:95,73))); axis image; view(-90,90); axis off; title('b = 6750 s/mm^2', FontSize=14)
nexttile
imagesc(squeeze(S14(:,5:95,73))); axis image; view(-90,90); axis off; title('b = 9850 s/mm^2', FontSize=14)
nexttile
imagesc(squeeze(S15(:,5:95,73))); axis image; view(-90,90); axis off; title('b = 13500 s/mm^2', FontSize=14)
nexttile
imagesc(squeeze(S16(:,5:95,73))); axis image; view(-90,90); axis off; title('b = 17800 s/mm^2', FontSize=14)
colormap("gray")


%% fitting model parameters

% mono-exponential model: D
% DKI model: D and K
% Sub-diffusion model: D and beta, computing D* and K* from D and beta;
% lsqcurvefit

[X,Y,L] = size(S0);
D_MONO = nan(X,Y,L);
error_MONO = nan(X,Y,L);

D_guess = 1e-3; Dmin = 0; Dmax = 3e-3; % mm^2/s
K_guess = 0; Kmin = 0; Kmax = 3;
beta_guess = 0.5; betamin = 0; betamax = 1; % time fractional index

% fit each voxel 
tic
for l= 73 % select one axial slice
    for x = 1:X
        parfor y = 1:Y
            %Skip voxels outside for brain mask
            if S0(x,y,l) <= 0
                D_MONO(x,y,l) = 0;
                error_MONO(x,y,l) = 0;
                continue; 
            end
            
            Svec = zeros(9,1);
            Svec(1) = S0(x,y,l);
            Svec(2) = S9(x,y,l);
            Svec(3) = S10(x,y,l);
            Svec(4) = S11(x,y,l);
            Svec(5) = S12(x,y,l);
            Svec(6) = S13(x,y,l);
            Svec(7) = S14(x,y,l);
            Svec(8) = S15(x,y,l);
            Svec(9) = S16(x,y,l);

            Svec = Svec/S0(x,y,l);

            [D_MONO(x,y,l), error_MONO(x,y,l)] = fit_voxel_mono(Svec(1:3), bvec49(1:3), Dmin, Dmax);
            [D_DKI(x,y,l), K_DKI(x,y,l), error_DKI(x,y,l)] = fit_voxel_dki(Svec(1:4), bvec49(1:4), D_guess, K_guess, Dmin, Dmax, Kmin, Kmax);
            [D_SUB(x,y,l), BETA_SUB(x,y,l), error_SUB(x,y,l)] = fit_voxel_sub_ml(Svec, bvec49, D_guess, beta_guess, Dmin, Dmax, betamin, betamax);

        end
    end
end
toc

%% Parameter maps
l = 73; % select an axial slice

figure; t = tiledlayout(3,4,'TileSpacing','tight','Padding','compact');
title(t,'Parameter maps','Color',[1, 1 ,1])
set(gcf, 'InvertHardCopy', 'off'); set(gcf,'color', [0 0 0]);

nexttile
imagesc(squeeze(D_MONO(:,:,l))), axis image, title('D_{mono}','Color',[1, 1 ,1]);clim([0 3e-3]);
colorbar('Color',[1, 1 ,1]); view(-90,90); axis off; colormap hot

nexttile
imagesc(squeeze(error_MONO(:,:,l))), axis image,clim([0 0.1]); title('RMSE','Color',[1, 1 ,1]);
colorbar('Color',[1, 1 ,1]); view(-90,90); axis off;

nexttile
nexttile

nexttile
imagesc(squeeze(D_DKI(:,:,l))), axis image, title('D_{DKI}','Color',[1,1,1]); clim([0 3e-3]);
colorbar('Color',[1,1,1]); view(-90,90); axis off; colormap hot

nexttile
imagesc(squeeze(error_DKI(:,:,l))), axis image, clim([0 0.1]); title('RMSE','Color',[1,1,1]);
colorbar('Color',[1,1,1]); view(-90,90); axis off;

nexttile
imagesc(squeeze(K_DKI(:,:,l))), axis image, clim([0 3]); title('K','Color',[1,1,1]);
colorbar('Color',[1,1,1]); view(-90,90); axis off;

nexttile

nexttile
imagesc(squeeze(D_SUB(:,:,l))), axis image, title('D_{SUB}','Color',[1,1,1]); clim([0 3e-3]);
colorbar('Color',[1,1,1]); view(-90,90); axis off; colormap hot

nexttile
imagesc(squeeze(error_SUB(:,:,l))), axis image, clim([0 0.1]); title('RMSE','Color',[1,1,1]);
colorbar('Color',[1,1,1]); view(-90,90); axis off;

nexttile
imagesc(squeeze(BETA_SUB(:,:,l))), axis image, clim([0 1]); title('Beta','Color',[1,1,1]);
colorbar('Color',[1,1,1]); view(-90,90); axis off;

k_star = 6 .* gamma(1+BETA_SUB).^2./gamma(1+2*BETA_SUB) - 3;

nexttile
imagesc(squeeze(k_star(:,:,l))), axis image, clim([0 3]); title('K_star','Color',[1,1,1]);
colorbar('Color',[1,1,1]); view(-90,90); axis off;


%% pick a voxel to check the fitting
x = 43; y = 55; l = 73;
% plot S against b or q with data point and fitted curve
