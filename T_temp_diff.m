fid = fopen('breastinfo_simple.txt','r');
breastID = str2double(fgetl(fid));
s1 = str2double(fgetl(fid));
s2 = str2double(fgetl(fid));
s3 = str2double(fgetl(fid));
class = str2double(fgetl(fid));
fclose(fid);

load mtype.mat;
load pval.mat;

muscle_wall = 153;
skin_start = 138;

% Convert vector into cube
mtype_cube = zeros(s1,s2,s3); % each voxel is .5mmx.5mmx.5mm
pval_cube = zeros(s1,s2,s3);
cur_pos = 1;
for k=1:s3
    for j=1:s2
        for i= 1:s1
            mtype_cube(i,j,k) = mtype(cur_pos);
            pval_cube(i,j,k) = pval(cur_pos);
            cur_pos = cur_pos + 1;
        end 
    end
end

% subsample cubes in order to solve sparse matrix
s1_ss = floor(s1/2); % voxels are now 1mmx1mmx1mm
s2_ss = floor(s2/2);
s3_ss = floor(s3/2);
xi = 1; yi = 1; zi = 1;
mtype_cube_subsamp = zeros(s1_ss,s2_ss,s3_ss);
pval_cube_subsamp = zeros(s1_ss,s2_ss,s3_ss);
for z=1:2:s3-1
    for y = 1:2:s2
        for x = 1:2:s1
            mid = mtype_cube(x,y,z);
            pid = pval_cube(x,y,z);
            mtype_cube_subsamp(xi,yi,zi) = mid;
            pval_cube_subsamp(xi,yi,zi) = pid;
            xi = xi+1;
        end
        xi = 1;
        yi = yi + 1;
    end
    yi = 1;
    zi = zi + 1;
end
% some voxels not converted to muscle during subsampling
% so do that now
% still need to figure out how to get pval converted for fdtd
for x=1:s1_ss
    for y=1:s2_ss
        for z=1:s3_ss
            if x > 153
               mtype_cube_subsamp(x,y,z) = -4;
                pval_cube_subsamp(x,y,z) = 1;
            end
        end
    end
end

% save mtype_cube_subsamp; save pval_cube_subsamp.mat;
% load mtype_cube_subsamp.mat;
[s1_ss, s2_ss, s3_ss] = size(mtype_cube_subsamp);
model = mtype_cube_subsamp; air_id = -1;

% file_name = 'scat_fib_1.csv';
% create_csv(file_name,model,s1_ss,s2_ss,s3_ss,air_id,mtype_cube_subsamp,pval_cube_subsamp);

figure; colormap(gray);
contourf(model(:,:,floor(s3_ss/2))); 
% model_pval = pval_cube_subsamp;
figure; contourf((0:s2_ss-1)*.1,(0:s1_ss-1)*.1,model(:,:,floor(s3_ss/2)),'LineStyle','none');
% colormap(flipud(gray));brighten(.4);
xlabel('Distance (cm)','FontSize',14); ylabel('Distance (cm)','FontSize',14);
hcb = colorbar; 
set(hcb,'YTick',[0,1,2,3,4,5,6,7,8],'YTickLabel',{'Air','Skin','Gland-1','Gland-2','Gland-3','Fat-1','Fat-2','Fat-3','Muscle'},'FontSize',14);

% create temperature model
tumor_on = 0; tumor_depth = 8;  tumor_radius = 10; Tambient = 27; Tart = 37; 
% muscle_wall = 154; skin_start = 0; % can't remember why I put this in
[T_3d_nom,tissue_3d_nom] = gen_breast_therm_model(model,s1_ss,s2_ss,s3_ss,tumor_on,tumor_depth,tumor_radius,Tambient,Tart,muscle_wall,skin_start);
figure; contourf(T_3d_nom(:,:,floor(s3_ss/2)));
colorbar;
xlabel('Distance (mm)','FontSize',14); ylabel('Distance (mm)','FontSize',14);
title('Normal Temperature Profile (\circC)','FontSize',14);

% Generate temperature anomalies with radius = 10
tumor_on = 1; tum_y_cen = 90; tum_z_cen = floor(s3_ss/2);

tumor_radius = 10; tumor_depth = 10; tum_x_cen = 20 + tumor_depth;% tumor dept is 1cm
[T_3d_abn1,tissue_3d_abn1] = gen_breast_therm_model(model,s1_ss,s2_ss,s3_ss,tumor_on,tumor_depth,tumor_radius,Tambient,Tart,muscle_wall,skin_start,tum_x_cen,tum_y_cen,tum_z_cen);

T_3d_diff = T_3d_abn1 - T_3d_nom;
figure
plot(T_3d_diff(20:100,90,80),'k','LineWidth',1);



