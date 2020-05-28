clear all;
rr_out = zeros(500,1);
N_tm_out = zeros(500,1);
N_t1_out = zeros(500,1);
N_t2_out = zeros(500,1);
N_a1_out = zeros(500,1);
N_chi_out = zeros(500,1);
F_tm_out = zeros(500,1);
F_t1_out = zeros(500,1);
F_t2_out = zeros(500,1);
F_a1_out = zeros(500,1);
F_chi_out = zeros(500,1);
Cellpix_out = zeros(500,1);
Cytopix_out = zeros(500,1);
Cellnum_out = zeros(500,1);
Imnum_out = zeros(500,1);
D_out = zeros(500,1);
for a = 1:7;
    try
    mask_image = imread(strcat('E:\020718\Dish 1\Mask_image_0',num2str(a+0),'.tif'));
    catch
        continue
    end
    mask_image = imread(strcat('E:\020718\Dish 1\Mask_image_0',num2str(a+0),'.tif'));
    cell_image = imread(strcat('E:\020718\Dish 1\Cell_mask_0',num2str(a+0),'.tif'));
    N_photons = imread(strcat('E:\020718\Dish 1\Tcells-00',num2str(a),'_n_photons.tiff'));
    F_photons = imread(strcat('E:\020718\Dish 1\Tcells-00',num2str(a),'_f_photons.tiff'));
    N_t1 = imread(strcat('E:\020718\Dish 1\Tcells-00',num2str(a),'_n_t1.tiff'));
    N_t2 = imread(strcat('E:\020718\Dish 1\Tcells-00',num2str(a),'_n_t2.tiff'));
    N_a1 = imread(strcat('E:\020718\Dish 1\Tcells-00',num2str(a),'_n_a1[%].tiff'));
    N_a2 = imread(strcat('E:\020718\Dish 1\Tcells-00',num2str(a),'_n_a2[%].tiff'));
    N_chi = imread(strcat('E:\020718\Dish 1\Tcells-00',num2str(a),'_n_chi.tiff'));
    F_t1 = imread(strcat('E:\020718\Dish 1\Tcells-00',num2str(a),'_f_t1.tiff'));
    F_t2 = imread(strcat('E:\020718\Dish 1\Tcells-00',num2str(a),'_f_t2.tiff'));
    F_a1 = imread(strcat('E:\020718\Dish 1\Tcells-00',num2str(a),'_f_a1[%].tiff'));
    F_a2 = imread(strcat('E:\020718\Dish 1\Tcells-00',num2str(a),'_f_a2[%].tiff'));
    F_chi = imread(strcat('E:\020718\Dish 1\Tcells-00',num2str(a),'_f_chi.tiff'));
    rr_image = N_photons./(N_photons+F_photons);
    N_tm_image = N_t1.*N_a1./100+N_t2.*N_a2./100;
    F_tm_image = F_t1.*F_a1./100+F_t2.*F_a2./100;
    
    rr_cell = zeros(1,max(max(mask_image)));
    N_tm_cell = zeros(1,max(max(mask_image)));
    N_t1_cell = zeros(1,max(max(mask_image)));
    N_t2_cell = zeros(1,max(max(mask_image)));
    N_a1_cell = zeros(1,max(max(mask_image)));
    F_tm_cell = zeros(1,max(max(mask_image)));
    F_t1_cell = zeros(1,max(max(mask_image)));
    F_t2_cell = zeros(1,max(max(mask_image)));
    F_a1_cell = zeros(1,max(max(mask_image)));
    F_chi_cell = zeros(1,max(max(mask_image)));
    N_chi_cell = zeros(1,max(max(mask_image)));
    Cellpix_cell = zeros(1,max(max(mask_image)));
Cytopix_cell = zeros(1,max(max(mask_image)));
Cellnum_cell = zeros(1,max(max(mask_image)));
Imnum_cell = zeros(1,max(max(mask_image)));
     D_cell = zeros(1,max(max(mask_image)));  
     
    for i = 1:max(max(mask_image));
        
        [x y] = find(mask_image == i);
        if length(x)>0;
        rr_pix = zeros(1,length(x));
        N_tm_pix = zeros(1,length(x));
        N_t1_pix = zeros(1,length(x));
        N_t2_pix = zeros(1,length(x));
        N_a1_pix = zeros(1,length(x));
        F_tm_pix = zeros(1,length(x));
        F_t1_pix = zeros(1,length(x));
        F_t2_pix = zeros(1,length(x));
        F_a1_pix = zeros(1,length(x));
        F_chi_pix = zeros(1,length(x));
        N_chi_pix = zeros(1,length(x));
        
        for j = 1:length(x);
            if F_tm_image(x(j),y(j)) > 500;
            else
            rr_pix(j) = rr_image(x(j),y(j));
            N_tm_pix(j) = N_tm_image(x(j),y(j));
            N_t1_pix(j) = N_t1(x(j),y(j));
            N_t2_pix(j) = N_t2(x(j),y(j));
            N_a1_pix(j) = N_a1(x(j),y(j));
            F_tm_pix(j) = F_tm_image(x(j),y(j));
            F_t1_pix(j) = F_t1(x(j),y(j));
            F_t2_pix(j) = F_t2(x(j),y(j));
            F_a1_pix(j) = F_a1(x(j),y(j));
            N_chi_pix(j) = N_chi(x(j),y(j));
            F_chi_pix(j) = F_chi(x(j),y(j));
            end
        end
        if mean(N_tm_pix(N_tm_pix>00))<200; %filter to remove red blood cells
        elseif mean(F_tm_pix(F_tm_pix>00))<10; % filter to remove cells with fewer than 10 FAD pixels
        elseif mean(N_tm_pix)==0;
        elseif mean(F_tm_pix)==0;
        else
        rr_cell(i)= mean(rr_pix(rr_pix>0));
        N_tm_cell(i) = mean(N_tm_pix(N_tm_pix>0));
        N_t1_cell(i) = mean(N_t1_pix(N_t1_pix>0));
        N_t2_cell(i) = mean(N_t2_pix(N_t2_pix>0));
        N_a1_cell(i) = mean(N_a1_pix(N_a1_pix>0));
        F_tm_cell(i) = mean(F_tm_pix(F_tm_pix>0));
        F_t1_cell(i) = mean(F_t1_pix(F_t1_pix>0));
        F_t2_cell(i) = mean(F_t2_pix(F_t2_pix>0));
        F_a1_cell(i) = mean(F_a1_pix(F_a1_pix>0));
        
        [cell_x cell_y] = find(cell_image==i);
        Cellpix_cell(i) = length(cell_x);
        Cytopix_cell(i) = length(x);
        Cellnum_cell(i) = i;
        Imnum_cell(i)=a;
        D_cell(i) = pi*(mean([(max(x)-min(x)), (max(y)-min(y))]))^2;
       F_chi_cell(i) = mean(F_chi_pix(F_chi_pix<10));
        N_chi_cell(i) = mean(N_chi_pix(N_chi_pix<10));
        end
        end
    end
    ind1 = length(rr_out(rr_out>0))+1;
    ind2 = length(rr_out(rr_out>0))+length(rr_cell);
    rr_out(ind1:ind2) = rr_cell;
    N_tm_out(ind1:ind2) = N_tm_cell;
    N_t1_out(ind1:ind2) = N_t1_cell;
    N_t2_out(ind1:ind2) = N_t2_cell;
    N_a1_out(ind1:ind2) = N_a1_cell;
    F_tm_out(ind1:ind2) = F_tm_cell;
    F_t1_out(ind1:ind2) = F_t1_cell;
    F_t2_out(ind1:ind2) = F_t2_cell;
    F_a1_out(ind1:ind2) = F_a1_cell;
    F_chi_out(ind1:ind2) = F_chi_cell;
    N_chi_out(ind1:ind2) = N_chi_cell;
   Cellpix_out(ind1:ind2) = Cellpix_cell;
   Cytopix_out(ind1:ind2) = Cytopix_cell;
   Cellnum_out(ind1:ind2) = Cellnum_cell;
   Imnum_out(ind1:ind2)=Imnum_cell;
       D_out(ind1:ind2) = D_cell;
end
rr_out = rr_out(rr_out>0);
N_tm_out = N_tm_out(N_tm_out>0);
N_t1_out = N_t1_out(N_t1_out>0);
N_t2_out = N_t2_out(N_t2_out>0);
N_a1_out = N_a1_out(N_a1_out>0);
F_tm_out = F_tm_out(F_tm_out>0);
F_t1_out = F_t1_out(F_t1_out>0);
F_t2_out = F_t2_out(F_t2_out>0);
F_a1_out = F_a1_out(F_a1_out>0);
N_chi_out = N_chi_out(N_chi_out>0);
F_chi_out = F_chi_out(F_chi_out>0);
Cellpix_out = Cellpix_out(Cellpix_out>0);
Cytopix_out = Cytopix_out(Cytopix_out>0);
Cellnum_out = Cellnum_out(Cellnum_out>0);
Imnum_out=Imnum_out(Imnum_out>0);
D_out = D_out(D_out>0);

save('E:\020718\Dish 1\FLIM_out.mat','rr_out','N_tm_out','N_t1_out','N_t2_out','N_a1_out','F_tm_out','F_t1_out','F_t2_out','F_a1_out','N_chi_out','F_chi_out','Cellpix_out','Cytopix_out','Cellnum_out','Imnum_out','D_out');

