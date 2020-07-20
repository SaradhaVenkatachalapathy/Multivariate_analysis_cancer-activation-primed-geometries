%% 20190110_NIH3t3_p27_CM_1mg_col_Actingreen_DAPI_asma_647_mkl_568 control
cd 'E:\nih3t3_3d_cm\20190110_NIH3t3_p27_CM_1mg_col_Actingreen_DAPI_asma_647_mkl_568\control\zproject\indivisual_nuclei_ch1'
clear all
sample='*';
filenames = vertcat(dir([sample,'.nd2']),dir([sample,'.tif']));

'loading images'

XYZ={};
MetaData={};

parfor reader_count=1:size(filenames,1);
    filename = filenames(reader_count).name;
    [XYZ{1,reader_count},MetaData{1,reader_count}]=imread(filename,1);
end
'Files imported'


cd 'E:\nih3t3_3d_cm\matlab_pipeline\'

MetaData=[];
MetaData.Voxel_Size_X = 0.21; % in µm
MetaData.Voxel_Size_Y =MetaData.Voxel_Size_X;
MetaData.Filename=filenames.name;

Selected_Segmented_2d_Raw=XYZ;
Basic_Measurements=cell(size(Selected_Segmented_2d_Raw));
Boundary_RadOfCurva_K=cell(size(XYZ));
Norm_Pix_Intensities=cell(size(Selected_Segmented_2d_Raw));
num_of_pix_in_largest_nuc=zeros(size(Selected_Segmented_2d_Raw));
num_of_nuc=zeros(size(Selected_Segmented_2d_Raw));
total_num_of_nuc=zeros(size(Selected_Segmented_2d_Raw,1),1);

'Basic Measurements'

[Basic_Measurements,Boundary_RadOfCurva_K,Norm_Pix_Intensities]=basic_measurements(Selected_Segmented_2d_Raw,MetaData);

header={'Filename' 'Nuc label' 'Centroid X(um)' 'Centroid Y(um)' 'Pro. Area(um^2)' 'Perimeter(um)' 'A.R.' 'Shape Factor' 'PDI' 'Centre Mismatch' 'I80_by_I20' 'nHigh_by_nLow' 'Mean of Normalized Int' 'Median of Normalized Int' 'S.D of Normalized Int' 'Mode of Normalized Int' 'Entropy' 'Relative concavity'...
    'Num. of times curvature changes polarity' 'Avg Curvature(um^-1)' 'SD in curvature(um^-1)'...
    'Max positive curvature(um^-1)' 'Total positive curvature(um^-1)' 'Avg positive curvature(um^-1)' 'Total Positive Curvature>0.2(um^-1)' 'Avg Positive Curvature>0.2(um^-1)' 'Length with positive curvature>0.2(um)' 'Fraction of perimeter with positive curvature>0.2'...
    'Num positive curvatures > 0.25' 'Avg prominance of positive curvatures > 0.25(um)' 'Avg width of positive curvatures > 0.25(um)'...
    'Max negative curvature(um^-1)' 'Total negative curvature(um^-1)' 'Avg negative curvature(um^-1)' 'Total negative Curvature<-0.005(um^-1)' 'Avg negative Curvature<-0.005(um^-1)' 'Length with negative Curvature<-0.005(um^-1)' 'Fraction of perimeter with negative Curvature<-0.005(um^-1)'...
    'Num negative curvatures<-0.05' 'Avg prominance of negative curvatures<-0.05(um)' 'Avg width of negative curvatures<-0.05(um)'...
    'Average curvature x Perimeter'  'Avg positive curvature x Perimeter' 'Avg negative curvature x Perimeter'};


Result_Basic_Measurements=vertcat(header,Basic_Measurements);

clearvars condition_count file_count file_count1 upper_bound lower_bound num_of_nuc nuc_count Norm_Pix_Intensities Basic_Measurements header Boundary_RadOfCurva_K num_of_pix_in_largest_nuc number


% Lengthscales, Moment Invariants and Intensity Analyses
LengthScale_header={'Modal_LengthScale_FullNuc' 'Min_LengthScale_FullNuc' 'Avg_LengthScale_FullNuc' 'LengthScale_OfAvgCorr_FullNuc' ...
    'Modal_LengthScale_InCircle' 'Min_LengthScale_InCircle' 'Avg_LengthScale_InCircle' 'LengthScale_OfAvgCorr_InCircle' ...
    'SmallestNuclearSide' 'LargestNuclearSide' 'LargestToSmallestSide' 'avg_boundary_dist' 'SD_in_boundary_dist'};

MomentInvariants_header={'MomInv1' 'MomInv2' 'MomInv3' 'MomInv4' 'MomInv5' 'MomInv6' 'MomInv7' 'MomInv8'};

'Lengthscales and Moment Invariants'

CorrLengthScales=[];
MomentInv=[];
TotalCh1=[];
 for img_count=1:size(Selected_Segmented_2d_Raw,2)
        if ~isempty(Selected_Segmented_2d_Raw{img_count})
            CorrLengthScales_temp=[];
            MomentInv_temp=[];
            TotalCh1_temp=[];
            CorrLengthScales_temp(1,:)=Correlation_Lengthscales(Selected_Segmented_2d_Raw{img_count},MetaData.Voxel_Size_X);
            MomentInv_temp(1,:)=MomentInvariants(uint16(Selected_Segmented_2d_Raw{img_count}),MetaData.Voxel_Size_X);
            CorrLengthScales=vertcat(CorrLengthScales,CorrLengthScales_temp);
            MomentInv=vertcat(MomentInv,MomentInv_temp);
            
         end
    end


Result_Basic_Measurements=horzcat(Result_Basic_Measurements,...
        vertcat(LengthScale_header,num2cell(CorrLengthScales)),...
        vertcat(MomentInvariants_header,num2cell(MomentInv)));

cd 'E:\nih3t3_3d_cm\20190110_NIH3t3_p27_CM_1mg_col_Actingreen_DAPI_asma_647_mkl_568\control\2D measures\'   
T = cell2table(Result_Basic_Measurements(1:end,:));
writetable(T,'kamal_pipeline.csv')

close all

%% 20190110_NIH3t3_p27_CM_1mg_col_Actingreen_DAPI_asma_647_mkl_568 T1
cd 'E:\nih3t3_3d_cm\20190110_NIH3t3_p27_CM_1mg_col_Actingreen_DAPI_asma_647_mkl_568\T1\zproject\indivisual_nuclei_ch1'
clear all
sample='*';
filenames = vertcat(dir([sample,'.nd2']),dir([sample,'.tif']));

'loading images'

XYZ={};
MetaData={};

parfor reader_count=1:size(filenames,1);
    filename = filenames(reader_count).name;
    [XYZ{1,reader_count},MetaData{1,reader_count}]=imread(filename,1);
end
'Files imported'


cd 'E:\nih3t3_3d_cm\matlab_pipeline\'

MetaData=[];
MetaData.Voxel_Size_X = 0.21; % in µm
MetaData.Voxel_Size_Y =MetaData.Voxel_Size_X;
MetaData.Filename=filenames.name;

Selected_Segmented_2d_Raw=XYZ;
Basic_Measurements=cell(size(Selected_Segmented_2d_Raw));
Boundary_RadOfCurva_K=cell(size(XYZ));
Norm_Pix_Intensities=cell(size(Selected_Segmented_2d_Raw));
num_of_pix_in_largest_nuc=zeros(size(Selected_Segmented_2d_Raw));
num_of_nuc=zeros(size(Selected_Segmented_2d_Raw));
total_num_of_nuc=zeros(size(Selected_Segmented_2d_Raw,1),1);

'Basic Measurements'

[Basic_Measurements,Boundary_RadOfCurva_K,Norm_Pix_Intensities]=basic_measurements(Selected_Segmented_2d_Raw,MetaData);

header={'Filename' 'Nuc label' 'Centroid X(um)' 'Centroid Y(um)' 'Pro. Area(um^2)' 'Perimeter(um)' 'A.R.' 'Shape Factor' 'PDI' 'Centre Mismatch' 'I80_by_I20' 'nHigh_by_nLow' 'Mean of Normalized Int' 'Median of Normalized Int' 'S.D of Normalized Int' 'Mode of Normalized Int' 'Entropy' 'Relative concavity'...
    'Num. of times curvature changes polarity' 'Avg Curvature(um^-1)' 'SD in curvature(um^-1)'...
    'Max positive curvature(um^-1)' 'Total positive curvature(um^-1)' 'Avg positive curvature(um^-1)' 'Total Positive Curvature>0.2(um^-1)' 'Avg Positive Curvature>0.2(um^-1)' 'Length with positive curvature>0.2(um)' 'Fraction of perimeter with positive curvature>0.2'...
    'Num positive curvatures > 0.25' 'Avg prominance of positive curvatures > 0.25(um)' 'Avg width of positive curvatures > 0.25(um)'...
    'Max negative curvature(um^-1)' 'Total negative curvature(um^-1)' 'Avg negative curvature(um^-1)' 'Total negative Curvature<-0.005(um^-1)' 'Avg negative Curvature<-0.005(um^-1)' 'Length with negative Curvature<-0.005(um^-1)' 'Fraction of perimeter with negative Curvature<-0.005(um^-1)'...
    'Num negative curvatures<-0.05' 'Avg prominance of negative curvatures<-0.05(um)' 'Avg width of negative curvatures<-0.05(um)'...
    'Average curvature x Perimeter'  'Avg positive curvature x Perimeter' 'Avg negative curvature x Perimeter'};


Result_Basic_Measurements=vertcat(header,Basic_Measurements);

clearvars condition_count file_count file_count1 upper_bound lower_bound num_of_nuc nuc_count Norm_Pix_Intensities Basic_Measurements header Boundary_RadOfCurva_K num_of_pix_in_largest_nuc number


% Lengthscales, Moment Invariants and Intensity Analyses
LengthScale_header={'Modal_LengthScale_FullNuc' 'Min_LengthScale_FullNuc' 'Avg_LengthScale_FullNuc' 'LengthScale_OfAvgCorr_FullNuc' ...
    'Modal_LengthScale_InCircle' 'Min_LengthScale_InCircle' 'Avg_LengthScale_InCircle' 'LengthScale_OfAvgCorr_InCircle' ...
    'SmallestNuclearSide' 'LargestNuclearSide' 'LargestToSmallestSide' 'avg_boundary_dist' 'SD_in_boundary_dist'};

MomentInvariants_header={'MomInv1' 'MomInv2' 'MomInv3' 'MomInv4' 'MomInv5' 'MomInv6' 'MomInv7' 'MomInv8'};

'Lengthscales and Moment Invariants'

CorrLengthScales=[];
MomentInv=[];
TotalCh1=[];
 for img_count=1:size(Selected_Segmented_2d_Raw,2)
        if ~isempty(Selected_Segmented_2d_Raw{img_count})
            CorrLengthScales_temp=[];
            MomentInv_temp=[];
            TotalCh1_temp=[];
            CorrLengthScales_temp(1,:)=Correlation_Lengthscales(Selected_Segmented_2d_Raw{img_count},MetaData.Voxel_Size_X);
            MomentInv_temp(1,:)=MomentInvariants(uint16(Selected_Segmented_2d_Raw{img_count}),MetaData.Voxel_Size_X);
            CorrLengthScales=vertcat(CorrLengthScales,CorrLengthScales_temp);
            MomentInv=vertcat(MomentInv,MomentInv_temp);
            
         end
    end


Result_Basic_Measurements=horzcat(Result_Basic_Measurements,...
        vertcat(LengthScale_header,num2cell(CorrLengthScales)),...
        vertcat(MomentInvariants_header,num2cell(MomentInv)));

cd 'E:\nih3t3_3d_cm\20190110_NIH3t3_p27_CM_1mg_col_Actingreen_DAPI_asma_647_mkl_568\T1\2D measures\'   
T = cell2table(Result_Basic_Measurements(1:end,:));
writetable(T,'kamal_pipeline.csv')
close all


%% 20190110_NIH3t3_p27_CM_1mg_col_Actingreen_DAPI_asma_647_mkl_568 T2
cd 'E:\nih3t3_3d_cm\20190110_NIH3t3_p27_CM_1mg_col_Actingreen_DAPI_asma_647_mkl_568\T2\zproject\indivisual_nuclei_ch1'
clear all
sample='*';
filenames = vertcat(dir([sample,'.nd2']),dir([sample,'.tif']));

'loading images'

XYZ={};
MetaData={};

parfor reader_count=1:size(filenames,1);
    filename = filenames(reader_count).name;
    [XYZ{1,reader_count},MetaData{1,reader_count}]=imread(filename,1);
end
'Files imported'


cd 'E:\nih3t3_3d_cm\matlab_pipeline\'

MetaData=[];
MetaData.Voxel_Size_X = 0.21; % in µm
MetaData.Voxel_Size_Y =MetaData.Voxel_Size_X;
MetaData.Filename=filenames.name;

Selected_Segmented_2d_Raw=XYZ;
Basic_Measurements=cell(size(Selected_Segmented_2d_Raw));
Boundary_RadOfCurva_K=cell(size(XYZ));
Norm_Pix_Intensities=cell(size(Selected_Segmented_2d_Raw));
num_of_pix_in_largest_nuc=zeros(size(Selected_Segmented_2d_Raw));
num_of_nuc=zeros(size(Selected_Segmented_2d_Raw));
total_num_of_nuc=zeros(size(Selected_Segmented_2d_Raw,1),1);

'Basic Measurements'

[Basic_Measurements,Boundary_RadOfCurva_K,Norm_Pix_Intensities]=basic_measurements(Selected_Segmented_2d_Raw,MetaData);

header={'Filename' 'Nuc label' 'Centroid X(um)' 'Centroid Y(um)' 'Pro. Area(um^2)' 'Perimeter(um)' 'A.R.' 'Shape Factor' 'PDI' 'Centre Mismatch' 'I80_by_I20' 'nHigh_by_nLow' 'Mean of Normalized Int' 'Median of Normalized Int' 'S.D of Normalized Int' 'Mode of Normalized Int' 'Entropy' 'Relative concavity'...
    'Num. of times curvature changes polarity' 'Avg Curvature(um^-1)' 'SD in curvature(um^-1)'...
    'Max positive curvature(um^-1)' 'Total positive curvature(um^-1)' 'Avg positive curvature(um^-1)' 'Total Positive Curvature>0.2(um^-1)' 'Avg Positive Curvature>0.2(um^-1)' 'Length with positive curvature>0.2(um)' 'Fraction of perimeter with positive curvature>0.2'...
    'Num positive curvatures > 0.25' 'Avg prominance of positive curvatures > 0.25(um)' 'Avg width of positive curvatures > 0.25(um)'...
    'Max negative curvature(um^-1)' 'Total negative curvature(um^-1)' 'Avg negative curvature(um^-1)' 'Total negative Curvature<-0.005(um^-1)' 'Avg negative Curvature<-0.005(um^-1)' 'Length with negative Curvature<-0.005(um^-1)' 'Fraction of perimeter with negative Curvature<-0.005(um^-1)'...
    'Num negative curvatures<-0.05' 'Avg prominance of negative curvatures<-0.05(um)' 'Avg width of negative curvatures<-0.05(um)'...
    'Average curvature x Perimeter'  'Avg positive curvature x Perimeter' 'Avg negative curvature x Perimeter'};


Result_Basic_Measurements=vertcat(header,Basic_Measurements);

clearvars condition_count file_count file_count1 upper_bound lower_bound num_of_nuc nuc_count Norm_Pix_Intensities Basic_Measurements header Boundary_RadOfCurva_K num_of_pix_in_largest_nuc number


% Lengthscales, Moment Invariants and Intensity Analyses
LengthScale_header={'Modal_LengthScale_FullNuc' 'Min_LengthScale_FullNuc' 'Avg_LengthScale_FullNuc' 'LengthScale_OfAvgCorr_FullNuc' ...
    'Modal_LengthScale_InCircle' 'Min_LengthScale_InCircle' 'Avg_LengthScale_InCircle' 'LengthScale_OfAvgCorr_InCircle' ...
    'SmallestNuclearSide' 'LargestNuclearSide' 'LargestToSmallestSide' 'avg_boundary_dist' 'SD_in_boundary_dist'};

MomentInvariants_header={'MomInv1' 'MomInv2' 'MomInv3' 'MomInv4' 'MomInv5' 'MomInv6' 'MomInv7' 'MomInv8'};

'Lengthscales and Moment Invariants'

CorrLengthScales=[];
MomentInv=[];
TotalCh1=[];
 for img_count=1:size(Selected_Segmented_2d_Raw,2)
        if ~isempty(Selected_Segmented_2d_Raw{img_count})
            CorrLengthScales_temp=[];
            MomentInv_temp=[];
            TotalCh1_temp=[];
            CorrLengthScales_temp(1,:)=Correlation_Lengthscales(Selected_Segmented_2d_Raw{img_count},MetaData.Voxel_Size_X);
            MomentInv_temp(1,:)=MomentInvariants(uint16(Selected_Segmented_2d_Raw{img_count}),MetaData.Voxel_Size_X);
            CorrLengthScales=vertcat(CorrLengthScales,CorrLengthScales_temp);
            MomentInv=vertcat(MomentInv,MomentInv_temp);
            
         end
    end


Result_Basic_Measurements=horzcat(Result_Basic_Measurements,...
        vertcat(LengthScale_header,num2cell(CorrLengthScales)),...
        vertcat(MomentInvariants_header,num2cell(MomentInv)));

cd 'E:\nih3t3_3d_cm\20190110_NIH3t3_p27_CM_1mg_col_Actingreen_DAPI_asma_647_mkl_568\T2\2D measures\'   
T = cell2table(Result_Basic_Measurements(1:end,:));
writetable(T,'kamal_pipeline.csv')

close all

%% 20190117_NIH3t3_3d_1mg_col_3days_MCF7_cm_asma_647_mkl_568_actin_488_DAPI control
cd 'E:\nih3t3_3d_cm\20190117_NIH3t3_3d_1mg_col_3days_MCF7_cm_asma_647_mkl_568_actin_488_DAPI\control\zproject\indivisual_nuclei_ch1'
clear all
sample='*';
filenames = vertcat(dir([sample,'.nd2']),dir([sample,'.tif']));

'loading images'

XYZ={};
MetaData={};

parfor reader_count=1:size(filenames,1);
    filename = filenames(reader_count).name;
    [XYZ{1,reader_count},MetaData{1,reader_count}]=imread(filename,1);
end
'Files imported'


cd 'E:\nih3t3_3d_cm\matlab_pipeline\'

MetaData=[];
MetaData.Voxel_Size_X = 0.21; % in µm
MetaData.Voxel_Size_Y =MetaData.Voxel_Size_X;
MetaData.Filename=filenames.name;

Selected_Segmented_2d_Raw=XYZ;
Basic_Measurements=cell(size(Selected_Segmented_2d_Raw));
Boundary_RadOfCurva_K=cell(size(XYZ));
Norm_Pix_Intensities=cell(size(Selected_Segmented_2d_Raw));
num_of_pix_in_largest_nuc=zeros(size(Selected_Segmented_2d_Raw));
num_of_nuc=zeros(size(Selected_Segmented_2d_Raw));
total_num_of_nuc=zeros(size(Selected_Segmented_2d_Raw,1),1);

'Basic Measurements'

[Basic_Measurements,Boundary_RadOfCurva_K,Norm_Pix_Intensities]=basic_measurements(Selected_Segmented_2d_Raw,MetaData);

header={'Filename' 'Nuc label' 'Centroid X(um)' 'Centroid Y(um)' 'Pro. Area(um^2)' 'Perimeter(um)' 'A.R.' 'Shape Factor' 'PDI' 'Centre Mismatch' 'I80_by_I20' 'nHigh_by_nLow' 'Mean of Normalized Int' 'Median of Normalized Int' 'S.D of Normalized Int' 'Mode of Normalized Int' 'Entropy' 'Relative concavity'...
    'Num. of times curvature changes polarity' 'Avg Curvature(um^-1)' 'SD in curvature(um^-1)'...
    'Max positive curvature(um^-1)' 'Total positive curvature(um^-1)' 'Avg positive curvature(um^-1)' 'Total Positive Curvature>0.2(um^-1)' 'Avg Positive Curvature>0.2(um^-1)' 'Length with positive curvature>0.2(um)' 'Fraction of perimeter with positive curvature>0.2'...
    'Num positive curvatures > 0.25' 'Avg prominance of positive curvatures > 0.25(um)' 'Avg width of positive curvatures > 0.25(um)'...
    'Max negative curvature(um^-1)' 'Total negative curvature(um^-1)' 'Avg negative curvature(um^-1)' 'Total negative Curvature<-0.005(um^-1)' 'Avg negative Curvature<-0.005(um^-1)' 'Length with negative Curvature<-0.005(um^-1)' 'Fraction of perimeter with negative Curvature<-0.005(um^-1)'...
    'Num negative curvatures<-0.05' 'Avg prominance of negative curvatures<-0.05(um)' 'Avg width of negative curvatures<-0.05(um)'...
    'Average curvature x Perimeter'  'Avg positive curvature x Perimeter' 'Avg negative curvature x Perimeter'};


Result_Basic_Measurements=vertcat(header,Basic_Measurements);

clearvars condition_count file_count file_count1 upper_bound lower_bound num_of_nuc nuc_count Norm_Pix_Intensities Basic_Measurements header Boundary_RadOfCurva_K num_of_pix_in_largest_nuc number


% Lengthscales, Moment Invariants and Intensity Analyses
LengthScale_header={'Modal_LengthScale_FullNuc' 'Min_LengthScale_FullNuc' 'Avg_LengthScale_FullNuc' 'LengthScale_OfAvgCorr_FullNuc' ...
    'Modal_LengthScale_InCircle' 'Min_LengthScale_InCircle' 'Avg_LengthScale_InCircle' 'LengthScale_OfAvgCorr_InCircle' ...
    'SmallestNuclearSide' 'LargestNuclearSide' 'LargestToSmallestSide' 'avg_boundary_dist' 'SD_in_boundary_dist'};

MomentInvariants_header={'MomInv1' 'MomInv2' 'MomInv3' 'MomInv4' 'MomInv5' 'MomInv6' 'MomInv7' 'MomInv8'};

'Lengthscales and Moment Invariants'

CorrLengthScales=[];
MomentInv=[];
TotalCh1=[];
 for img_count=1:size(Selected_Segmented_2d_Raw,2)
        if ~isempty(Selected_Segmented_2d_Raw{img_count})
            CorrLengthScales_temp=[];
            MomentInv_temp=[];
            TotalCh1_temp=[];
            CorrLengthScales_temp(1,:)=Correlation_Lengthscales(Selected_Segmented_2d_Raw{img_count},MetaData.Voxel_Size_X);
            MomentInv_temp(1,:)=MomentInvariants(uint16(Selected_Segmented_2d_Raw{img_count}),MetaData.Voxel_Size_X);
            CorrLengthScales=vertcat(CorrLengthScales,CorrLengthScales_temp);
            MomentInv=vertcat(MomentInv,MomentInv_temp);
            
         end
    end


Result_Basic_Measurements=horzcat(Result_Basic_Measurements,...
        vertcat(LengthScale_header,num2cell(CorrLengthScales)),...
        vertcat(MomentInvariants_header,num2cell(MomentInv)));

cd 'E:\nih3t3_3d_cm\20190117_NIH3t3_3d_1mg_col_3days_MCF7_cm_asma_647_mkl_568_actin_488_DAPI\control\2D measures\'   
T = cell2table(Result_Basic_Measurements(1:end,:));
writetable(T,'kamal_pipeline.csv')

close all

%% 20190117_NIH3t3_3d_1mg_col_3days_MCF7_cm_asma_647_mkl_568_actin_488_DAPI T1
cd 'E:\nih3t3_3d_cm\20190117_NIH3t3_3d_1mg_col_3days_MCF7_cm_asma_647_mkl_568_actin_488_DAPI\T1\zproject\indivisual_nuclei_ch1'
clear all
sample='*';
filenames = vertcat(dir([sample,'.nd2']),dir([sample,'.tif']));

'loading images'

XYZ={};
MetaData={};

parfor reader_count=1:size(filenames,1);
    filename = filenames(reader_count).name;
    [XYZ{1,reader_count},MetaData{1,reader_count}]=imread(filename,1);
end
'Files imported'


cd 'E:\nih3t3_3d_cm\matlab_pipeline\'

MetaData=[];
MetaData.Voxel_Size_X = 0.21; % in µm
MetaData.Voxel_Size_Y =MetaData.Voxel_Size_X;
MetaData.Filename=filenames.name;

Selected_Segmented_2d_Raw=XYZ;
Basic_Measurements=cell(size(Selected_Segmented_2d_Raw));
Boundary_RadOfCurva_K=cell(size(XYZ));
Norm_Pix_Intensities=cell(size(Selected_Segmented_2d_Raw));
num_of_pix_in_largest_nuc=zeros(size(Selected_Segmented_2d_Raw));
num_of_nuc=zeros(size(Selected_Segmented_2d_Raw));
total_num_of_nuc=zeros(size(Selected_Segmented_2d_Raw,1),1);

'Basic Measurements'

[Basic_Measurements,Boundary_RadOfCurva_K,Norm_Pix_Intensities]=basic_measurements(Selected_Segmented_2d_Raw,MetaData);

header={'Filename' 'Nuc label' 'Centroid X(um)' 'Centroid Y(um)' 'Pro. Area(um^2)' 'Perimeter(um)' 'A.R.' 'Shape Factor' 'PDI' 'Centre Mismatch' 'I80_by_I20' 'nHigh_by_nLow' 'Mean of Normalized Int' 'Median of Normalized Int' 'S.D of Normalized Int' 'Mode of Normalized Int' 'Entropy' 'Relative concavity'...
    'Num. of times curvature changes polarity' 'Avg Curvature(um^-1)' 'SD in curvature(um^-1)'...
    'Max positive curvature(um^-1)' 'Total positive curvature(um^-1)' 'Avg positive curvature(um^-1)' 'Total Positive Curvature>0.2(um^-1)' 'Avg Positive Curvature>0.2(um^-1)' 'Length with positive curvature>0.2(um)' 'Fraction of perimeter with positive curvature>0.2'...
    'Num positive curvatures > 0.25' 'Avg prominance of positive curvatures > 0.25(um)' 'Avg width of positive curvatures > 0.25(um)'...
    'Max negative curvature(um^-1)' 'Total negative curvature(um^-1)' 'Avg negative curvature(um^-1)' 'Total negative Curvature<-0.005(um^-1)' 'Avg negative Curvature<-0.005(um^-1)' 'Length with negative Curvature<-0.005(um^-1)' 'Fraction of perimeter with negative Curvature<-0.005(um^-1)'...
    'Num negative curvatures<-0.05' 'Avg prominance of negative curvatures<-0.05(um)' 'Avg width of negative curvatures<-0.05(um)'...
    'Average curvature x Perimeter'  'Avg positive curvature x Perimeter' 'Avg negative curvature x Perimeter'};


Result_Basic_Measurements=vertcat(header,Basic_Measurements);

clearvars condition_count file_count file_count1 upper_bound lower_bound num_of_nuc nuc_count Norm_Pix_Intensities Basic_Measurements header Boundary_RadOfCurva_K num_of_pix_in_largest_nuc number


% Lengthscales, Moment Invariants and Intensity Analyses
LengthScale_header={'Modal_LengthScale_FullNuc' 'Min_LengthScale_FullNuc' 'Avg_LengthScale_FullNuc' 'LengthScale_OfAvgCorr_FullNuc' ...
    'Modal_LengthScale_InCircle' 'Min_LengthScale_InCircle' 'Avg_LengthScale_InCircle' 'LengthScale_OfAvgCorr_InCircle' ...
    'SmallestNuclearSide' 'LargestNuclearSide' 'LargestToSmallestSide' 'avg_boundary_dist' 'SD_in_boundary_dist'};

MomentInvariants_header={'MomInv1' 'MomInv2' 'MomInv3' 'MomInv4' 'MomInv5' 'MomInv6' 'MomInv7' 'MomInv8'};

'Lengthscales and Moment Invariants'

CorrLengthScales=[];
MomentInv=[];
TotalCh1=[];
 for img_count=1:size(Selected_Segmented_2d_Raw,2)
        if ~isempty(Selected_Segmented_2d_Raw{img_count})
            CorrLengthScales_temp=[];
            MomentInv_temp=[];
            TotalCh1_temp=[];
            CorrLengthScales_temp(1,:)=Correlation_Lengthscales(Selected_Segmented_2d_Raw{img_count},MetaData.Voxel_Size_X);
            MomentInv_temp(1,:)=MomentInvariants(uint16(Selected_Segmented_2d_Raw{img_count}),MetaData.Voxel_Size_X);
            CorrLengthScales=vertcat(CorrLengthScales,CorrLengthScales_temp);
            MomentInv=vertcat(MomentInv,MomentInv_temp);
            
         end
    end


Result_Basic_Measurements=horzcat(Result_Basic_Measurements,...
        vertcat(LengthScale_header,num2cell(CorrLengthScales)),...
        vertcat(MomentInvariants_header,num2cell(MomentInv)));

cd 'E:\nih3t3_3d_cm\20190117_NIH3t3_3d_1mg_col_3days_MCF7_cm_asma_647_mkl_568_actin_488_DAPI\T1\2D measures\'   
T = cell2table(Result_Basic_Measurements(1:end,:));
writetable(T,'kamal_pipeline.csv')
close all


%% 20190117_NIH3t3_3d_1mg_col_3days_MCF7_cm_asma_647_mkl_568_actin_488_DAPI T2
cd 'E:\nih3t3_3d_cm\20190117_NIH3t3_3d_1mg_col_3days_MCF7_cm_asma_647_mkl_568_actin_488_DAPI\T2\zproject\indivisual_nuclei_ch1'
clear all
sample='*';
filenames = vertcat(dir([sample,'.nd2']),dir([sample,'.tif']));

'loading images'

XYZ={};
MetaData={};

parfor reader_count=1:size(filenames,1);
    filename = filenames(reader_count).name;
    [XYZ{1,reader_count},MetaData{1,reader_count}]=imread(filename,1);
end
'Files imported'


cd 'E:\nih3t3_3d_cm\matlab_pipeline\'

MetaData=[];
MetaData.Voxel_Size_X = 0.21; % in µm
MetaData.Voxel_Size_Y =MetaData.Voxel_Size_X;
MetaData.Filename=filenames.name;

Selected_Segmented_2d_Raw=XYZ;
Basic_Measurements=cell(size(Selected_Segmented_2d_Raw));
Boundary_RadOfCurva_K=cell(size(XYZ));
Norm_Pix_Intensities=cell(size(Selected_Segmented_2d_Raw));
num_of_pix_in_largest_nuc=zeros(size(Selected_Segmented_2d_Raw));
num_of_nuc=zeros(size(Selected_Segmented_2d_Raw));
total_num_of_nuc=zeros(size(Selected_Segmented_2d_Raw,1),1);

'Basic Measurements'

[Basic_Measurements,Boundary_RadOfCurva_K,Norm_Pix_Intensities]=basic_measurements(Selected_Segmented_2d_Raw,MetaData);

header={'Filename' 'Nuc label' 'Centroid X(um)' 'Centroid Y(um)' 'Pro. Area(um^2)' 'Perimeter(um)' 'A.R.' 'Shape Factor' 'PDI' 'Centre Mismatch' 'I80_by_I20' 'nHigh_by_nLow' 'Mean of Normalized Int' 'Median of Normalized Int' 'S.D of Normalized Int' 'Mode of Normalized Int' 'Entropy' 'Relative concavity'...
    'Num. of times curvature changes polarity' 'Avg Curvature(um^-1)' 'SD in curvature(um^-1)'...
    'Max positive curvature(um^-1)' 'Total positive curvature(um^-1)' 'Avg positive curvature(um^-1)' 'Total Positive Curvature>0.2(um^-1)' 'Avg Positive Curvature>0.2(um^-1)' 'Length with positive curvature>0.2(um)' 'Fraction of perimeter with positive curvature>0.2'...
    'Num positive curvatures > 0.25' 'Avg prominance of positive curvatures > 0.25(um)' 'Avg width of positive curvatures > 0.25(um)'...
    'Max negative curvature(um^-1)' 'Total negative curvature(um^-1)' 'Avg negative curvature(um^-1)' 'Total negative Curvature<-0.005(um^-1)' 'Avg negative Curvature<-0.005(um^-1)' 'Length with negative Curvature<-0.005(um^-1)' 'Fraction of perimeter with negative Curvature<-0.005(um^-1)'...
    'Num negative curvatures<-0.05' 'Avg prominance of negative curvatures<-0.05(um)' 'Avg width of negative curvatures<-0.05(um)'...
    'Average curvature x Perimeter'  'Avg positive curvature x Perimeter' 'Avg negative curvature x Perimeter'};


Result_Basic_Measurements=vertcat(header,Basic_Measurements);

clearvars condition_count file_count file_count1 upper_bound lower_bound num_of_nuc nuc_count Norm_Pix_Intensities Basic_Measurements header Boundary_RadOfCurva_K num_of_pix_in_largest_nuc number


% Lengthscales, Moment Invariants and Intensity Analyses
LengthScale_header={'Modal_LengthScale_FullNuc' 'Min_LengthScale_FullNuc' 'Avg_LengthScale_FullNuc' 'LengthScale_OfAvgCorr_FullNuc' ...
    'Modal_LengthScale_InCircle' 'Min_LengthScale_InCircle' 'Avg_LengthScale_InCircle' 'LengthScale_OfAvgCorr_InCircle' ...
    'SmallestNuclearSide' 'LargestNuclearSide' 'LargestToSmallestSide' 'avg_boundary_dist' 'SD_in_boundary_dist'};

MomentInvariants_header={'MomInv1' 'MomInv2' 'MomInv3' 'MomInv4' 'MomInv5' 'MomInv6' 'MomInv7' 'MomInv8'};

'Lengthscales and Moment Invariants'

CorrLengthScales=[];
MomentInv=[];
TotalCh1=[];
 for img_count=1:size(Selected_Segmented_2d_Raw,2)
        if ~isempty(Selected_Segmented_2d_Raw{img_count})
            CorrLengthScales_temp=[];
            MomentInv_temp=[];
            TotalCh1_temp=[];
            CorrLengthScales_temp(1,:)=Correlation_Lengthscales(Selected_Segmented_2d_Raw{img_count},MetaData.Voxel_Size_X);
            MomentInv_temp(1,:)=MomentInvariants(uint16(Selected_Segmented_2d_Raw{img_count}),MetaData.Voxel_Size_X);
            CorrLengthScales=vertcat(CorrLengthScales,CorrLengthScales_temp);
            MomentInv=vertcat(MomentInv,MomentInv_temp);
            
         end
    end


Result_Basic_Measurements=horzcat(Result_Basic_Measurements,...
        vertcat(LengthScale_header,num2cell(CorrLengthScales)),...
        vertcat(MomentInvariants_header,num2cell(MomentInv)));

cd 'E:\nih3t3_3d_cm\20190117_NIH3t3_3d_1mg_col_3days_MCF7_cm_asma_647_mkl_568_actin_488_DAPI\T2\2D measures\'   
T = cell2table(Result_Basic_Measurements(1:end,:));
writetable(T,'kamal_pipeline.csv')
close all

%% 20190227_NIH3T3_3d_reprgrammed_spheroid_cm_aSMA_647_DAPI control
cd 'E:\nih3t3_3d_cm\20190227_NIH3T3_3d_reprgrammed_spheroid_cm_aSMA_647_DAPI\control\zproject\indivisual_nuclei_ch1'
clear all
sample='*';
filenames = vertcat(dir([sample,'.nd2']),dir([sample,'.tif']));

'loading images'

XYZ={};
MetaData={};

parfor reader_count=1:size(filenames,1);
    filename = filenames(reader_count).name;
    [XYZ{1,reader_count},MetaData{1,reader_count}]=imread(filename,1);
end
'Files imported'


cd 'E:\nih3t3_3d_cm\matlab_pipeline\'

MetaData=[];
MetaData.Voxel_Size_X = 0.21; % in µm
MetaData.Voxel_Size_Y =MetaData.Voxel_Size_X;
MetaData.Filename=filenames.name;

Selected_Segmented_2d_Raw=XYZ;
Basic_Measurements=cell(size(Selected_Segmented_2d_Raw));
Boundary_RadOfCurva_K=cell(size(XYZ));
Norm_Pix_Intensities=cell(size(Selected_Segmented_2d_Raw));
num_of_pix_in_largest_nuc=zeros(size(Selected_Segmented_2d_Raw));
num_of_nuc=zeros(size(Selected_Segmented_2d_Raw));
total_num_of_nuc=zeros(size(Selected_Segmented_2d_Raw,1),1);

'Basic Measurements'

[Basic_Measurements,Boundary_RadOfCurva_K,Norm_Pix_Intensities]=basic_measurements(Selected_Segmented_2d_Raw,MetaData);

header={'Filename' 'Nuc label' 'Centroid X(um)' 'Centroid Y(um)' 'Pro. Area(um^2)' 'Perimeter(um)' 'A.R.' 'Shape Factor' 'PDI' 'Centre Mismatch' 'I80_by_I20' 'nHigh_by_nLow' 'Mean of Normalized Int' 'Median of Normalized Int' 'S.D of Normalized Int' 'Mode of Normalized Int' 'Entropy' 'Relative concavity'...
    'Num. of times curvature changes polarity' 'Avg Curvature(um^-1)' 'SD in curvature(um^-1)'...
    'Max positive curvature(um^-1)' 'Total positive curvature(um^-1)' 'Avg positive curvature(um^-1)' 'Total Positive Curvature>0.2(um^-1)' 'Avg Positive Curvature>0.2(um^-1)' 'Length with positive curvature>0.2(um)' 'Fraction of perimeter with positive curvature>0.2'...
    'Num positive curvatures > 0.25' 'Avg prominance of positive curvatures > 0.25(um)' 'Avg width of positive curvatures > 0.25(um)'...
    'Max negative curvature(um^-1)' 'Total negative curvature(um^-1)' 'Avg negative curvature(um^-1)' 'Total negative Curvature<-0.005(um^-1)' 'Avg negative Curvature<-0.005(um^-1)' 'Length with negative Curvature<-0.005(um^-1)' 'Fraction of perimeter with negative Curvature<-0.005(um^-1)'...
    'Num negative curvatures<-0.05' 'Avg prominance of negative curvatures<-0.05(um)' 'Avg width of negative curvatures<-0.05(um)'...
    'Average curvature x Perimeter'  'Avg positive curvature x Perimeter' 'Avg negative curvature x Perimeter'};


Result_Basic_Measurements=vertcat(header,Basic_Measurements);

clearvars condition_count file_count file_count1 upper_bound lower_bound num_of_nuc nuc_count Norm_Pix_Intensities Basic_Measurements header Boundary_RadOfCurva_K num_of_pix_in_largest_nuc number


% Lengthscales, Moment Invariants and Intensity Analyses
LengthScale_header={'Modal_LengthScale_FullNuc' 'Min_LengthScale_FullNuc' 'Avg_LengthScale_FullNuc' 'LengthScale_OfAvgCorr_FullNuc' ...
    'Modal_LengthScale_InCircle' 'Min_LengthScale_InCircle' 'Avg_LengthScale_InCircle' 'LengthScale_OfAvgCorr_InCircle' ...
    'SmallestNuclearSide' 'LargestNuclearSide' 'LargestToSmallestSide' 'avg_boundary_dist' 'SD_in_boundary_dist'};

MomentInvariants_header={'MomInv1' 'MomInv2' 'MomInv3' 'MomInv4' 'MomInv5' 'MomInv6' 'MomInv7' 'MomInv8'};

'Lengthscales and Moment Invariants'

CorrLengthScales=[];
MomentInv=[];
TotalCh1=[];
 for img_count=1:size(Selected_Segmented_2d_Raw,2)
        if ~isempty(Selected_Segmented_2d_Raw{img_count})
            CorrLengthScales_temp=[];
            MomentInv_temp=[];
            TotalCh1_temp=[];
            CorrLengthScales_temp(1,:)=Correlation_Lengthscales(Selected_Segmented_2d_Raw{img_count},MetaData.Voxel_Size_X);
            MomentInv_temp(1,:)=MomentInvariants(uint16(Selected_Segmented_2d_Raw{img_count}),MetaData.Voxel_Size_X);
            CorrLengthScales=vertcat(CorrLengthScales,CorrLengthScales_temp);
            MomentInv=vertcat(MomentInv,MomentInv_temp);
            
         end
    end


Result_Basic_Measurements=horzcat(Result_Basic_Measurements,...
        vertcat(LengthScale_header,num2cell(CorrLengthScales)),...
        vertcat(MomentInvariants_header,num2cell(MomentInv)));

cd 'E:\nih3t3_3d_cm\20190227_NIH3T3_3d_reprgrammed_spheroid_cm_aSMA_647_DAPI\control\2D measures\'   
T = cell2table(Result_Basic_Measurements(1:end,:));
writetable(T,'kamal_pipeline.csv')
close all


%% 20190227_NIH3T3_3d_reprgrammed_spheroid_cm_aSMA_647_DAPI T1
cd 'E:\nih3t3_3d_cm\20190227_NIH3T3_3d_reprgrammed_spheroid_cm_aSMA_647_DAPI\CM\zproject\indivisual_nuclei_ch1'
clear all
sample='*';
filenames = vertcat(dir([sample,'.nd2']),dir([sample,'.tif']));

'loading images'

XYZ={};
MetaData={};

parfor reader_count=1:size(filenames,1);
    filename = filenames(reader_count).name;
    [XYZ{1,reader_count},MetaData{1,reader_count}]=imread(filename,1);
end
'Files imported'


cd 'E:\nih3t3_3d_cm\matlab_pipeline\'

MetaData=[];
MetaData.Voxel_Size_X = 0.21; % in µm
MetaData.Voxel_Size_Y =MetaData.Voxel_Size_X;
MetaData.Filename=filenames.name;

Selected_Segmented_2d_Raw=XYZ;
Basic_Measurements=cell(size(Selected_Segmented_2d_Raw));
Boundary_RadOfCurva_K=cell(size(XYZ));
Norm_Pix_Intensities=cell(size(Selected_Segmented_2d_Raw));
num_of_pix_in_largest_nuc=zeros(size(Selected_Segmented_2d_Raw));
num_of_nuc=zeros(size(Selected_Segmented_2d_Raw));
total_num_of_nuc=zeros(size(Selected_Segmented_2d_Raw,1),1);

'Basic Measurements'

[Basic_Measurements,Boundary_RadOfCurva_K,Norm_Pix_Intensities]=basic_measurements(Selected_Segmented_2d_Raw,MetaData);

header={'Filename' 'Nuc label' 'Centroid X(um)' 'Centroid Y(um)' 'Pro. Area(um^2)' 'Perimeter(um)' 'A.R.' 'Shape Factor' 'PDI' 'Centre Mismatch' 'I80_by_I20' 'nHigh_by_nLow' 'Mean of Normalized Int' 'Median of Normalized Int' 'S.D of Normalized Int' 'Mode of Normalized Int' 'Entropy' 'Relative concavity'...
    'Num. of times curvature changes polarity' 'Avg Curvature(um^-1)' 'SD in curvature(um^-1)'...
    'Max positive curvature(um^-1)' 'Total positive curvature(um^-1)' 'Avg positive curvature(um^-1)' 'Total Positive Curvature>0.2(um^-1)' 'Avg Positive Curvature>0.2(um^-1)' 'Length with positive curvature>0.2(um)' 'Fraction of perimeter with positive curvature>0.2'...
    'Num positive curvatures > 0.25' 'Avg prominance of positive curvatures > 0.25(um)' 'Avg width of positive curvatures > 0.25(um)'...
    'Max negative curvature(um^-1)' 'Total negative curvature(um^-1)' 'Avg negative curvature(um^-1)' 'Total negative Curvature<-0.005(um^-1)' 'Avg negative Curvature<-0.005(um^-1)' 'Length with negative Curvature<-0.005(um^-1)' 'Fraction of perimeter with negative Curvature<-0.005(um^-1)'...
    'Num negative curvatures<-0.05' 'Avg prominance of negative curvatures<-0.05(um)' 'Avg width of negative curvatures<-0.05(um)'...
    'Average curvature x Perimeter'  'Avg positive curvature x Perimeter' 'Avg negative curvature x Perimeter'};


Result_Basic_Measurements=vertcat(header,Basic_Measurements);

clearvars condition_count file_count file_count1 upper_bound lower_bound num_of_nuc nuc_count Norm_Pix_Intensities Basic_Measurements header Boundary_RadOfCurva_K num_of_pix_in_largest_nuc number


% Lengthscales, Moment Invariants and Intensity Analyses
LengthScale_header={'Modal_LengthScale_FullNuc' 'Min_LengthScale_FullNuc' 'Avg_LengthScale_FullNuc' 'LengthScale_OfAvgCorr_FullNuc' ...
    'Modal_LengthScale_InCircle' 'Min_LengthScale_InCircle' 'Avg_LengthScale_InCircle' 'LengthScale_OfAvgCorr_InCircle' ...
    'SmallestNuclearSide' 'LargestNuclearSide' 'LargestToSmallestSide' 'avg_boundary_dist' 'SD_in_boundary_dist'};

MomentInvariants_header={'MomInv1' 'MomInv2' 'MomInv3' 'MomInv4' 'MomInv5' 'MomInv6' 'MomInv7' 'MomInv8'};

'Lengthscales and Moment Invariants'

CorrLengthScales=[];
MomentInv=[];
TotalCh1=[];
 for img_count=1:size(Selected_Segmented_2d_Raw,2)
        if ~isempty(Selected_Segmented_2d_Raw{img_count})
            CorrLengthScales_temp=[];
            MomentInv_temp=[];
            TotalCh1_temp=[];
            CorrLengthScales_temp(1,:)=Correlation_Lengthscales(Selected_Segmented_2d_Raw{img_count},MetaData.Voxel_Size_X);
            MomentInv_temp(1,:)=MomentInvariants(uint16(Selected_Segmented_2d_Raw{img_count}),MetaData.Voxel_Size_X);
            CorrLengthScales=vertcat(CorrLengthScales,CorrLengthScales_temp);
            MomentInv=vertcat(MomentInv,MomentInv_temp);
            
         end
    end


Result_Basic_Measurements=horzcat(Result_Basic_Measurements,...
        vertcat(LengthScale_header,num2cell(CorrLengthScales)),...
        vertcat(MomentInvariants_header,num2cell(MomentInv)));

cd 'E:\nih3t3_3d_cm\20190227_NIH3T3_3d_reprgrammed_spheroid_cm_aSMA_647_DAPI\CM\2D measures\'   
T = cell2table(Result_Basic_Measurements(1:end,:));
writetable(T,'kamal_pipeline.csv')
close all

%% 20190227_NIH3T3_3d_reprogrammed_spheroid_cm_asma_647_mkl_568_actin_488_dapi control
cd 'E:\nih3t3_3d_cm\20190227_NIH3T3_3d_reprogrammed_spheroid_cm_asma_647_mkl_568_actin_488_dapi\control\zproject\indivisual_nuclei_ch1'
clear all
sample='*';
filenames = vertcat(dir([sample,'.nd2']),dir([sample,'.tif']));

'loading images'

XYZ={};
MetaData={};

parfor reader_count=1:size(filenames,1);
    filename = filenames(reader_count).name;
    [XYZ{1,reader_count},MetaData{1,reader_count}]=imread(filename,1);
end
'Files imported'


cd 'E:\nih3t3_3d_cm\matlab_pipeline\'

MetaData=[];
MetaData.Voxel_Size_X = 0.21; % in µm
MetaData.Voxel_Size_Y =MetaData.Voxel_Size_X;
MetaData.Filename=filenames.name;

Selected_Segmented_2d_Raw=XYZ;
Basic_Measurements=cell(size(Selected_Segmented_2d_Raw));
Boundary_RadOfCurva_K=cell(size(XYZ));
Norm_Pix_Intensities=cell(size(Selected_Segmented_2d_Raw));
num_of_pix_in_largest_nuc=zeros(size(Selected_Segmented_2d_Raw));
num_of_nuc=zeros(size(Selected_Segmented_2d_Raw));
total_num_of_nuc=zeros(size(Selected_Segmented_2d_Raw,1),1);

'Basic Measurements'

[Basic_Measurements,Boundary_RadOfCurva_K,Norm_Pix_Intensities]=basic_measurements(Selected_Segmented_2d_Raw,MetaData);

header={'Filename' 'Nuc label' 'Centroid X(um)' 'Centroid Y(um)' 'Pro. Area(um^2)' 'Perimeter(um)' 'A.R.' 'Shape Factor' 'PDI' 'Centre Mismatch' 'I80_by_I20' 'nHigh_by_nLow' 'Mean of Normalized Int' 'Median of Normalized Int' 'S.D of Normalized Int' 'Mode of Normalized Int' 'Entropy' 'Relative concavity'...
    'Num. of times curvature changes polarity' 'Avg Curvature(um^-1)' 'SD in curvature(um^-1)'...
    'Max positive curvature(um^-1)' 'Total positive curvature(um^-1)' 'Avg positive curvature(um^-1)' 'Total Positive Curvature>0.2(um^-1)' 'Avg Positive Curvature>0.2(um^-1)' 'Length with positive curvature>0.2(um)' 'Fraction of perimeter with positive curvature>0.2'...
    'Num positive curvatures > 0.25' 'Avg prominance of positive curvatures > 0.25(um)' 'Avg width of positive curvatures > 0.25(um)'...
    'Max negative curvature(um^-1)' 'Total negative curvature(um^-1)' 'Avg negative curvature(um^-1)' 'Total negative Curvature<-0.005(um^-1)' 'Avg negative Curvature<-0.005(um^-1)' 'Length with negative Curvature<-0.005(um^-1)' 'Fraction of perimeter with negative Curvature<-0.005(um^-1)'...
    'Num negative curvatures<-0.05' 'Avg prominance of negative curvatures<-0.05(um)' 'Avg width of negative curvatures<-0.05(um)'...
    'Average curvature x Perimeter'  'Avg positive curvature x Perimeter' 'Avg negative curvature x Perimeter'};


Result_Basic_Measurements=vertcat(header,Basic_Measurements);

clearvars condition_count file_count file_count1 upper_bound lower_bound num_of_nuc nuc_count Norm_Pix_Intensities Basic_Measurements header Boundary_RadOfCurva_K num_of_pix_in_largest_nuc number


% Lengthscales, Moment Invariants and Intensity Analyses
LengthScale_header={'Modal_LengthScale_FullNuc' 'Min_LengthScale_FullNuc' 'Avg_LengthScale_FullNuc' 'LengthScale_OfAvgCorr_FullNuc' ...
    'Modal_LengthScale_InCircle' 'Min_LengthScale_InCircle' 'Avg_LengthScale_InCircle' 'LengthScale_OfAvgCorr_InCircle' ...
    'SmallestNuclearSide' 'LargestNuclearSide' 'LargestToSmallestSide' 'avg_boundary_dist' 'SD_in_boundary_dist'};

MomentInvariants_header={'MomInv1' 'MomInv2' 'MomInv3' 'MomInv4' 'MomInv5' 'MomInv6' 'MomInv7' 'MomInv8'};

'Lengthscales and Moment Invariants'

CorrLengthScales=[];
MomentInv=[];
TotalCh1=[];
 for img_count=1:size(Selected_Segmented_2d_Raw,2)
        if ~isempty(Selected_Segmented_2d_Raw{img_count})
            CorrLengthScales_temp=[];
            MomentInv_temp=[];
            TotalCh1_temp=[];
            CorrLengthScales_temp(1,:)=Correlation_Lengthscales(Selected_Segmented_2d_Raw{img_count},MetaData.Voxel_Size_X);
            MomentInv_temp(1,:)=MomentInvariants(uint16(Selected_Segmented_2d_Raw{img_count}),MetaData.Voxel_Size_X);
            CorrLengthScales=vertcat(CorrLengthScales,CorrLengthScales_temp);
            MomentInv=vertcat(MomentInv,MomentInv_temp);
            
         end
    end


Result_Basic_Measurements=horzcat(Result_Basic_Measurements,...
        vertcat(LengthScale_header,num2cell(CorrLengthScales)),...
        vertcat(MomentInvariants_header,num2cell(MomentInv)));

cd 'E:\nih3t3_3d_cm\20190227_NIH3T3_3d_reprogrammed_spheroid_cm_asma_647_mkl_568_actin_488_dapi\control\2D measures\'   
T = cell2table(Result_Basic_Measurements(1:end,:));
writetable(T,'kamal_pipeline.csv')

close all

%% 20190227_NIH3T3_3d_reprogrammed_spheroid_cm_asma_647_mkl_568_actin_488_dapi CM
cd 'E:\nih3t3_3d_cm\20190227_NIH3T3_3d_reprogrammed_spheroid_cm_asma_647_mkl_568_actin_488_dapi\CM\zproject\indivisual_nuclei_ch1'
clear all
sample='*';
filenames = vertcat(dir([sample,'.nd2']),dir([sample,'.tif']));

'loading images'

XYZ={};
MetaData={};

parfor reader_count=1:size(filenames,1);
    filename = filenames(reader_count).name;
    [XYZ{1,reader_count},MetaData{1,reader_count}]=imread(filename,1);
end
'Files imported'


cd 'E:\nih3t3_3d_cm\matlab_pipeline\'

MetaData=[];
MetaData.Voxel_Size_X = 0.21; % in µm
MetaData.Voxel_Size_Y =MetaData.Voxel_Size_X;
MetaData.Filename=filenames.name;

Selected_Segmented_2d_Raw=XYZ;
Basic_Measurements=cell(size(Selected_Segmented_2d_Raw));
Boundary_RadOfCurva_K=cell(size(XYZ));
Norm_Pix_Intensities=cell(size(Selected_Segmented_2d_Raw));
num_of_pix_in_largest_nuc=zeros(size(Selected_Segmented_2d_Raw));
num_of_nuc=zeros(size(Selected_Segmented_2d_Raw));
total_num_of_nuc=zeros(size(Selected_Segmented_2d_Raw,1),1);

'Basic Measurements'

[Basic_Measurements,Boundary_RadOfCurva_K,Norm_Pix_Intensities]=basic_measurements(Selected_Segmented_2d_Raw,MetaData);

header={'Filename' 'Nuc label' 'Centroid X(um)' 'Centroid Y(um)' 'Pro. Area(um^2)' 'Perimeter(um)' 'A.R.' 'Shape Factor' 'PDI' 'Centre Mismatch' 'I80_by_I20' 'nHigh_by_nLow' 'Mean of Normalized Int' 'Median of Normalized Int' 'S.D of Normalized Int' 'Mode of Normalized Int' 'Entropy' 'Relative concavity'...
    'Num. of times curvature changes polarity' 'Avg Curvature(um^-1)' 'SD in curvature(um^-1)'...
    'Max positive curvature(um^-1)' 'Total positive curvature(um^-1)' 'Avg positive curvature(um^-1)' 'Total Positive Curvature>0.2(um^-1)' 'Avg Positive Curvature>0.2(um^-1)' 'Length with positive curvature>0.2(um)' 'Fraction of perimeter with positive curvature>0.2'...
    'Num positive curvatures > 0.25' 'Avg prominance of positive curvatures > 0.25(um)' 'Avg width of positive curvatures > 0.25(um)'...
    'Max negative curvature(um^-1)' 'Total negative curvature(um^-1)' 'Avg negative curvature(um^-1)' 'Total negative Curvature<-0.005(um^-1)' 'Avg negative Curvature<-0.005(um^-1)' 'Length with negative Curvature<-0.005(um^-1)' 'Fraction of perimeter with negative Curvature<-0.005(um^-1)'...
    'Num negative curvatures<-0.05' 'Avg prominance of negative curvatures<-0.05(um)' 'Avg width of negative curvatures<-0.05(um)'...
    'Average curvature x Perimeter'  'Avg positive curvature x Perimeter' 'Avg negative curvature x Perimeter'};


Result_Basic_Measurements=vertcat(header,Basic_Measurements);

clearvars condition_count file_count file_count1 upper_bound lower_bound num_of_nuc nuc_count Norm_Pix_Intensities Basic_Measurements header Boundary_RadOfCurva_K num_of_pix_in_largest_nuc number


% Lengthscales, Moment Invariants and Intensity Analyses
LengthScale_header={'Modal_LengthScale_FullNuc' 'Min_LengthScale_FullNuc' 'Avg_LengthScale_FullNuc' 'LengthScale_OfAvgCorr_FullNuc' ...
    'Modal_LengthScale_InCircle' 'Min_LengthScale_InCircle' 'Avg_LengthScale_InCircle' 'LengthScale_OfAvgCorr_InCircle' ...
    'SmallestNuclearSide' 'LargestNuclearSide' 'LargestToSmallestSide' 'avg_boundary_dist' 'SD_in_boundary_dist'};

MomentInvariants_header={'MomInv1' 'MomInv2' 'MomInv3' 'MomInv4' 'MomInv5' 'MomInv6' 'MomInv7' 'MomInv8'};

'Lengthscales and Moment Invariants'

CorrLengthScales=[];
MomentInv=[];
TotalCh1=[];
 for img_count=1:size(Selected_Segmented_2d_Raw,2)
        if ~isempty(Selected_Segmented_2d_Raw{img_count})
            CorrLengthScales_temp=[];
            MomentInv_temp=[];
            TotalCh1_temp=[];
            CorrLengthScales_temp(1,:)=Correlation_Lengthscales(Selected_Segmented_2d_Raw{img_count},MetaData.Voxel_Size_X);
            MomentInv_temp(1,:)=MomentInvariants(uint16(Selected_Segmented_2d_Raw{img_count}),MetaData.Voxel_Size_X);
            CorrLengthScales=vertcat(CorrLengthScales,CorrLengthScales_temp);
            MomentInv=vertcat(MomentInv,MomentInv_temp);
            
         end
    end


Result_Basic_Measurements=horzcat(Result_Basic_Measurements,...
        vertcat(LengthScale_header,num2cell(CorrLengthScales)),...
        vertcat(MomentInvariants_header,num2cell(MomentInv)));

cd 'E:\nih3t3_3d_cm\20190227_NIH3T3_3d_reprogrammed_spheroid_cm_asma_647_mkl_568_actin_488_dapi\CM\2D measures\'   
T = cell2table(Result_Basic_Measurements(1:end,:));
writetable(T,'kamal_pipeline.csv')
close all




