% %% 20171130_NIH3t3_mcf7_co_culture_phalloidin_647_pol2_488_MKL_568_DAPI control
% cd 'E:\co_culture\20171130_NIH3t3_mcf7_co_culture_phalloidin_647_pol2_488_MKL_568_DAPI\control\zproject\indivisual_nuclei_ch1'
% clear all
% sample='*';
% filenames = vertcat(dir([sample,'.nd2']),dir([sample,'.tif']));
% 
% 'loading images'
% 
% XYZ={};
% MetaData={};
% 
% parfor reader_count=1:size(filenames,1);
%     filename = filenames(reader_count).name;
%     [XYZ{1,reader_count},MetaData{1,reader_count}]=imread(filename,1);
% end
% 'Files imported'
% 
% 
% cd 'E:\co_culture\matlab_pipeline\'
% 
% MetaData=[];
% MetaData.Voxel_Size_X = 0.21; % in µm
% MetaData.Voxel_Size_Y =MetaData.Voxel_Size_X;
% MetaData.Filename=filenames.name;
% 
% Selected_Segmented_2d_Raw=XYZ;
% Basic_Measurements=cell(size(Selected_Segmented_2d_Raw));
% Boundary_RadOfCurva_K=cell(size(XYZ));
% Norm_Pix_Intensities=cell(size(Selected_Segmented_2d_Raw));
% num_of_pix_in_largest_nuc=zeros(size(Selected_Segmented_2d_Raw));
% num_of_nuc=zeros(size(Selected_Segmented_2d_Raw));
% total_num_of_nuc=zeros(size(Selected_Segmented_2d_Raw,1),1);
% 
% 'Basic Measurements'
% 
% [Basic_Measurements,Boundary_RadOfCurva_K,Norm_Pix_Intensities]=basic_measurements(Selected_Segmented_2d_Raw,MetaData);
% 
% header={'Filename' 'Nuc label' 'Centroid X(um)' 'Centroid Y(um)' 'Pro. Area(um^2)' 'Perimeter(um)' 'A.R.' 'Shape Factor' 'PDI' 'Centre Mismatch' 'I80_by_I20' 'nHigh_by_nLow' 'Mean of Normalized Int' 'Median of Normalized Int' 'S.D of Normalized Int' 'Mode of Normalized Int' 'Entropy' 'Relative concavity'...
%     'Num. of times curvature changes polarity' 'Avg Curvature(um^-1)' 'SD in curvature(um^-1)'...
%     'Max positive curvature(um^-1)' 'Total positive curvature(um^-1)' 'Avg positive curvature(um^-1)' 'Total Positive Curvature>0.2(um^-1)' 'Avg Positive Curvature>0.2(um^-1)' 'Length with positive curvature>0.2(um)' 'Fraction of perimeter with positive curvature>0.2'...
%     'Num positive curvatures > 0.25' 'Avg prominance of positive curvatures > 0.25(um)' 'Avg width of positive curvatures > 0.25(um)'...
%     'Max negative curvature(um^-1)' 'Total negative curvature(um^-1)' 'Avg negative curvature(um^-1)' 'Total negative Curvature<-0.005(um^-1)' 'Avg negative Curvature<-0.005(um^-1)' 'Length with negative Curvature<-0.005(um^-1)' 'Fraction of perimeter with negative Curvature<-0.005(um^-1)'...
%     'Num negative curvatures<-0.05' 'Avg prominance of negative curvatures<-0.05(um)' 'Avg width of negative curvatures<-0.05(um)'...
%     'Average curvature x Perimeter'  'Avg positive curvature x Perimeter' 'Avg negative curvature x Perimeter'};
% 
% 
% Result_Basic_Measurements=vertcat(header,Basic_Measurements);
% 
% clearvars condition_count file_count file_count1 upper_bound lower_bound num_of_nuc nuc_count Norm_Pix_Intensities Basic_Measurements header Boundary_RadOfCurva_K num_of_pix_in_largest_nuc number
% 
% 
% % Lengthscales, Moment Invariants and Intensity Analyses
% LengthScale_header={'Modal_LengthScale_FullNuc' 'Min_LengthScale_FullNuc' 'Avg_LengthScale_FullNuc' 'LengthScale_OfAvgCorr_FullNuc' ...
%     'Modal_LengthScale_InCircle' 'Min_LengthScale_InCircle' 'Avg_LengthScale_InCircle' 'LengthScale_OfAvgCorr_InCircle' ...
%     'SmallestNuclearSide' 'LargestNuclearSide' 'LargestToSmallestSide' 'avg_boundary_dist' 'SD_in_boundary_dist'};
% 
% MomentInvariants_header={'MomInv1' 'MomInv2' 'MomInv3' 'MomInv4' 'MomInv5' 'MomInv6' 'MomInv7' 'MomInv8'};
% 
% 'Lengthscales and Moment Invariants'
% 
% CorrLengthScales=[];
% MomentInv=[];
% TotalCh1=[];
%  for img_count=1:size(Selected_Segmented_2d_Raw,2)
%         if ~isempty(Selected_Segmented_2d_Raw{img_count})
%             CorrLengthScales_temp=[];
%             MomentInv_temp=[];
%             TotalCh1_temp=[];
%             CorrLengthScales_temp(1,:)=Correlation_Lengthscales(Selected_Segmented_2d_Raw{img_count},MetaData.Voxel_Size_X);
%             MomentInv_temp(1,:)=MomentInvariants(uint16(Selected_Segmented_2d_Raw{img_count}),MetaData.Voxel_Size_X);
%             CorrLengthScales=vertcat(CorrLengthScales,CorrLengthScales_temp);
%             MomentInv=vertcat(MomentInv,MomentInv_temp);
%             
%          end
%     end
% 
% 
% Result_Basic_Measurements=horzcat(Result_Basic_Measurements,...
%         vertcat(LengthScale_header,num2cell(CorrLengthScales)),...
%         vertcat(MomentInvariants_header,num2cell(MomentInv)));
% 
% cd 'E:\co_culture\20171130_NIH3t3_mcf7_co_culture_phalloidin_647_pol2_488_MKL_568_DAPI\control\2D measures\'   
% T = cell2table(Result_Basic_Measurements(1:end,:));
% writetable(T,'kamal_pipeline.csv')
% 
% %% 20171130_NIH3t3_mcf7_co_culture_phalloidin_647_pol2_488_MKL_568_DAPI co_culture
% cd 'E:\co_culture\20171130_NIH3t3_mcf7_co_culture_phalloidin_647_pol2_488_MKL_568_DAPI\co_culture\zproject\indivisual_nuclei_ch1'
% clear all
% sample='*';
% filenames = vertcat(dir([sample,'.nd2']),dir([sample,'.tif']));
% 
% 'loading images'
% 
% XYZ={};
% MetaData={};
% 
% parfor reader_count=1:size(filenames,1);
%     filename = filenames(reader_count).name;
%     [XYZ{1,reader_count},MetaData{1,reader_count}]=imread(filename,1);
% end
% 'Files imported'
% 
% 
% cd 'E:\co_culture\matlab_pipeline\'
% 
% MetaData=[];
% MetaData.Voxel_Size_X = 0.21; % in µm
% MetaData.Voxel_Size_Y =MetaData.Voxel_Size_X;
% MetaData.Filename=filenames.name;
% 
% Selected_Segmented_2d_Raw=XYZ;
% Basic_Measurements=cell(size(Selected_Segmented_2d_Raw));
% Boundary_RadOfCurva_K=cell(size(XYZ));
% Norm_Pix_Intensities=cell(size(Selected_Segmented_2d_Raw));
% num_of_pix_in_largest_nuc=zeros(size(Selected_Segmented_2d_Raw));
% num_of_nuc=zeros(size(Selected_Segmented_2d_Raw));
% total_num_of_nuc=zeros(size(Selected_Segmented_2d_Raw,1),1);
% 
% 'Basic Measurements'
% 
% [Basic_Measurements,Boundary_RadOfCurva_K,Norm_Pix_Intensities]=basic_measurements(Selected_Segmented_2d_Raw,MetaData);
% 
% header={'Filename' 'Nuc label' 'Centroid X(um)' 'Centroid Y(um)' 'Pro. Area(um^2)' 'Perimeter(um)' 'A.R.' 'Shape Factor' 'PDI' 'Centre Mismatch' 'I80_by_I20' 'nHigh_by_nLow' 'Mean of Normalized Int' 'Median of Normalized Int' 'S.D of Normalized Int' 'Mode of Normalized Int' 'Entropy' 'Relative concavity'...
%     'Num. of times curvature changes polarity' 'Avg Curvature(um^-1)' 'SD in curvature(um^-1)'...
%     'Max positive curvature(um^-1)' 'Total positive curvature(um^-1)' 'Avg positive curvature(um^-1)' 'Total Positive Curvature>0.2(um^-1)' 'Avg Positive Curvature>0.2(um^-1)' 'Length with positive curvature>0.2(um)' 'Fraction of perimeter with positive curvature>0.2'...
%     'Num positive curvatures > 0.25' 'Avg prominance of positive curvatures > 0.25(um)' 'Avg width of positive curvatures > 0.25(um)'...
%     'Max negative curvature(um^-1)' 'Total negative curvature(um^-1)' 'Avg negative curvature(um^-1)' 'Total negative Curvature<-0.005(um^-1)' 'Avg negative Curvature<-0.005(um^-1)' 'Length with negative Curvature<-0.005(um^-1)' 'Fraction of perimeter with negative Curvature<-0.005(um^-1)'...
%     'Num negative curvatures<-0.05' 'Avg prominance of negative curvatures<-0.05(um)' 'Avg width of negative curvatures<-0.05(um)'...
%     'Average curvature x Perimeter'  'Avg positive curvature x Perimeter' 'Avg negative curvature x Perimeter'};
% 
% 
% Result_Basic_Measurements=vertcat(header,Basic_Measurements);
% 
% clearvars condition_count file_count file_count1 upper_bound lower_bound num_of_nuc nuc_count Norm_Pix_Intensities Basic_Measurements header Boundary_RadOfCurva_K num_of_pix_in_largest_nuc number
% 
% 
% % Lengthscales, Moment Invariants and Intensity Analyses
% LengthScale_header={'Modal_LengthScale_FullNuc' 'Min_LengthScale_FullNuc' 'Avg_LengthScale_FullNuc' 'LengthScale_OfAvgCorr_FullNuc' ...
%     'Modal_LengthScale_InCircle' 'Min_LengthScale_InCircle' 'Avg_LengthScale_InCircle' 'LengthScale_OfAvgCorr_InCircle' ...
%     'SmallestNuclearSide' 'LargestNuclearSide' 'LargestToSmallestSide' 'avg_boundary_dist' 'SD_in_boundary_dist'};
% 
% MomentInvariants_header={'MomInv1' 'MomInv2' 'MomInv3' 'MomInv4' 'MomInv5' 'MomInv6' 'MomInv7' 'MomInv8'};
% 
% 'Lengthscales and Moment Invariants'
% 
% CorrLengthScales=[];
% MomentInv=[];
% TotalCh1=[];
%  for img_count=1:size(Selected_Segmented_2d_Raw,2)
%         if ~isempty(Selected_Segmented_2d_Raw{img_count})
%             CorrLengthScales_temp=[];
%             MomentInv_temp=[];
%             TotalCh1_temp=[];
%             CorrLengthScales_temp(1,:)=Correlation_Lengthscales(Selected_Segmented_2d_Raw{img_count},MetaData.Voxel_Size_X);
%             MomentInv_temp(1,:)=MomentInvariants(uint16(Selected_Segmented_2d_Raw{img_count}),MetaData.Voxel_Size_X);
%             CorrLengthScales=vertcat(CorrLengthScales,CorrLengthScales_temp);
%             MomentInv=vertcat(MomentInv,MomentInv_temp);
%             
%          end
%     end
% 
% 
% Result_Basic_Measurements=horzcat(Result_Basic_Measurements,...
%         vertcat(LengthScale_header,num2cell(CorrLengthScales)),...
%         vertcat(MomentInvariants_header,num2cell(MomentInv)));
% 
% cd 'E:\co_culture\20171130_NIH3t3_mcf7_co_culture_phalloidin_647_pol2_488_MKL_568_DAPI\co_culture\2D measures\'   
% T = cell2table(Result_Basic_Measurements(1:end,:));
% writetable(T,'kamal_pipeline.csv')

% %% 20180310_MKL488_aSMA_568_DAPI_EpCAM_647 control
% cd 'E:\co_culture\20180310_MKL488_aSMA_568_DAPI_EpCAM_647\control\zproject\indivisual_nuclei_ch1'
% clear all
% sample='*';
% filenames = vertcat(dir([sample,'.nd2']),dir([sample,'.tif']));
% 
% 'loading images'
% 
% XYZ={};
% MetaData={};
% 
% parfor reader_count=1:size(filenames,1);
%     filename = filenames(reader_count).name;
%     [XYZ{1,reader_count},MetaData{1,reader_count}]=imread(filename,1);
% end
% 'Files imported'
% 
% 
% cd 'E:\co_culture\matlab_pipeline\'
% 
% MetaData=[];
% MetaData.Voxel_Size_X = 0.21; % in µm
% MetaData.Voxel_Size_Y =MetaData.Voxel_Size_X;
% MetaData.Filename=filenames.name;
% 
% Selected_Segmented_2d_Raw=XYZ;
% Basic_Measurements=cell(size(Selected_Segmented_2d_Raw));
% Boundary_RadOfCurva_K=cell(size(XYZ));
% Norm_Pix_Intensities=cell(size(Selected_Segmented_2d_Raw));
% num_of_pix_in_largest_nuc=zeros(size(Selected_Segmented_2d_Raw));
% num_of_nuc=zeros(size(Selected_Segmented_2d_Raw));
% total_num_of_nuc=zeros(size(Selected_Segmented_2d_Raw,1),1);
% 
% 'Basic Measurements'
% 
% [Basic_Measurements,Boundary_RadOfCurva_K,Norm_Pix_Intensities]=basic_measurements(Selected_Segmented_2d_Raw,MetaData);
% 
% header={'Filename' 'Nuc label' 'Centroid X(um)' 'Centroid Y(um)' 'Pro. Area(um^2)' 'Perimeter(um)' 'A.R.' 'Shape Factor' 'PDI' 'Centre Mismatch' 'I80_by_I20' 'nHigh_by_nLow' 'Mean of Normalized Int' 'Median of Normalized Int' 'S.D of Normalized Int' 'Mode of Normalized Int' 'Entropy' 'Relative concavity'...
%     'Num. of times curvature changes polarity' 'Avg Curvature(um^-1)' 'SD in curvature(um^-1)'...
%     'Max positive curvature(um^-1)' 'Total positive curvature(um^-1)' 'Avg positive curvature(um^-1)' 'Total Positive Curvature>0.2(um^-1)' 'Avg Positive Curvature>0.2(um^-1)' 'Length with positive curvature>0.2(um)' 'Fraction of perimeter with positive curvature>0.2'...
%     'Num positive curvatures > 0.25' 'Avg prominance of positive curvatures > 0.25(um)' 'Avg width of positive curvatures > 0.25(um)'...
%     'Max negative curvature(um^-1)' 'Total negative curvature(um^-1)' 'Avg negative curvature(um^-1)' 'Total negative Curvature<-0.005(um^-1)' 'Avg negative Curvature<-0.005(um^-1)' 'Length with negative Curvature<-0.005(um^-1)' 'Fraction of perimeter with negative Curvature<-0.005(um^-1)'...
%     'Num negative curvatures<-0.05' 'Avg prominance of negative curvatures<-0.05(um)' 'Avg width of negative curvatures<-0.05(um)'...
%     'Average curvature x Perimeter'  'Avg positive curvature x Perimeter' 'Avg negative curvature x Perimeter'};
% 
% 
% Result_Basic_Measurements=vertcat(header,Basic_Measurements);
% 
% clearvars condition_count file_count file_count1 upper_bound lower_bound num_of_nuc nuc_count Norm_Pix_Intensities Basic_Measurements header Boundary_RadOfCurva_K num_of_pix_in_largest_nuc number
% 
% 
% % Lengthscales, Moment Invariants and Intensity Analyses
% LengthScale_header={'Modal_LengthScale_FullNuc' 'Min_LengthScale_FullNuc' 'Avg_LengthScale_FullNuc' 'LengthScale_OfAvgCorr_FullNuc' ...
%     'Modal_LengthScale_InCircle' 'Min_LengthScale_InCircle' 'Avg_LengthScale_InCircle' 'LengthScale_OfAvgCorr_InCircle' ...
%     'SmallestNuclearSide' 'LargestNuclearSide' 'LargestToSmallestSide' 'avg_boundary_dist' 'SD_in_boundary_dist'};
% 
% MomentInvariants_header={'MomInv1' 'MomInv2' 'MomInv3' 'MomInv4' 'MomInv5' 'MomInv6' 'MomInv7' 'MomInv8'};
% 
% 'Lengthscales and Moment Invariants'
% 
% CorrLengthScales=[];
% MomentInv=[];
% TotalCh1=[];
%  for img_count=1:size(Selected_Segmented_2d_Raw,2)
%         if ~isempty(Selected_Segmented_2d_Raw{img_count})
%             CorrLengthScales_temp=[];
%             MomentInv_temp=[];
%             TotalCh1_temp=[];
%             CorrLengthScales_temp(1,:)=Correlation_Lengthscales(Selected_Segmented_2d_Raw{img_count},MetaData.Voxel_Size_X);
%             MomentInv_temp(1,:)=MomentInvariants(uint16(Selected_Segmented_2d_Raw{img_count}),MetaData.Voxel_Size_X);
%             CorrLengthScales=vertcat(CorrLengthScales,CorrLengthScales_temp);
%             MomentInv=vertcat(MomentInv,MomentInv_temp);
%             
%          end
%     end
% 
% 
% Result_Basic_Measurements=horzcat(Result_Basic_Measurements,...
%         vertcat(LengthScale_header,num2cell(CorrLengthScales)),...
%         vertcat(MomentInvariants_header,num2cell(MomentInv)));
% 
% cd 'E:\co_culture\20180310_MKL488_aSMA_568_DAPI_EpCAM_647\control\2D measures\'   
% T = cell2table(Result_Basic_Measurements(1:end,:));
% writetable(T,'kamal_pipeline.csv')
% 
% %% 20180310_MKL488_aSMA_568_DAPI_EpCAM_647 co_culture
% cd 'E:\co_culture\20180310_MKL488_aSMA_568_DAPI_EpCAM_647\co_culture\zproject\indivisual_nuclei_ch1'
% clear all
% sample='*';
% filenames = vertcat(dir([sample,'.nd2']),dir([sample,'.tif']));
% 
% 'loading images'
% 
% XYZ={};
% MetaData={};
% 
% parfor reader_count=1:size(filenames,1);
%     filename = filenames(reader_count).name;
%     [XYZ{1,reader_count},MetaData{1,reader_count}]=imread(filename,1);
% end
% 'Files imported'
% 
% 
% cd 'E:\co_culture\matlab_pipeline\'
% 
% MetaData=[];
% MetaData.Voxel_Size_X = 0.21; % in µm
% MetaData.Voxel_Size_Y =MetaData.Voxel_Size_X;
% MetaData.Filename=filenames.name;
% 
% Selected_Segmented_2d_Raw=XYZ;
% Basic_Measurements=cell(size(Selected_Segmented_2d_Raw));
% Boundary_RadOfCurva_K=cell(size(XYZ));
% Norm_Pix_Intensities=cell(size(Selected_Segmented_2d_Raw));
% num_of_pix_in_largest_nuc=zeros(size(Selected_Segmented_2d_Raw));
% num_of_nuc=zeros(size(Selected_Segmented_2d_Raw));
% total_num_of_nuc=zeros(size(Selected_Segmented_2d_Raw,1),1);
% 
% 'Basic Measurements'
% 
% [Basic_Measurements,Boundary_RadOfCurva_K,Norm_Pix_Intensities]=basic_measurements(Selected_Segmented_2d_Raw,MetaData);
% 
% header={'Filename' 'Nuc label' 'Centroid X(um)' 'Centroid Y(um)' 'Pro. Area(um^2)' 'Perimeter(um)' 'A.R.' 'Shape Factor' 'PDI' 'Centre Mismatch' 'I80_by_I20' 'nHigh_by_nLow' 'Mean of Normalized Int' 'Median of Normalized Int' 'S.D of Normalized Int' 'Mode of Normalized Int' 'Entropy' 'Relative concavity'...
%     'Num. of times curvature changes polarity' 'Avg Curvature(um^-1)' 'SD in curvature(um^-1)'...
%     'Max positive curvature(um^-1)' 'Total positive curvature(um^-1)' 'Avg positive curvature(um^-1)' 'Total Positive Curvature>0.2(um^-1)' 'Avg Positive Curvature>0.2(um^-1)' 'Length with positive curvature>0.2(um)' 'Fraction of perimeter with positive curvature>0.2'...
%     'Num positive curvatures > 0.25' 'Avg prominance of positive curvatures > 0.25(um)' 'Avg width of positive curvatures > 0.25(um)'...
%     'Max negative curvature(um^-1)' 'Total negative curvature(um^-1)' 'Avg negative curvature(um^-1)' 'Total negative Curvature<-0.005(um^-1)' 'Avg negative Curvature<-0.005(um^-1)' 'Length with negative Curvature<-0.005(um^-1)' 'Fraction of perimeter with negative Curvature<-0.005(um^-1)'...
%     'Num negative curvatures<-0.05' 'Avg prominance of negative curvatures<-0.05(um)' 'Avg width of negative curvatures<-0.05(um)'...
%     'Average curvature x Perimeter'  'Avg positive curvature x Perimeter' 'Avg negative curvature x Perimeter'};
% 
% 
% Result_Basic_Measurements=vertcat(header,Basic_Measurements);
% 
% clearvars condition_count file_count file_count1 upper_bound lower_bound num_of_nuc nuc_count Norm_Pix_Intensities Basic_Measurements header Boundary_RadOfCurva_K num_of_pix_in_largest_nuc number
% 
% 
% % Lengthscales, Moment Invariants and Intensity Analyses
% LengthScale_header={'Modal_LengthScale_FullNuc' 'Min_LengthScale_FullNuc' 'Avg_LengthScale_FullNuc' 'LengthScale_OfAvgCorr_FullNuc' ...
%     'Modal_LengthScale_InCircle' 'Min_LengthScale_InCircle' 'Avg_LengthScale_InCircle' 'LengthScale_OfAvgCorr_InCircle' ...
%     'SmallestNuclearSide' 'LargestNuclearSide' 'LargestToSmallestSide' 'avg_boundary_dist' 'SD_in_boundary_dist'};
% 
% MomentInvariants_header={'MomInv1' 'MomInv2' 'MomInv3' 'MomInv4' 'MomInv5' 'MomInv6' 'MomInv7' 'MomInv8'};
% 
% 'Lengthscales and Moment Invariants'
% 
% CorrLengthScales=[];
% MomentInv=[];
% TotalCh1=[];
%  for img_count=1:size(Selected_Segmented_2d_Raw,2)
%         if ~isempty(Selected_Segmented_2d_Raw{img_count})
%             CorrLengthScales_temp=[];
%             MomentInv_temp=[];
%             TotalCh1_temp=[];
%             CorrLengthScales_temp(1,:)=Correlation_Lengthscales(Selected_Segmented_2d_Raw{img_count},MetaData.Voxel_Size_X);
%             MomentInv_temp(1,:)=MomentInvariants(uint16(Selected_Segmented_2d_Raw{img_count}),MetaData.Voxel_Size_X);
%             CorrLengthScales=vertcat(CorrLengthScales,CorrLengthScales_temp);
%             MomentInv=vertcat(MomentInv,MomentInv_temp);
%             
%          end
%     end
% 
% 
% Result_Basic_Measurements=horzcat(Result_Basic_Measurements,...
%         vertcat(LengthScale_header,num2cell(CorrLengthScales)),...
%         vertcat(MomentInvariants_header,num2cell(MomentInv)));
% 
% cd 'E:\co_culture\20180310_MKL488_aSMA_568_DAPI_EpCAM_647\co_culture\2D measures\'   
% T = cell2table(Result_Basic_Measurements(1:end,:));
% writetable(T,'kamal_pipeline.csv')

% %% 20180811_MCF7_3T3_co_culture_3Days_MKL_gt__647_pol2_rb_568_mcf7_488 control
% cd 'E:\co_culture\20180811_MCF7_3T3_co_culture_3Days_MKL_gt__647_pol2_rb_568_mcf7_488\control\zproject\indivisual_nuclei_ch1'
% clear all
% sample='*';
% filenames = vertcat(dir([sample,'.nd2']),dir([sample,'.tif']));
% 
% 'loading images'
% 
% XYZ={};
% MetaData={};
% 
% parfor reader_count=1:size(filenames,1);
%     filename = filenames(reader_count).name;
%     [XYZ{1,reader_count},MetaData{1,reader_count}]=imread(filename,1);
% end
% 'Files imported'
% 
% 
% cd 'E:\co_culture\matlab_pipeline\'
% 
% MetaData=[];
% MetaData.Voxel_Size_X = 0.21; % in µm
% MetaData.Voxel_Size_Y =MetaData.Voxel_Size_X;
% MetaData.Filename=filenames.name;
% 
% Selected_Segmented_2d_Raw=XYZ;
% Basic_Measurements=cell(size(Selected_Segmented_2d_Raw));
% Boundary_RadOfCurva_K=cell(size(XYZ));
% Norm_Pix_Intensities=cell(size(Selected_Segmented_2d_Raw));
% num_of_pix_in_largest_nuc=zeros(size(Selected_Segmented_2d_Raw));
% num_of_nuc=zeros(size(Selected_Segmented_2d_Raw));
% total_num_of_nuc=zeros(size(Selected_Segmented_2d_Raw,1),1);
% 
% 'Basic Measurements'
% 
% [Basic_Measurements,Boundary_RadOfCurva_K,Norm_Pix_Intensities]=basic_measurements(Selected_Segmented_2d_Raw,MetaData);
% 
% header={'Filename' 'Nuc label' 'Centroid X(um)' 'Centroid Y(um)' 'Pro. Area(um^2)' 'Perimeter(um)' 'A.R.' 'Shape Factor' 'PDI' 'Centre Mismatch' 'I80_by_I20' 'nHigh_by_nLow' 'Mean of Normalized Int' 'Median of Normalized Int' 'S.D of Normalized Int' 'Mode of Normalized Int' 'Entropy' 'Relative concavity'...
%     'Num. of times curvature changes polarity' 'Avg Curvature(um^-1)' 'SD in curvature(um^-1)'...
%     'Max positive curvature(um^-1)' 'Total positive curvature(um^-1)' 'Avg positive curvature(um^-1)' 'Total Positive Curvature>0.2(um^-1)' 'Avg Positive Curvature>0.2(um^-1)' 'Length with positive curvature>0.2(um)' 'Fraction of perimeter with positive curvature>0.2'...
%     'Num positive curvatures > 0.25' 'Avg prominance of positive curvatures > 0.25(um)' 'Avg width of positive curvatures > 0.25(um)'...
%     'Max negative curvature(um^-1)' 'Total negative curvature(um^-1)' 'Avg negative curvature(um^-1)' 'Total negative Curvature<-0.005(um^-1)' 'Avg negative Curvature<-0.005(um^-1)' 'Length with negative Curvature<-0.005(um^-1)' 'Fraction of perimeter with negative Curvature<-0.005(um^-1)'...
%     'Num negative curvatures<-0.05' 'Avg prominance of negative curvatures<-0.05(um)' 'Avg width of negative curvatures<-0.05(um)'...
%     'Average curvature x Perimeter'  'Avg positive curvature x Perimeter' 'Avg negative curvature x Perimeter'};
% 
% 
% Result_Basic_Measurements=vertcat(header,Basic_Measurements);
% 
% clearvars condition_count file_count file_count1 upper_bound lower_bound num_of_nuc nuc_count Norm_Pix_Intensities Basic_Measurements header Boundary_RadOfCurva_K num_of_pix_in_largest_nuc number
% 
% 
% % Lengthscales, Moment Invariants and Intensity Analyses
% LengthScale_header={'Modal_LengthScale_FullNuc' 'Min_LengthScale_FullNuc' 'Avg_LengthScale_FullNuc' 'LengthScale_OfAvgCorr_FullNuc' ...
%     'Modal_LengthScale_InCircle' 'Min_LengthScale_InCircle' 'Avg_LengthScale_InCircle' 'LengthScale_OfAvgCorr_InCircle' ...
%     'SmallestNuclearSide' 'LargestNuclearSide' 'LargestToSmallestSide' 'avg_boundary_dist' 'SD_in_boundary_dist'};
% 
% MomentInvariants_header={'MomInv1' 'MomInv2' 'MomInv3' 'MomInv4' 'MomInv5' 'MomInv6' 'MomInv7' 'MomInv8'};
% 
% 'Lengthscales and Moment Invariants'
% 
% CorrLengthScales=[];
% MomentInv=[];
% TotalCh1=[];
%  for img_count=1:size(Selected_Segmented_2d_Raw,2)
%         if ~isempty(Selected_Segmented_2d_Raw{img_count})
%             CorrLengthScales_temp=[];
%             MomentInv_temp=[];
%             TotalCh1_temp=[];
%             CorrLengthScales_temp(1,:)=Correlation_Lengthscales(Selected_Segmented_2d_Raw{img_count},MetaData.Voxel_Size_X);
%             MomentInv_temp(1,:)=MomentInvariants(uint16(Selected_Segmented_2d_Raw{img_count}),MetaData.Voxel_Size_X);
%             CorrLengthScales=vertcat(CorrLengthScales,CorrLengthScales_temp);
%             MomentInv=vertcat(MomentInv,MomentInv_temp);
%             
%          end
%     end
% 
% 
% Result_Basic_Measurements=horzcat(Result_Basic_Measurements,...
%         vertcat(LengthScale_header,num2cell(CorrLengthScales)),...
%         vertcat(MomentInvariants_header,num2cell(MomentInv)));
% 
% cd 'E:\co_culture\20180811_MCF7_3T3_co_culture_3Days_MKL_gt__647_pol2_rb_568_mcf7_488\control\2D measures\'   
% T = cell2table(Result_Basic_Measurements(1:end,:));
% writetable(T,'kamal_pipeline.csv')

% %% 20180811_MCF7_3T3_co_culture_3Days_MKL_gt__647_pol2_rb_568_mcf7_488 co_culture
% cd 'E:\co_culture\20180811_MCF7_3T3_co_culture_3Days_MKL_gt__647_pol2_rb_568_mcf7_488\co_culture\zproject\indivisual_nuclei_ch1'
% clear all
% sample='*';
% filenames = vertcat(dir([sample,'.nd2']),dir([sample,'.tif']));
% 
% 'loading images'
% 
% XYZ={};
% MetaData={};
% 
% parfor reader_count=1:size(filenames,1);
%     filename = filenames(reader_count).name;
%     [XYZ{1,reader_count},MetaData{1,reader_count}]=imread(filename,1);
% end
% 'Files imported'
% 
% 
% cd 'E:\co_culture\matlab_pipeline\'
% 
% MetaData=[];
% MetaData.Voxel_Size_X = 0.21; % in µm
% MetaData.Voxel_Size_Y =MetaData.Voxel_Size_X;
% MetaData.Filename=filenames.name;
% 
% Selected_Segmented_2d_Raw=XYZ;
% Basic_Measurements=cell(size(Selected_Segmented_2d_Raw));
% Boundary_RadOfCurva_K=cell(size(XYZ));
% Norm_Pix_Intensities=cell(size(Selected_Segmented_2d_Raw));
% num_of_pix_in_largest_nuc=zeros(size(Selected_Segmented_2d_Raw));
% num_of_nuc=zeros(size(Selected_Segmented_2d_Raw));
% total_num_of_nuc=zeros(size(Selected_Segmented_2d_Raw,1),1);
% 
% 'Basic Measurements'
% 
% [Basic_Measurements,Boundary_RadOfCurva_K,Norm_Pix_Intensities]=basic_measurements(Selected_Segmented_2d_Raw,MetaData);
% 
% header={'Filename' 'Nuc label' 'Centroid X(um)' 'Centroid Y(um)' 'Pro. Area(um^2)' 'Perimeter(um)' 'A.R.' 'Shape Factor' 'PDI' 'Centre Mismatch' 'I80_by_I20' 'nHigh_by_nLow' 'Mean of Normalized Int' 'Median of Normalized Int' 'S.D of Normalized Int' 'Mode of Normalized Int' 'Entropy' 'Relative concavity'...
%     'Num. of times curvature changes polarity' 'Avg Curvature(um^-1)' 'SD in curvature(um^-1)'...
%     'Max positive curvature(um^-1)' 'Total positive curvature(um^-1)' 'Avg positive curvature(um^-1)' 'Total Positive Curvature>0.2(um^-1)' 'Avg Positive Curvature>0.2(um^-1)' 'Length with positive curvature>0.2(um)' 'Fraction of perimeter with positive curvature>0.2'...
%     'Num positive curvatures > 0.25' 'Avg prominance of positive curvatures > 0.25(um)' 'Avg width of positive curvatures > 0.25(um)'...
%     'Max negative curvature(um^-1)' 'Total negative curvature(um^-1)' 'Avg negative curvature(um^-1)' 'Total negative Curvature<-0.005(um^-1)' 'Avg negative Curvature<-0.005(um^-1)' 'Length with negative Curvature<-0.005(um^-1)' 'Fraction of perimeter with negative Curvature<-0.005(um^-1)'...
%     'Num negative curvatures<-0.05' 'Avg prominance of negative curvatures<-0.05(um)' 'Avg width of negative curvatures<-0.05(um)'...
%     'Average curvature x Perimeter'  'Avg positive curvature x Perimeter' 'Avg negative curvature x Perimeter'};
% 
% 
% Result_Basic_Measurements=vertcat(header,Basic_Measurements);
% 
% clearvars condition_count file_count file_count1 upper_bound lower_bound num_of_nuc nuc_count Norm_Pix_Intensities Basic_Measurements header Boundary_RadOfCurva_K num_of_pix_in_largest_nuc number
% 

%Lengthscales, Moment Invariants and Intensity Analyses
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

cd 'E:\co_culture\20180811_MCF7_3T3_co_culture_3Days_MKL_gt__647_pol2_rb_568_mcf7_488\co_culture\2D measures\'   
T = cell2table(Result_Basic_Measurements(1:end,:));
writetable(T,'kamal_pipeline.csv')

%% 20180812_MCF7_3T3_co_culture_3Days_MKL_rb_568_asma_ms_647_mcf7_488 control
cd 'E:\co_culture\20180812_MCF7_3T3_co_culture_3Days_MKL_rb_568_asma_ms_647_mcf7_488\control\zproject\indivisual_nuclei_ch1'
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


cd 'E:\co_culture\matlab_pipeline\'

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

cd 'E:\co_culture\20180812_MCF7_3T3_co_culture_3Days_MKL_rb_568_asma_ms_647_mcf7_488\control\2D measures\'   
T = cell2table(Result_Basic_Measurements(1:end,:));
writetable(T,'kamal_pipeline.csv')

%% 20180812_MCF7_3T3_co_culture_3Days_MKL_rb_568_asma_ms_647_mcf7_488 co_culture
cd 'E:\co_culture\20180812_MCF7_3T3_co_culture_3Days_MKL_rb_568_asma_ms_647_mcf7_488\co_culture\zproject\indivisual_nuclei_ch1'
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


cd 'E:\co_culture\matlab_pipeline\'

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

cd 'E:\co_culture\20180812_MCF7_3T3_co_culture_3Days_MKL_rb_568_asma_ms_647_mcf7_488\co_culture\2D measures\'   
T = cell2table(Result_Basic_Measurements(1:end,:));
writetable(T,'kamal_pipeline.csv')

%% 20180830_MCF7_3T3_co_culture_3Days_MKLrb_568_asma_ms_647_mcf7_488 control
cd 'E:\co_culture\20180830_MCF7_3T3_co_culture_3Days_MKLrb_568_asma_ms_647_mcf7_488\control\zproject\indivisual_nuclei_ch1'
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


cd 'E:\co_culture\matlab_pipeline\'

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

cd 'E:\co_culture\20180830_MCF7_3T3_co_culture_3Days_MKLrb_568_asma_ms_647_mcf7_488\control\2D measures\'   
T = cell2table(Result_Basic_Measurements(1:end,:));
writetable(T,'kamal_pipeline.csv')

%% 20180830_MCF7_3T3_co_culture_3Days_MKLrb_568_asma_ms_647_mcf7_488 co_culture
cd 'E:\co_culture\20180830_MCF7_3T3_co_culture_3Days_MKLrb_568_asma_ms_647_mcf7_488\co_culture\zproject\indivisual_nuclei_ch1'
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


cd 'E:\co_culture\matlab_pipeline\'

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

cd 'E:\co_culture\20180830_MCF7_3T3_co_culture_3Days_MKLrb_568_asma_ms_647_mcf7_488\co_culture\2D measures\'   
T = cell2table(Result_Basic_Measurements(1:end,:));
writetable(T,'kamal_pipeline.csv')

%% 20181220_MCF7_3T3_Co_culture_3Days_MKLrb_568_asma_ms_647_mcf7_488 control
cd 'E:\co_culture\20181220_MCF7_3T3_Co_culture_3Days_MKLrb_568_asma_ms_647_mcf7_488\control\zproject\indivisual_nuclei_ch1'
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


cd 'E:\co_culture\matlab_pipeline\'

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

cd 'E:\co_culture\20181220_MCF7_3T3_Co_culture_3Days_MKLrb_568_asma_ms_647_mcf7_488\control\2D measures\'   
T = cell2table(Result_Basic_Measurements(1:end,:));
writetable(T,'kamal_pipeline.csv')

%% 20181220_MCF7_3T3_Co_culture_3Days_MKLrb_568_asma_ms_647_mcf7_488 co_culture
cd 'E:\co_culture\20181220_MCF7_3T3_Co_culture_3Days_MKLrb_568_asma_ms_647_mcf7_488\co_culture\zproject\indivisual_nuclei_ch1'
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


cd 'E:\co_culture\matlab_pipeline\'

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

cd 'E:\co_culture\20181220_MCF7_3T3_Co_culture_3Days_MKLrb_568_asma_ms_647_mcf7_488\co_culture\2D measures\'   
T = cell2table(Result_Basic_Measurements(1:end,:));
writetable(T,'kamal_pipeline.csv')
