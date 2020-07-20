%%% written by D.S.JOKHUN on 10/05/2018


clear all


sample='*'

filenames = vertcat(dir([sample,'.nd2']),dir([sample,'.tif']));





%%
'loading images'
tic
XYZ={};
XYZ_Ch2={};
XYZ_Ch3={};
XYZ_Ch4={};
MetaData={};

parfor reader_count=1:size(filenames,1);
    filename = filenames(reader_count).name;
    [XYZ{1,reader_count},MetaData{1,reader_count}]=img_reader_3D(filename,1);
    [XYZ_Ch2{1,reader_count}]=img_reader_3D(filename,2);
    [XYZ_Ch3{1,reader_count}]=img_reader_3D(filename,3);
    [XYZ_Ch4{1,reader_count}]=img_reader_3D(filename,4);
end
'Files imported'

toc
%%




%%

% 'Segmentation, Flattening and Cropping'
% tic
% Raw_cropped_2D=cell(size(XYZ));
% Segmented_2D=cell(size(XYZ));
% parfor file_count=1:numel(XYZ)
%     file_count
%     [Raw_cropped_2D{file_count}, Segmented_2D{file_count}]=SegFlattenCrop3D(XYZ{file_count},MetaData{file_count})
% end
% '2D segmentation complete'



'Segmentation, Flattening and Cropping'
tic
Raw_cropped_2D=cell(size(XYZ));
Segmented_2D=cell(size(XYZ));
Segmented_2D_Ch2=cell(size(XYZ));
Segmented_2D_Ch3=cell(size(XYZ));
Segmented_2D_Ch4=cell(size(XYZ));

parfor file_count=1:numel(XYZ)
    file_count
    [Raw_cropped_2D{file_count}, Segmented_2D{file_count},Segmented_2D_Ch2{file_count},Segmented_2D_Ch3{file_count},Segmented_2D_Ch4{file_count}]...
        =SegFlattenCrop3D_4Ch(XYZ{file_count},XYZ_Ch2{file_count},XYZ_Ch3{file_count},XYZ_Ch4{file_count},MetaData{file_count})
end
'2D segmentation complete'


toc
%%



%%
% tic
% 'Manual Selection and AutoLabelling'
% Selected_Segmented_2d_Raw=cell(size(Segmented_2D));
% Selected_Segmented_2d_AutoLabelled=cell(size(Segmented_2D));
% 
% for condition_count=1:size(Segmented_2D,1)
%     for file_count=1:size(Segmented_2D,2)
%         ['file ',num2str(file_count),' of ',num2str(size(Segmented_2D,2))]
%         [Selected_Segmented_2d_Raw{condition_count,file_count},Selected_Segmented_2d_AutoLabelled{condition_count,file_count},MetaData{condition_count,file_count}]= manual_obj_selection(Segmented_2D{condition_count,file_count}, Raw_cropped_2D{condition_count,file_count}, MetaData{condition_count,file_count});
%         
%     end
% end
% 
% clearvars condition_count file_count num_of_2d_fields
% 'Manual Selection and AutoLabelling complete'
% toc


tic
'Manual Selection and AutoLabelling'
Selected_Segmented_2d_Raw=cell(size(Segmented_2D));
Selected_Segmented_2d_AutoLabelled=cell(size(Segmented_2D));

Selected_Segmented_2d_Raw_Ch2=cell(size(Segmented_2D));
Selected_Segmented_2d_Raw_Ch3=cell(size(Segmented_2D));
Selected_Segmented_2d_Raw_Ch4=cell(size(Segmented_2D));

for condition_count=1:size(Segmented_2D,1)
    for file_count=1:size(Segmented_2D,2)
        ['file ',num2str(file_count),' of ',num2str(size(Segmented_2D,2))]
        [Selected_Segmented_2d_Raw{condition_count,file_count},Selected_Segmented_2d_Raw_Ch2{condition_count,file_count},Selected_Segmented_2d_Raw_Ch3{condition_count,file_count},Selected_Segmented_2d_Raw_Ch4{condition_count,file_count},Selected_Segmented_2d_AutoLabelled{condition_count,file_count},MetaData{condition_count,file_count}]...
            = manual_obj_selection_4Ch(...
            Segmented_2D{condition_count,file_count},Segmented_2D_Ch2{condition_count,file_count},Segmented_2D_Ch3{condition_count,file_count},Segmented_2D_Ch4{condition_count,file_count}, Raw_cropped_2D{condition_count,file_count}, MetaData{condition_count,file_count});
        
    end
end

clearvars condition_count file_count num_of_2d_fields
'Manual Selection and AutoLabelling complete'
toc

%%





%%
tic

MetaData=[];
MetaData.Voxel_Size_X = 0.2; % in µm
MetaData.Voxel_Size_Y =MetaData.Voxel_Size_X;
MetaData.Filename='xxx';

Selected_Segmented_2d_Raw={};
Selected_Segmented_2d_Raw{1, 1}=A;
Selected_Segmented_2d_Raw{1, 2}=B;
Selected_Segmented_2d_Raw{1, 3}=C;
Selected_Segmented_2d_Raw{1, 4}=D;
Selected_Segmented_2d_Raw{1, 5}=E;


Basic_Measurements=cell(size(Selected_Segmented_2d_Raw));
% Boundary_RadOfCurva_K=cell(size(XYZ));
Norm_Pix_Intensities=cell(size(Selected_Segmented_2d_Raw));
num_of_pix_in_largest_nuc=zeros(size(Selected_Segmented_2d_Raw));
num_of_nuc=zeros(size(Selected_Segmented_2d_Raw));
total_num_of_nuc=zeros(size(Selected_Segmented_2d_Raw,1),1);

'Basic Measurements'
for condition_count=1:size(Selected_Segmented_2d_Raw,1)
    parfor file_count=1:size(Selected_Segmented_2d_Raw,2)
        if size(Selected_Segmented_2d_Raw{condition_count,file_count},1)*size(Selected_Segmented_2d_Raw{condition_count,file_count},2)>0
            
            ['File ',num2str(file_count),' of ',num2str(size(Selected_Segmented_2d_Raw,2)),' from condition ',num2str(condition_count),' of ',num2str(size(Selected_Segmented_2d_Raw,1))]
            
            %             [Basic_Measurements{condition_count,file_count},Boundary_RadOfCurva_K{condition_count,file_count},Norm_Pix_Intensities{condition_count,file_count}]=basic_measurements(Selected_proper_2d_Raw{condition_count,file_count},MetaData{condition_count,file_count});
            [Basic_Measurements{condition_count,file_count},~,Norm_Pix_Intensities{condition_count,file_count}]=basic_measurements(Selected_Segmented_2d_Raw{condition_count,file_count},MetaData);
            
        end
    end
    
    %Reorganizing Boundary_RadOfCurva_K into Result_Boundary_RadOfCurva_K
    %adding nuc_label to the Result_Basic_Meansurements sheet
    %adding nuc_label to the Result_Norm_Pix_Intensities and reorganizing it.
    'Reorganizing output'
    
    for file_count1=1:size(Selected_Segmented_2d_Raw,2)
        number=size(Norm_Pix_Intensities{condition_count,file_count1},1);
        if number~=0
            number=number-2;
        end
        num_of_pix_in_largest_nuc(condition_count,file_count1)=number;
        num_of_nuc(condition_count,file_count1)=size(Norm_Pix_Intensities{condition_count,file_count1},2);
    end
    
    total_num_of_nuc(condition_count,1)=sum(num_of_nuc(condition_count,:));
    
    Result_Norm_Pix_Intensities=cell(max(num_of_pix_in_largest_nuc(condition_count,:)),total_num_of_nuc(condition_count,1));
    %     Result_Boundary_RadOfCurva_K=cell(size(Selected_proper_2d_Raw,1),1);
    lower_bound=1;
    for file_count=1:size(Selected_Segmented_2d_Raw,2)
        if size(Selected_Segmented_2d_Raw{condition_count,file_count},1)*size(Selected_Segmented_2d_Raw{condition_count,file_count},2)>0
            %             Result_Boundary_RadOfCurva_K{condition_count,1}(1:size(Boundary_RadOfCurva_K{condition_count,file_count},1),file_count)=Boundary_RadOfCurva_K{condition_count,file_count}(:,1);
            for nuc_count=1:size (Selected_Segmented_2d_Raw{condition_count,file_count},3)
                Norm_Pix_Intensities{condition_count,file_count}(2,nuc_count)=num2cell(1);%num2cell(mean(nonzeros(Selected_Segmented_2d_AutoLabelled{condition_count,file_count}(:,:,nuc_count))));
                Basic_Measurements{condition_count,file_count}{nuc_count,2}=1;%mean(nonzeros(Selected_Segmented_2d_AutoLabelled{condition_count,file_count}(:,:,nuc_count)));
            end
            
            upper_bound=lower_bound-1+size(Norm_Pix_Intensities{condition_count,file_count},2);
            Result_Norm_Pix_Intensities(1:size(Norm_Pix_Intensities{condition_count,file_count},1),lower_bound:upper_bound)=Norm_Pix_Intensities{condition_count,file_count};
            lower_bound=upper_bound+1;
        end
    end
    
    
end

header={'Filename' 'Nuc label' 'Centroid X(um)' 'Centroid Y(um)' 'Pro. Area(um^2)' 'Perimeter(um)' 'A.R.' 'Shape Factor' 'PDI' 'Centre Mismatch' 'I80_by_I20' 'nHigh_by_nLow' 'Mean of Normalized Int' 'Median of Normalized Int' 'S.D of Normalized Int' 'Mode of Normalized Int' 'Entropy' 'Relative concavity'...
    'Num. of times curvature changes polarity' 'Avg Curvature(um^-1)' 'SD in curvature(um^-1)'...
    'Max positive curvature(um^-1)' 'Total positive curvature(um^-1)' 'Avg positive curvature(um^-1)' 'Total Positive Curvature>0.2(um^-1)' 'Avg Positive Curvature>0.2(um^-1)' 'Length with positive curvature>0.2(um)' 'Fraction of perimeter with positive curvature>0.2'...
    'Num positive curvatures > 0.25' 'Avg prominance of positive curvatures > 0.25(um)' 'Avg width of positive curvatures > 0.25(um)'...
    'Max negative curvature(um^-1)' 'Total negative curvature(um^-1)' 'Avg negative curvature(um^-1)' 'Total negative Curvature<-0.005(um^-1)' 'Avg negative Curvature<-0.005(um^-1)' 'Length with negative Curvature<-0.005(um^-1)' 'Fraction of perimeter with negative Curvature<-0.005(um^-1)'...
    'Num negative curvatures<-0.05' 'Avg prominance of negative curvatures<-0.05(um)' 'Avg width of negative curvatures<-0.05(um)'...
    'Average curvature x Perimeter'  'Avg positive curvature x Perimeter' 'Avg negative curvature x Perimeter'};


Result_Basic_Measurements=vertcat(header,Basic_Measurements{:});

Result_Combined_Norm_Pix_Intensities=Result_Norm_Pix_Intensities(3:end,:);
Result_Combined_Norm_Pix_Intensities=vertcat(Result_Combined_Norm_Pix_Intensities{:});

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
TotalCh2=[];
TotalCh3=[];
TotalCh4=[];
Ch2DNAnormalized=[];
Ch3DNAnormalized=[];
Ch4DNAnormalized=[];

for condition_count=1:size(Selected_Segmented_2d_Raw,1)
    for file_count=1:size(Selected_Segmented_2d_Raw,2)
        if ~isempty(Selected_Segmented_2d_Raw{condition_count,file_count})
            CorrLengthScales_temp=[];
            MomentInv_temp=[];
            
            TotalCh1_temp=[];
            TotalCh2_temp=[];
            TotalCh3_temp=[];
            TotalCh4_temp=[];
            Ch2DNAnormalized_temp=[];
            Ch3DNAnormalized_temp=[];
            Ch4DNAnormalized_temp=[];
            
            
            for img_count=1:size(Selected_Segmented_2d_Raw{condition_count,file_count},3)
                CorrLengthScales_temp(img_count,:)=Correlation_Lengthscales(Selected_Segmented_2d_Raw{condition_count,file_count}(:,:,img_count),MetaData.Voxel_Size_X);
                MomentInv_temp(img_count,:)=MomentInvariants(Selected_Segmented_2d_Raw{condition_count,file_count}(:,:,img_count),MetaData.Voxel_Size_X);
                
%                 TotalCh1_temp(img_count,1)=sum(sum(Selected_Segmented_2d_Raw{condition_count,file_count}(:,:,img_count)));
%                 TotalCh2_temp(img_count,1)=sum(sum(Selected_Segmented_2d_Raw_Ch2{condition_count,file_count}(:,:,img_count)));
%                 TotalCh3_temp(img_count,1)=sum(sum(Selected_Segmented_2d_Raw_Ch3{condition_count,file_count}(:,:,img_count)));
%                 TotalCh4_temp(img_count,1)=sum(sum(Selected_Segmented_2d_Raw_Ch4{condition_count,file_count}(:,:,img_count)));
%                 Ch2DNAnormalized_temp(img_count,:)=TotalCh2_temp(img_count,1)/TotalCh1_temp(img_count,1);
%                 Ch3DNAnormalized_temp(img_count,:)=TotalCh3_temp(img_count,1)/TotalCh1_temp(img_count,1);
%                 Ch4DNAnormalized_temp(img_count,:)=TotalCh4_temp(img_count,1)/TotalCh1_temp(img_count,1);
                
            end
            CorrLengthScales=vertcat(CorrLengthScales,CorrLengthScales_temp);
            MomentInv=vertcat(MomentInv,MomentInv_temp);
            
%             TotalCh1=vertcat(TotalCh1,TotalCh1_temp);
%             TotalCh2=vertcat(TotalCh2,TotalCh2_temp);
%             TotalCh3=vertcat(TotalCh3,TotalCh3_temp);
%             TotalCh4=vertcat(TotalCh4,TotalCh4_temp);
%             Ch2DNAnormalized=vertcat(Ch2DNAnormalized,Ch2DNAnormalized_temp);
%             Ch3DNAnormalized=vertcat(Ch3DNAnormalized,Ch3DNAnormalized_temp);
%             Ch4DNAnormalized=vertcat(Ch4DNAnormalized,Ch4DNAnormalized_temp);
            
            
        end
    end
    
    Result_Basic_Measurements=horzcat(Result_Basic_Measurements,...
        vertcat(LengthScale_header,num2cell(CorrLengthScales)),...
        vertcat(MomentInvariants_header,num2cell(MomentInv)));
%     ,...
%         vertcat('TotalCh1',num2cell(TotalCh1)),...
%         vertcat('TotalCh2',num2cell(TotalCh2)),...
%         vertcat('TotalCh3',num2cell(TotalCh3)),...
%         vertcat('TotalCh4',num2cell(TotalCh4)),...
%         vertcat('Ch2DNAnormalized',num2cell(Ch2DNAnormalized)),...
%         vertcat('Ch3DNAnormalized',num2cell(Ch3DNAnormalized)),...
%         vertcat('Ch4DNAnormalized',num2cell(Ch4DNAnormalized)));
%     
end



clearvars condition_count file_count img_count *header *temp



'Basic Measurements complete'




toc
%%











