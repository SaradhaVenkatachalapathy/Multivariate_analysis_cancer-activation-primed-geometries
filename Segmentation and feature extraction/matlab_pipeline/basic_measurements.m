

%%% written by D.S.JOKHUN on 04/09/2018



function [Basic_Measurements,Boundary_RadOfCurva_K,Norm_Pix_Intensities]=basic_measurements(segmented_raw_2d,MtaData)

% segmented_raw_2d=Selected_Segmented_2d_Raw{3}(:,:,2);
% MtaData=MetaData{1};


Filename=cell(size(segmented_raw_2d,2),1);
Nuc_label=zeros(size(segmented_raw_2d,2),1);
Centroid=zeros(size(segmented_raw_2d,2),2);
Pro_area=zeros(size(segmented_raw_2d,2),1);
Perimeter=zeros(size(segmented_raw_2d,2),1);
AR=zeros(size(segmented_raw_2d,2),1);
Shape_factor=zeros(size(segmented_raw_2d,2),1);
PDI=zeros(size(segmented_raw_2d,2),1);
Centre_mismatch=zeros(size(segmented_raw_2d,2),1);
I80_by_I20=zeros(size(segmented_raw_2d,2),1);
nHigh_by_nLow=zeros(size(segmented_raw_2d,2),1);
Mean_Norm_Int=zeros(size(segmented_raw_2d,2),1);
Median_Norm_Int=zeros(size(segmented_raw_2d,2),1);
SD_Norm_Int=zeros(size(segmented_raw_2d,2),1);
Mode_Norm_Int=zeros(size(segmented_raw_2d,2),1);
Entropy=zeros(size(segmented_raw_2d,2),1);

Relative_concavity=zeros(size(segmented_raw_2d,2),1);
AvgCurva=zeros(size(segmented_raw_2d,2),1);
AvgCurva_X_Peri=zeros(size(segmented_raw_2d,2),1);
SDinCurva=zeros(size(segmented_raw_2d,2),1);
AvgPosCurva=zeros(size(segmented_raw_2d,2),1);
TotPosCurva=zeros(size(segmented_raw_2d,2),1);
TotPosCurva_MoreThn0pt2=zeros(size(segmented_raw_2d,2),1);
AvgPosCurva_MoreThn0pt2=zeros(size(segmented_raw_2d,2),1);
LengthWithPosCurva_MoreThn0pt2=zeros(size(segmented_raw_2d,2),1);
FracOfPeriWithPosCurva_MoreThn0pt2=zeros(size(segmented_raw_2d,2),1);
AvgPosCurva_X_Peri=zeros(size(segmented_raw_2d,2),1);
AvgNegCurva=zeros(size(segmented_raw_2d,2),1);
TotNegCurva=zeros(size(segmented_raw_2d,2),1);
TotNegCurva_LessThnNeg0pt005=zeros(size(segmented_raw_2d,2),1);
AvgNegCurva_LessThnNeg0pt005=zeros(size(segmented_raw_2d,2),1);
LengthWithNegCurva_LessThnNeg0pt005=zeros(size(segmented_raw_2d,2),1);
FracOfPeriWithNegCurva_LessThnNeg0pt005=zeros(size(segmented_raw_2d,2),1);
AvgNegCurva_X_Peri=zeros(size(segmented_raw_2d,2),1);
n_CurvaChangesSign=zeros(size(segmented_raw_2d,2),1);
MaxPosCurva=zeros(size(segmented_raw_2d,2),1);
nPeaks_WithPosCurva_MoreThn0pt25=zeros(size(segmented_raw_2d,2),1);
AvgProminance_of_PeaksWithPosCurva_MoreThn0pt25=zeros(size(segmented_raw_2d,2),1);
AvgWidth_of_PeaksWithPosCurva_MoreThn0pt25=zeros(size(segmented_raw_2d,2),1);
MaxNegCurva=zeros(size(segmented_raw_2d,2),1);
nPeaks_WithNegCurva_LessThnNeg0pt05=zeros(size(segmented_raw_2d,2),1);
AvgProminance_of_PeaksWithNegCurva_LessThnNeg0pt05=zeros(size(segmented_raw_2d,2),1);
AvgWidth_of_PeaksWithNegCurva_LessThnNeg0pt05=zeros(size(segmented_raw_2d,2),1);


Boundary_RadOfCurva_K=cell(size(segmented_raw_2d,2),1);
%

for nuc_count=1:size(segmented_raw_2d,2);
    try
        Filename{nuc_count,1}=MtaData.Filename;

        raw=segmented_raw_2d{nuc_count};
        bw=raw>0;
        stats=regionprops(bw,raw,'PixelValues','Area','MajorAxisLength','MinorAxisLength','Perimeter', 'Centroid','WeightedCentroid', 'ConvexArea');

        Norm_Pix_Intensities{1,nuc_count}=MtaData.Filename;
        Norm_Pix_Intensities(3:2+size(stats.PixelValues,1),nuc_count)=num2cell(mat2gray(stats.PixelValues));

        Mean_Norm_Int(nuc_count,1)=mean(mat2gray(stats.PixelValues));
        Median_Norm_Int(nuc_count,1)=median(mat2gray(stats.PixelValues));
        SD_Norm_Int(nuc_count,1)=std(mat2gray(stats.PixelValues));

        figure
        h_IntDistri=histfit(mat2gray(stats.PixelValues),255,'kernel');
        [pks_IntDistri,locs_IntDistri] = findpeaks(h_IntDistri(2, 1).YData,h_IntDistri(2, 1).XData);
        h_IntDistri(close);
        Mode_Norm_Int(nuc_count,1)=min(locs_IntDistri(pks_IntDistri==max(pks_IntDistri)));


        Centroid(nuc_count,1)=stats.Centroid(1)*MtaData.Voxel_Size_X;
        Centroid(nuc_count,2)=stats.Centroid(2)*MtaData.Voxel_Size_Y;
        Pro_area(nuc_count,1)=(stats.Area * (MtaData.Voxel_Size_X*MtaData.Voxel_Size_Y));
        Perimeter(nuc_count,1)=(stats.Perimeter * MtaData.Voxel_Size_X);
        AR(nuc_count,1)=(stats.MajorAxisLength/stats.MinorAxisLength);
        Shape_factor(nuc_count,1)=((stats.Perimeter^2)/(4*pi*stats.Area));


        %% PDI
        NumOfPix=size(raw);
        distance_frm_cen_squared=zeros(NumOfPix);

        for countX=1:NumOfPix(2)
            for countY=1:NumOfPix(1)
                distance_frm_cen_squared(countY,countX)=((sqrt(((countX-stats.Centroid(1))^2)...
                    +(((NumOfPix(1)-countY+1)-stats.Centroid(2))^2)))...
                    *MtaData.Voxel_Size_X).^2;
            end
        end
        actual_IntMoment_2=distance_frm_cen_squared.*mat2gray(raw);
        total_NormInt=sum(sum(mat2gray(raw)));
        uniform_IntMoment_2=distance_frm_cen_squared.*(double(bw)*(total_NormInt/stats.Area));
        PDI(nuc_count,1)=sum(sum(actual_IntMoment_2))/sum(sum(uniform_IntMoment_2));  %peripheral distributon index

        %% Centre mismatch
        Centre_mismatch(nuc_count,1)=(sqrt(((stats.WeightedCentroid(1)-stats.Centroid(1))^2)...
            +((stats.WeightedCentroid(2)-stats.Centroid(2))^2)))...
            *MtaData.Voxel_Size_X;


        %% Intensity Histogram Analysis
        I_exclude_percentiles=prctile(single(nonzeros(raw)),[0.1,99.9]);   %elimimating extreme values from the image (e.g saturated pixels etc.)
        aft_excl_extremes=uint16(raw).*uint16(raw>=I_exclude_percentiles(1)).*uint16(raw<=I_exclude_percentiles(2));

        I_percentiles=prctile(single(nonzeros(aft_excl_extremes)),[20,80]);
        I80_by_I20(nuc_count,1)=I_percentiles(2)/I_percentiles(1);

        normalize_aft_excl_extremes=mat2gray(nonzeros(aft_excl_extremes));
        nHigh_by_nLow(nuc_count,1)=sum(normalize_aft_excl_extremes>=0.8)/sum(normalize_aft_excl_extremes<=0.2);


        Entropy(nuc_count,1)=entropy(normalize_aft_excl_extremes);



        %% global curvature analysis
        Relative_concavity(nuc_count,1)=(stats.ConvexArea-stats.Area)/stats.ConvexArea;



        %% local curvature analysis

        boundary=bwboundaries(bw);
        boundary_xy=[];
        boundary_xy(:,1)=boundary{1,1}(1:end-1,2);   %The last point is a repeat of the 1st point.
        boundary_xy(:,2)=boundary{1,1}(1:end-1,1);
        smooth_span=2.5/MtaData.Voxel_Size_X; %smooth span of 2um
        cond = mod(smooth_span,2)<1;  % =1 if remainder is less than 1 (rounding down will bring it to an even num) and =0 if remainder is more than 1 (rounding down will bring it to an odd number).
        smooth_span = floor(smooth_span);
        smooth_span = smooth_span+cond;  % if cond was 1, 1 will be added to floor(span). if cond was 0, floor(span) is already odd and nothing is added.
        %padding the x-y series for proper smoothing and proper subsequent fitting
        boundary_xy((2*smooth_span)+1:(2*smooth_span)+size(boundary_xy,1),:)=boundary_xy(:,:); %sifting the series by some rows
        boundary_xy(1:(2*smooth_span),:)=boundary_xy(end-(2*smooth_span)+1:end,:); %padding the top
        boundary_xy(end+1:end+(2*smooth_span),:)=boundary_xy((2*smooth_span)+1:(2*smooth_span)+(2*smooth_span),:); %padding the bottom
        boundary_xy(:,1)=smooth(boundary_xy(:,1),smooth_span,'lowess');
        boundary_xy(:,2)=smooth(boundary_xy(:,2),smooth_span,'lowess');
        boundary_xy(:,3:4)=0;

        [MtaData.Filename,' nuc No.',num2str(nuc_count),' of ',num2str(size(segmented_raw_2d,3)),' -boundary analysis']
        fit_span=1/MtaData.Voxel_Size_X; %%tangent 1 will be taken over 1um stretch before point P and tangent 2 will be taken over a 1um stretch after point P.
        cond = mod(fit_span,2)>=1;  % =1 if remainder is more or equal to 1 (rounding down will bring it to an odd num) and =0 if remainder is less than 1 (rounding down will bring it to an even number).
        fit_span = floor(fit_span);
        fit_span = fit_span+cond;  % if cond was 1, 1 will be added to floor(span). if cond was 0, floor(span) is already even and nothing is added.

        RadOfCur=zeros(size(boundary_xy,1),1);
        working_points=(smooth_span)+1:size(boundary_xy,1)-(smooth_span);
        syms x y   %declaring system of equations
        for point_count=working_points;
    %         point_count
            line_1=fit(boundary_xy(point_count-fit_span:point_count,1),boundary_xy(point_count-fit_span:point_count,2),'poly1');
            coeffvalues_1=coeffvalues(line_1);

            max_tan_grad=10^10;
            min_tan_grad=10^-10;

            if abs(coeffvalues_1(1,1))>max_tan_grad
                norm_eqn_1 = y==boundary_xy(point_count-(fit_span/2),2);
                norm_grad_1=0;
            end
            if abs(coeffvalues_1(1,1))<min_tan_grad
                norm_eqn_1 = x==boundary_xy(point_count-(fit_span/2),1);
                norm_grad_1=Inf;
            end
            if abs(coeffvalues_1(1,1)) <=max_tan_grad
                if abs(coeffvalues_1(1,1)) >=min_tan_grad
                    norm_grad_1=-1/(coeffvalues_1(1,1));
                    c1=(boundary_xy(point_count-(fit_span/2),2)) - (norm_grad_1*(boundary_xy(point_count-(fit_span/2),1))); %c=y-mx
                    norm_eqn_1 = y==(norm_grad_1*x)+c1;
                end
            end

            line_2=fit(boundary_xy(point_count:point_count+fit_span,1),boundary_xy(point_count:point_count+fit_span,2),'poly1');
            coeffvalues_2=coeffvalues(line_2);
            if abs(coeffvalues_2(1,1))>max_tan_grad
                norm_eqn_2 = y==boundary_xy(point_count+(fit_span/2),2);
                norm_grad_2=0;
            end
            if abs(coeffvalues_2(1,1))<min_tan_grad
                norm_eqn_2 = x==boundary_xy(point_count+(fit_span/2),1);
                norm_grad_2=Inf;
            end
            if abs(coeffvalues_2(1,1)) <=max_tan_grad
                if abs(coeffvalues_2(1,1)) >=min_tan_grad
                    norm_grad_2=-1/(coeffvalues_2(1,1));
                    c2=(boundary_xy(point_count+(fit_span/2),2)) - (norm_grad_2*(boundary_xy(point_count+(fit_span/2),1))); %c=y-mx
                    norm_eqn_2 = y==(norm_grad_2*x)+c2;
                end
            end

            if round(norm_grad_1,10) == round(norm_grad_2,10)
                radius_of_curvature=Inf;
            else
                sol = solve([norm_eqn_1, norm_eqn_2], [x, y]);
                centre_of_curvature=double([sol.x,sol.y]);
                radius_of_curvature=double((sqrt(((boundary_xy(point_count,1)-sol.x)^2)+((boundary_xy(point_count,2)-sol.y)^2))))*MtaData.Voxel_Size_X;

                cen_of_curva_to_pt=[boundary_xy(point_count,1)-centre_of_curvature(1,1),boundary_xy(point_count,2)-centre_of_curvature(1,2)];   %vector to move from centre of curvature to point P on the boundary
                unit_vec=cen_of_curva_to_pt./norm(cen_of_curva_to_pt);  % unit vector to move from centre of curvature to point P on the boundary
                query_point=[boundary_xy(point_count,1)+(2*unit_vec(1,1)),boundary_xy(point_count,2)+(2*unit_vec(1,2))];
                if inpolygon(query_point(1,1),query_point(1,2),boundary_xy(smooth_span+1:length(boundary_xy)-smooth_span,1),boundary_xy(smooth_span+1:length(boundary_xy)-smooth_span,2))==1    %If point P+2*(unit vector from the centre of curvature to point P) lies inside the cell, this statement will be true
                    radius_of_curvature=-radius_of_curvature;
                end
            end


            RadOfCur(point_count,1)=radius_of_curvature;


        end


        boundary_xy(working_points,3)=RadOfCur(working_points,1);
        boundary_xy(:,4)=smooth(1./boundary_xy(:,3),fit_span,'rlowess');

        pre_Boundary_RadOfCurva_K=boundary_xy((2*smooth_span)+1:end-(2*smooth_span),1:4);
    %     Boundary_RadOfCurva_K_header={'x(pix)' 'y(pix)' 'Radius of curvature(um)' 'curvature,K=1/r(um^-1)'};
    %     Boundary_RadOfCurva_K{nuc_count,1}=vertcat(Boundary_RadOfCurva_K_header,num2cell(pre_Boundary_RadOfCurva_K));

        AvgCurva(nuc_count,1)=mean(pre_Boundary_RadOfCurva_K(:,4));
        AvgCurva_X_Peri(nuc_count,1)=mean(pre_Boundary_RadOfCurva_K(:,4))*Perimeter(nuc_count,1);
        SDinCurva(nuc_count,1)=std(pre_Boundary_RadOfCurva_K(:,4));
        AvgPosCurva(nuc_count,1)=(mean(nonzeros(pre_Boundary_RadOfCurva_K(:,4).*(pre_Boundary_RadOfCurva_K(:,4)>0))));
        TotPosCurva(nuc_count,1)=(sum(nonzeros(pre_Boundary_RadOfCurva_K(:,4).*(pre_Boundary_RadOfCurva_K(:,4)>0))));
        TotPosCurva_MoreThn0pt2(nuc_count,1)=(sum(nonzeros(pre_Boundary_RadOfCurva_K(:,4).*(pre_Boundary_RadOfCurva_K(:,4)>0.2))));
        AvgPosCurva_MoreThn0pt2(nuc_count,1)=(mean(nonzeros(pre_Boundary_RadOfCurva_K(:,4).*(pre_Boundary_RadOfCurva_K(:,4)>0.2))));
        LengthWithPosCurva_MoreThn0pt2(nuc_count,1)=(sum(pre_Boundary_RadOfCurva_K(:,4)>0.2))*MtaData.Voxel_Size_X;
        FracOfPeriWithPosCurva_MoreThn0pt2(nuc_count,1)=((sum(pre_Boundary_RadOfCurva_K(:,4)>0.2))*MtaData.Voxel_Size_X)/Perimeter(nuc_count,1);
        AvgPosCurva_X_Peri(nuc_count,1)=(mean(nonzeros(pre_Boundary_RadOfCurva_K(:,4).*(pre_Boundary_RadOfCurva_K(:,4)>0)))) *Perimeter(nuc_count,1);
        AvgNegCurva(nuc_count,1)=(mean(nonzeros(pre_Boundary_RadOfCurva_K(:,4).*(pre_Boundary_RadOfCurva_K(:,4)<0))));
        TotNegCurva(nuc_count,1)=(sum(nonzeros(pre_Boundary_RadOfCurva_K(:,4).*(pre_Boundary_RadOfCurva_K(:,4)<0))));
        TotNegCurva_LessThnNeg0pt005(nuc_count,1)=(sum(nonzeros(pre_Boundary_RadOfCurva_K(:,4).*(pre_Boundary_RadOfCurva_K(:,4)<-0.005))));
        AvgNegCurva_LessThnNeg0pt005(nuc_count,1)=(mean(nonzeros(pre_Boundary_RadOfCurva_K(:,4).*(pre_Boundary_RadOfCurva_K(:,4)<-0.005))));
        LengthWithNegCurva_LessThnNeg0pt005(nuc_count,1)=(sum(pre_Boundary_RadOfCurva_K(:,4)<-0.005))*MtaData.Voxel_Size_X;
        FracOfPeriWithNegCurva_LessThnNeg0pt005(nuc_count,1)=((sum(pre_Boundary_RadOfCurva_K(:,4)<-0.005))*MtaData.Voxel_Size_X)/Perimeter(nuc_count,1);
        AvgNegCurva_X_Peri(nuc_count,1)=(mean(nonzeros(pre_Boundary_RadOfCurva_K(:,4).*(pre_Boundary_RadOfCurva_K(:,4)<0)))) *Perimeter(nuc_count,1);

        SignTest=diff((nonzeros(pre_Boundary_RadOfCurva_K(:,4)))>0);
        n_CurvaChangesSign(nuc_count,1)=numel(nonzeros(SignTest));

        MaxPosCurva(nuc_count,1)= max(findpeaks(pre_Boundary_RadOfCurva_K(:,4)));
        [pks,~,w,p] =findpeaks(pre_Boundary_RadOfCurva_K(:,4),'MinPeakDistance',round(3/MtaData.Voxel_Size_X),'MinPeakHeight',0.25);   %If 2 peaks are reparated by less than 3um along the boundary, they are counted as 1
        nPeaks_WithPosCurva_MoreThn0pt25(nuc_count,1)= numel(pks); 
        AvgProminance_of_PeaksWithPosCurva_MoreThn0pt25(nuc_count,1)=mean(p);
        AvgWidth_of_PeaksWithPosCurva_MoreThn0pt25(nuc_count,1)=mean(w);
        MaxNegCurva(nuc_count,1)= -1*(max(findpeaks(-1*(pre_Boundary_RadOfCurva_K(:,4)))));
        [pks,~,w,p] =findpeaks(-1*(pre_Boundary_RadOfCurva_K(:,4)),'MinPeakDistance',round(3/MtaData.Voxel_Size_X),'MinPeakHeight',0.05);
        nPeaks_WithNegCurva_LessThnNeg0pt05(nuc_count,1)= numel(pks);
        AvgProminance_of_PeaksWithNegCurva_LessThnNeg0pt05(nuc_count,1)=mean(p);
        AvgWidth_of_PeaksWithNegCurva_LessThnNeg0pt05(nuc_count,1)=mean(w);



    % 
    %     figure ('Name',['Curvature,K of ',MtaData.Filename])
    %     colormap jet
    %     patch(pre_Boundary_RadOfCurva_K(:,1),size(bw,1)-pre_Boundary_RadOfCurva_K(:,2)+1,pre_Boundary_RadOfCurva_K(:,4),'EdgeColor','interp','FaceColor','none','LineWidth',5)
    %     figure
    %     imshow(raw,[])
    %     
    %     figure ('Name',['Curvature,K of ',MtaData.Filename])
    %     plot(pre_Boundary_RadOfCurva_K(:,4))


   
    
    
    

    catch 
     %Nothing to do
    end
 end
Basic_Measurements=horzcat(Filename,num2cell(Nuc_label),num2cell(Centroid),num2cell(Pro_area),num2cell(Perimeter),num2cell(AR),num2cell(Shape_factor),num2cell(PDI),num2cell(Centre_mismatch),num2cell(I80_by_I20),num2cell(nHigh_by_nLow),num2cell(Mean_Norm_Int),num2cell(Median_Norm_Int),num2cell(SD_Norm_Int),num2cell(Mode_Norm_Int),num2cell(Entropy),num2cell(Relative_concavity),...
    num2cell(n_CurvaChangesSign),num2cell(AvgCurva),num2cell(SDinCurva),...
    num2cell(MaxPosCurva),num2cell(TotPosCurva),num2cell(AvgPosCurva),num2cell(TotPosCurva_MoreThn0pt2),num2cell(AvgPosCurva_MoreThn0pt2),num2cell(LengthWithPosCurva_MoreThn0pt2),num2cell(FracOfPeriWithPosCurva_MoreThn0pt2),...
    num2cell(nPeaks_WithPosCurva_MoreThn0pt25),num2cell(AvgProminance_of_PeaksWithPosCurva_MoreThn0pt25),num2cell(AvgWidth_of_PeaksWithPosCurva_MoreThn0pt25),...
    num2cell(MaxNegCurva),num2cell(TotNegCurva),num2cell(AvgNegCurva),num2cell(TotNegCurva_LessThnNeg0pt005),num2cell(AvgNegCurva_LessThnNeg0pt005),num2cell(LengthWithNegCurva_LessThnNeg0pt005),num2cell(FracOfPeriWithNegCurva_LessThnNeg0pt005),...
    num2cell(nPeaks_WithNegCurva_LessThnNeg0pt05),num2cell(AvgProminance_of_PeaksWithNegCurva_LessThnNeg0pt05),num2cell(AvgWidth_of_PeaksWithNegCurva_LessThnNeg0pt05),...
    num2cell(AvgCurva_X_Peri),num2cell(AvgPosCurva_X_Peri),num2cell(AvgNegCurva_X_Peri));
    
end










