#Functions to read in the data
{
  open_and_compile <- function(dir, exp, samplen, trial,file_type){
    setwd(dir)
    file_data<-list.files()
    for (j in 1:length(file_data)){
      if (!exists("dataset")){
        if(file_type=="tsv"){
          dataset <- read.table(file_data[j], header=TRUE, sep="\t",stringsAsFactors=FALSE)
        }
        else if(file_type=="csv"){
          dataset <- read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
        }
      }
      else if (exists("dataset")){
        if(file_type=="tsv"){
          temp_dataset <- read.table(file_data[j], header=TRUE, sep="\t",stringsAsFactors=FALSE)
        }
        else if(file_type=="csv"){
          temp_dataset <- read.table(file_data[j], header=TRUE, sep=",",stringsAsFactors=FALSE)
        }
        if(nrow(temp_dataset)>0){
          dataset<-rbind(dataset, temp_dataset)
        }
        rm(temp_dataset)
      }
    }
    dataset$Exp<-exp
    dataset$sample<-samplen
    dataset$trial<-trial
    return(dataset)
  }
  open_and_summary_hubs <- function(dir, exp, samplen, trial, cols){
    setwd(dir)
    file_data<-list.files()
    
    summary_pol2<-as.data.frame(matrix(nrow=length(file_data), ncol=(length(cols)*6+2)))
    colnames(summary_pol2)[1]<-"Label"
    colnames(summary_pol2)[2]<-"Number_pol2_hubs"
    k=3;
    for(i in 1:length(cols)){
      colnames(summary_pol2)[k]<-paste(cols[i],"mean",sep="_")
      colnames(summary_pol2)[k+1]<-paste(cols[i],"median",sep="_")
      colnames(summary_pol2)[k+2]<-paste(cols[i],"max",sep="_")
      colnames(summary_pol2)[k+3]<-paste(cols[i],"min",sep="_")
      colnames(summary_pol2)[k+4]<-paste(cols[i],"Stdev",sep="_")
      colnames(summary_pol2)[k+5]<-paste(cols[i],"sum",sep="_")
      k=k+6;
    }
    
    for (j in 1:length(file_data)){
      if (!exists("dataset")){
        dataset <- read.table(file_data[j], header=TRUE, sep="\t",stringsAsFactors=FALSE)
        dataset$ImageName<-sub('\\_pol2.*', '', dataset$Label[1])
        dataset[is.nan(as.matrix(dataset))]<-0
        summary_pol2$Label[j]<-sub('\\_pol2.*', '', dataset$Label[1])
        summary_pol2$Number_pol2_hubs[j]<-nrow(dataset)
        s=2;
        for(k in 1:length(cols)){
          a=1+s;
          d<-dataset[,which(colnames(dataset)==cols[k])]
          d<-d[is.finite(d)]
          summary_pol2[j,a]<-mean(d,na.rm=TRUE)
          summary_pol2[j,a+1]<-median(d,na.omit=TRUE)
          summary_pol2[j,a+2]<-max(d,na.rm=TRUE)
          summary_pol2[j,a+3]<-min(d,na.rm=TRUE)
          summary_pol2[j,a+4]<-sd(d,na.rm=TRUE)
          summary_pol2[j,a+5]<-sum(d,na.rm=TRUE)
          s=a+5;
        }
        
        
      }
      else if (exists("dataset")){
        temp_dataset <- read.table(file_data[j], header=TRUE, sep="\t",stringsAsFactors=FALSE)
        temp_dataset$ImageName<-sub('\\_pol2.*', '', temp_dataset$Label[1])
        temp_dataset[is.nan(as.matrix(temp_dataset))]<-0
        
        summary_pol2$Label[j]<-sub('\\_pol2.*', '', temp_dataset$Label[1])
        summary_pol2$Number_pol2_hubs[j]<-nrow(temp_dataset)
        s=2;
        for(k in 1:length(cols)){
          a=1+s;
          d<-temp_dataset[,which(colnames(temp_dataset)==cols[k])]
          d<-d[is.finite(d)]
          
          summary_pol2[j,a]<-mean(d,na.rm =TRUE)
          summary_pol2[j,a+1]<-median(d,na.omit=TRUE)
          summary_pol2[a+2]<-max(d,na.rm=TRUE)
          summary_pol2[j,a+3]<-min(d,na.rm=TRUE)
          summary_pol2[j,a+4]<-sd(d,na.rm=TRUE)
          summary_pol2[j,a+5]<-sum(d,na.rm=TRUE)
          s=a+5;
        }
        
        
        if(nrow(temp_dataset)>0){
          dataset<-rbind(dataset, temp_dataset)
        }
        rm(temp_dataset)
      }
    }
    dataset$Exp<-exp
    dataset$sample<-samplen
    dataset$trial<-trial
    summary_pol2$Exp<-exp
    summary_pol2$sample<-samplen
    summary_pol2$trial<-trial
    return(summary_pol2)
    
  }
  open_and_compile_hubs_nuclear_positioning_summary <- function(dir, exp, samplen, trial,cols){
    setwd(dir)
    file_data<-list.files()
    
    summary_pol2<-as.data.frame(matrix(nrow=length(file_data), ncol=(length(cols)*5+2)))
    colnames(summary_pol2)[1]<-"Label"
    colnames(summary_pol2)[2]<-"Number_pol2_hubs"
    k=3;
    for(i in 1:length(cols)){
      colnames(summary_pol2)[k]<-paste(cols[i],"mean_nuc",sep="_")
      colnames(summary_pol2)[k+1]<-paste(cols[i],"median_nuc",sep="_")
      colnames(summary_pol2)[k+2]<-paste(cols[i],"max_nuc",sep="_")
      colnames(summary_pol2)[k+3]<-paste(cols[i],"min_nuc",sep="_")
      colnames(summary_pol2)[k+4]<-paste(cols[i],"Stdev_nuc",sep="_")
      k=k+5;
    }
    
    for (j in 1:length(file_data)){
      if (!exists("dataset")){
        dataset <- read.table(file_data[j], header=TRUE, sep="\t",stringsAsFactors=FALSE)
        dataset$ImageName<-substring(sub('\\Distances.*', '', file_data[j]),3)
        dataset<-subset(dataset, dataset$Label1=="nucleus")
        summary_pol2$Label[j]<-substring(sub('\\Distances.*', '', file_data[j]),3)
        summary_pol2$Number_pol2_hubs[j]<-nrow(dataset)
        s=2;
        for(k in 1:length(cols)){
          a=1+s;
          d<-dataset[,which(colnames(dataset)==cols[k])]
          summary_pol2[j,a]<-mean(d,na.rm=TRUE)
          summary_pol2[j,a+1]<-median(d,na.rm=TRUE)
          summary_pol2[j,a+2]<-max(d,na.rm=TRUE)
          summary_pol2[j,a+3]<-min(d,na.rm=TRUE)
          summary_pol2[j,a+4]<-sd(d,na.rm=TRUE)
          s=a+4;
        }
      }
      else if (exists("dataset")){
        temp_dataset <- read.table(file_data[j], header=TRUE, sep="\t",stringsAsFactors=FALSE)
        temp_dataset$ImageName<-substring(sub('\\Distances.*', '', file_data[j]),3)
        temp_dataset<-subset(temp_dataset, temp_dataset$Label1=="nucleus")
        summary_pol2$Label[j]<-substring(sub('\\Distances.*', '', file_data[j]),3)
        summary_pol2$Number_pol2_hubs[j]<-nrow(temp_dataset)
        s=2;
        for(k in 1:length(cols)){
          a=1+s;
          d<-temp_dataset[,which(colnames(temp_dataset)==cols[k])]
          summary_pol2[j,a]<-mean(d,na.rm =TRUE)
          summary_pol2[j,a+1]<-median(d,na.rm=TRUE)
          summary_pol2[j,a+2]<-max(d,na.rm=TRUE)
          summary_pol2[j,a+3]<-min(d,na.rm=TRUE)
          summary_pol2[j,a+4]<-sd(d,na.rm=TRUE)
          s=a+4;
        }
        
        rm(temp_dataset)
      }
    }
    summary_pol2$Exp<-exp
    summary_pol2$sample<-samplen
    summary_pol2$trial<-trial
    return(summary_pol2)
    
  }
  open_and_compile_hubs_relative_distances_summary <- function(dir, exp, samplen, trial,cols){
    setwd(dir)
    file_data<-list.files()
    
    summary_pol2<-as.data.frame(matrix(nrow=length(file_data), ncol=(length(cols)*5+2)))
    colnames(summary_pol2)[1]<-"Label"
    colnames(summary_pol2)[2]<-"Number_pol2_hubs"
    k=3;
    for(i in 1:length(cols)){
      colnames(summary_pol2)[k]<-paste(cols[i],"mean_rel",sep="_")
      colnames(summary_pol2)[k+1]<-paste(cols[i],"median_rel",sep="_")
      colnames(summary_pol2)[k+2]<-paste(cols[i],"max_rel",sep="_")
      colnames(summary_pol2)[k+3]<-paste(cols[i],"min_rel",sep="_")
      colnames(summary_pol2)[k+4]<-paste(cols[i],"Stdev_rel",sep="_")
      k=k+5;
    }
    
    
    for (j in 1:length(file_data)){
      if (!exists("dataset")){
        dataset <- read.table(file_data[j], header=TRUE, sep="\t",stringsAsFactors=FALSE)
        dataset$ImageName<-substring(sub('\\Distances.*', '', file_data[j]),3)
        dataset<-subset(dataset, dataset$Label1!="nucleus")
        dataset<-subset(dataset, dataset$Label2!="nucleus")
        summary_pol2$Label[j]<-substring(sub('\\Distances.*', '', file_data[j]),3)
        summary_pol2$Number_pol2_hubs[j]<-nrow(dataset)
        s=2;
        for(k in 1:length(cols)){
          a=1+s;
          d<-dataset[,which(colnames(dataset)==cols[k])]
          summary_pol2[j,a]<-mean(d,na.rm=TRUE)
          summary_pol2[j,a+1]<-median(d,na.rm=TRUE)
          summary_pol2[j,a+2]<-max(d,na.rm=TRUE)
          summary_pol2[j,a+3]<-min(d,na.rm=TRUE)
          summary_pol2[j,a+4]<-sd(d,na.rm=TRUE)
          s=a+4;
        }
      }
      else if (exists("dataset")){
        temp_dataset <- read.table(file_data[j], header=TRUE, sep="\t",stringsAsFactors=FALSE)
        temp_dataset$ImageName<-substring(sub('\\Distances.*', '', file_data[j]),3)
        temp_dataset<-subset(temp_dataset, temp_dataset$Label1!="nucleus")
        temp_dataset<-subset(temp_dataset, temp_dataset$Label2!="nucleus")
        summary_pol2$Label[j]<-substring(sub('\\Distances.*', '', file_data[j]),3)
        summary_pol2$Number_pol2_hubs[j]<-nrow(dataset)
        s=2;
        for(k in 1:length(cols)){
          a=1+s;
          d<-temp_dataset[,which(colnames(temp_dataset)==cols[k])]
          summary_pol2[j,a]<-mean(d,na.rm=TRUE)
          summary_pol2[j,a+1]<-median(d,na.rm=TRUE)
          summary_pol2[j,a+2]<-max(d,na.rm=TRUE)
          summary_pol2[j,a+3]<-min(d,na.rm=TRUE)
          summary_pol2[j,a+4]<-sd(d,na.rm=TRUE)
          s=a+4;
        }
        
        rm(temp_dataset)
      }
    }
    summary_pol2$Exp<-exp
    summary_pol2$sample<-samplen
    summary_pol2$trial<-trial
    return(summary_pol2)
  }
  
}

path_to_sample="E:/co_culture/20171130_NIH3t3_mcf7_co_culture_phalloidin_647_pol2_488_MKL_568_DAPI/co_culture/"
experiment_id="20171130_NIH3t3_mcf7_co_culture_phalloidin_647_pol2_488_MKL_568_DAPI"
sample_id="Co_culture"
trial="T1"
filetype="tsv"
#Read in the measurement files
{
  geometrical_data<-open_and_compile(paste(path_to_sample,"3D geometrical data/",sep=""),experiment_id,sample_id,trial,filetype)
  geometrical_simple<-open_and_compile(paste(path_to_sample,"3D geometerical_simple/",sep=""),experiment_id,sample_id,trial,filetype)
  geometrical_shape<-open_and_compile(paste(path_to_sample,"3D shape measure/",sep=""),experiment_id,sample_id,trial,filetype)
  geometrical_ellipsoid<-open_and_compile(paste(path_to_sample,"3D ellipsoid/results/",sep=""),experiment_id,sample_id,trial,filetype)
  
  intensity_nuc<-open_and_compile(paste(path_to_sample,"3D int_data/DNA/",sep=""),experiment_id,sample_id,trial,filetype)
  intensity_mkl<-open_and_compile(paste(path_to_sample,"3D int_data/MKL/",sep=""),experiment_id,sample_id,trial,filetype)
  geo_int_2D<-read.csv(paste(path_to_sample,"/2D measures/2D_indivisual_nuclei_ch1.csv",sep=""), header=T, stringsAsFactors = F)
  geo_int_2D$Exp<-experiment_id
  geo_int_2D$sample<-sample_id
  geo_int_2D$trial<-trial
  compaction_measure<-read.csv(paste(path_to_sample,"/2D measures/2D_indivisual_nuclei_ch1.csv",sep=""), header=T, stringsAsFactors = F)
  compaction_measure$Exp<-experiment_id
  compaction_measure$sample<-sample_id
  compaction_measure$trial<-trial
  matlab_pipeline<-read.csv(paste(path_to_sample,"/kamal_pipeline.csv",sep=""), header=F, stringsAsFactors = F)
  colnames(matlab_pipeline)<-matlab_pipeline[2,]
  matlab_pipeline<-matlab_pipeline[-c(1,2),]
  matlab_pipeline$Exp<-experiment_id
  matlab_pipeline$sample<-sample_id
  matlab_pipeline$trial<-trial
  
  list_names<-c("Vol..unit." )
  geometrical_data_HC_summary<-open_and_summary_hubs(paste(path_to_sample,"3D geometrical data HC",sep=""),experiment_id,sample_id,trial,list_names)
  list_names<-c("cen.cen","cen.bor","bor.bor","excen","bor.rad","periph")
  nuclear_positioning_HC_summary<-open_and_compile_hubs_nuclear_positioning_summary(paste(path_to_sample,"3D distance data HC",sep=""),experiment_id,sample_id,trial,list_names)
  list_names<-c("cen.cen","bor.bor")
  relative_positioning_HC_summary<-open_and_compile_hubs_relative_distances_summary(paste(path_to_sample,"3D distance data HC",sep=""),experiment_id,sample_id,trial,list_names)
  
  #reading the qality control labels
  labels<-read.csv(paste(path_to_sample,"/data/Log.txt",sep=""), header=F, stringsAsFactors = F)
  labels$Exp<-experiment_id
  labels$sample<-sample_id
  labels$trial<-trial
  
}
#make column names unique
{
  colnames(geo_int_2D)<-c("X.1","Label","proj_nuc_Area","proj_nuc_Mean","proj_nuc_StdDev","proj_nuc_Mode","proj_nuc_Min","proj_nuc_Max","proj_nuc_X","proj_nuc_Y",
                          "proj_nuc_XM","proj_nuc_YM","proj_nuc_Perim.","proj_nuc_BX","proj_nuc_BY","proj_nuc_Width","proj_nuc_Height","proj_nuc_Major","proj_nuc_Minor",
                          "proj_nuc_Angle","proj_nuc_Circ.","proj_nuc_Feret","proj_nuc_IntDen","proj_nuc_Median","proj_nuc_Skew","proj_nuc_Kurt","proj_nuc_RawIntDen",
                          "proj_nuc_FeretX","proj_nuc_FeretY","proj_nuc_FeretAngle","proj_nuc_MinFeret","proj_nuc_AR","proj_nuc_Round","proj_nuc_Solidity","min Thresh","Max thresh","Exp","sample","trial"      
  )
  
  colnames(intensity_mkl)<-c("Nb","Obj","Label","mkl_AtCenter","mkl_CMx..pix.","mkl_CMy..pix.","mkl_CMy..pix..1","mkl_CMx..unit.","mkl_CMy..unit.","mkl_CMy..unit..1",
                             "mkl_IntDen","mkl_Min","mkl_Max","mkl_Mean","mkl_Sigma" ,"mkl_Mode","mkl_Mode.NonZero","X","Exp","sample","trial")
  colnames(intensity_nuc)<-c("Nb","Obj","Label","DNA_AtCenter","DNA_CMx..pix.","DNA_CMy..pix.","DNA_CMy..pix..1","DNA_CMx..unit.","DNA_CMy..unit.",
                             "DNA_CMy..unit..1","DNA_IntDen","DNA_Min","DNA_Max","DNA_Mean","DNA_Sigma" ,"DNA_Mode","DNA_Mode.NonZero","X","Exp","sample","trial")
}

#make labels uniform
{
  geo_int_2D$Label<-substring(geo_int_2D$Label,5, nchar(geo_int_2D$Label)-4)
  compaction_measure$Label<-substring(compaction_measure$Label1,1, nchar(compaction_measure$Label1)-4)
  geometrical_data_HC_summary$Label<- substring( sub("\\H.*", "", geometrical_data_HC_summary$Label), 1, nchar(sub("\\H.*", "", geometrical_data_HC_summary$Label))-1)
  labels$Label<-substring(labels$Label,1, nchar(labels$Label)-4)
  matlab_pipeline$Label<-substring(matlab_pipeline$Filename,5, nchar(matlab_pipeline$Filename)-4)
  
}
area_2d<-c("Label","proj_nuc_Area","proj_nuc_Perim.","proj_nuc_Width","proj_nuc_Height","proj_nuc_Major","proj_nuc_Minor","proj_nuc_Angle","proj_nuc_Circ.","proj_nuc_Feret","proj_nuc_MinFeret","proj_nuc_AR","proj_nuc_Round")
int_2d<-c("Label","proj_nuc_Mean","proj_nuc_StdDev","proj_nuc_Mode","proj_nuc_Min","proj_nuc_Max","proj_nuc_IntDen","proj_nuc_Median","proj_nuc_Skew","proj_nuc_Kurt","proj_nuc_RawIntDen")


# combine parameters for for informative datasets
{
  #combine 3d nuc_morph_data
  combined_nuc_morph_3d<-cbind(geometrical_data,geometrical_ellipsoid)
  combined_nuc_morph_3d<-cbind(combined_nuc_morph_3d,geometrical_shape[-2])
  combined_nuc_morph_3d<-cbind(combined_nuc_morph_3d,geometrical_simple[-2])
  
  HC_data<-merge(geometrical_data_HC_summary,nuclear_positioning_HC_summary[,-2], by="Label")
  HC_data<-merge(HC_data,relative_positioning_HC_summary[,-2], by="Label")
  colnames(HC_data)[-c(1,56,57)] <- paste("HC", colnames(HC_data)[-c(1,56,57)], sep = "_")
  
  #combined all
  all<-merge(combined_nuc_morph_3d,geo_int_2D, by="Label")
  all<-merge(all,compaction_measure, by="Label")
  all<-merge(all,intensity_nuc, by="Label")
  all<-merge(all,matlab_pipeline[,-c(1,2)], by="Label")
  all<-merge(all,HC_data, by="Label")
  colnames(all)[341:342]<-c("sample","trial")
  
  #remove imporper nuclei
  all<-subset(all,all$Label %in% proper)
  
  
}



