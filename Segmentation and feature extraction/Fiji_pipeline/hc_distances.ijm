// assumes TF is in channel 2 and pol2 is in channel 3 
args = split(getArgument()," "); 
dirsa=args[0];
index= args[1];
run("3D Manager Options", "volume surface compactness fit_ellipse 3d_moments integrated_density mean_grey_value std_dev_grey_value mode_grey_value feret minimum_grey_value maximum_grey_value centroid_(pix) centroid_(unit) distance_to_surface centre_of_mass_(pix) centre_of_mass_(unit) bounding_box radial_distance surface_contact closest use distance_between_centers=10 distance_max_contact=1.80");
run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value centroid mean_distance_to_surface std_dev_distance_to_surface median_distance_to_surface centre_of_mass bounding_box dots_size=5 font_size=20 redirect_to=none");

//dirsa=getDirectory("Please choose the source directory");
setBatchMode(true);
dirsa1=dirsa+"rawimages"+File.separator;
dir2 = dirsa + "HC and pol2 hubs" + File.separator;
newDir12 = dir2 + "HC_binary" + File.separator;
newDir32 = dirsa + "indivisual_nuclei_ch2" + File.separator;
newDir42 = dirsa + "indivisual_nuclei_ch1" + File.separator;
newDir52 = dirsa + "indivisual_nuclei_ch3" + File.separator;

Dir_res_hc=dirsa + "3D objects HC"+ File.separator;   
Dir_res_geo_hc=dirsa + "3D geometrical data HC"+ File.separator;  
Dir_res_int=dirsa + "3D HC int data"+ File.separator;  
Dir_res_ch1_int=Dir_res_int + "nucleus" +File.separator;  

run("3D Manager Options", "volume surface compactness fit_ellipse 3d_moments integrated_density mean_grey_value std_dev_grey_value mode_grey_value feret minimum_grey_value maximum_grey_value centroid_(pix) centroid_(unit) distance_to_surface centre_of_mass_(pix) centre_of_mass_(unit) bounding_box radial_distance surface_contact closest use distance_between_centers=10 distance_max_contact=1.80");

dir_nuc= dirsa + "indivisual_nuclei_ch1"+ File.separator;
Dir_res_hc=dirsa + "3D objects HC"+ File.separator;   
Dir_res_dist_hc=dirsa + "3D distance data hc"+ File.separator;  
File.makeDirectory(Dir_res_dist_hc); 
filenames_hc=getFileList(Dir_res_hc);
for(f=index;f<filenames_hc.length;f++){

	name=substring(filenames_hc[f], 0, lengthOf(filenames_hc[f])-1); 
	open(dir_nuc+name);
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd); 

	print(f);
	run("3D Objects Counter", "threshold=1 objects"); 
	//run("3D Convex Hull");//convex hull
	run("3D Manager");
	Ext.Manager3D_AddImage();
	Ext.Manager3D_Count(nb_obj);
	Ext.Manager3D_Select(0);
	Ext.Manager3D_Rename("nucleus");
	Ext.Manager3D_DeselectAll();

	open(Dir_res_hc+filenames_hc[f]);
	Ext.Manager3D_AddImage();
	Ext.Manager3D_Count(nb_obj2); 
	if(nb_obj2<100){
		for(p=nb_obj;p<nb_obj2;p++){
		Ext.Manager3D_Bounding3D(p,x0,x1,y0,y1,z0,z1);
		Ext.Manager3D_MonoSelect();
		Ext.Manager3D_Select(p);
		name=baseName+"_HC_"+x0+"_"+y0+"_"+z0+"_"+z1;
		Ext.Manager3D_Rename(name);
	}
	Ext.Manager3D_DeselectAll();

	
	Ext.Manager3D_DeselectAll();
	Ext.Manager3D_SelectAll();
    Ext.Manager3D_Distance();
    Ext.Manager3D_SaveResult("D",Dir_res_dist_hc+baseName+"Distances3D.tsv");
    Ext.Manager3D_CloseResult("D");
	}
	
	run("Close All");
	Ext.Manager3D_Close();
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
	run("Collect Garbage");
	
}













