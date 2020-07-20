// assumes TF is in channel 2 and pol2 is in channel 3 

dirsa=getArgument();

run("3D Manager Options", "volume surface compactness fit_ellipse 3d_moments integrated_density mean_grey_value std_dev_grey_value mode_grey_value feret minimum_grey_value maximum_grey_value centroid_(pix) centroid_(unit) distance_to_surface centre_of_mass_(pix) centre_of_mass_(unit) bounding_box radial_distance surface_contact closest use distance_between_centers=10 distance_max_contact=1.80");
run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value centroid mean_distance_to_surface std_dev_distance_to_surface median_distance_to_surface centre_of_mass bounding_box dots_size=5 font_size=20 redirect_to=none");
//dirsa=getDirectory("Please choose the source directory");
setBatchMode(true);
dirsa1=dirsa+"rawimages"+File.separator;
dir2 = dirsa + "HC and pol2 hubs" + File.separator;
newDir12 = dir2 + "HC_binary" + File.separator;

newDir42 = dirsa + "indivisual_nuclei_ch1" + File.separator;


Dir_res_hc=dirsa + "3D objects HC"+ File.separator;   
Dir_res_geo_hc=dirsa + "3D geometrical data HC"+ File.separator;  
Dir_res_int=dirsa + "3D HC int data"+ File.separator;  
Dir_res_ch1_int=Dir_res_int + "nucleus" +File.separator;  


File.makeDirectory(Dir_res_hc);
File.makeDirectory(Dir_res_geo_hc); 
File.makeDirectory(Dir_res_int); 
File.makeDirectory(Dir_res_ch1_int); 


filenames=getFileList(newDir12);
for(f=0;f<filenames.length;f++){
	open(newDir12+filenames[f]);
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tif"); 
	baseName=substring(obj_image, 0, baseNameEnd); 
	
	getVoxelSize(width, height, depth, unit);
	a=1000/(width*height*depth);// maximum volume is 100 cu.microns
	a_1=0/(width*height*depth);// minimum volume is 1 cu.microns
	
	run("3D Objects Counter", "threshold=128 min.=a_1 max.=a objects statistics"); //Size filter
	saveAs("Tiff", Dir_res_hc+baseName+".tiff"); 
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tiff"); 
	baseName=substring(obj_image, 0, baseNameEnd); 
	
	run("3D Manager");
	Ext.Manager3D_AddImage();
	Ext.Manager3D_Count(nb_obj); 
	for(p=0;p<nb_obj;p++){
		Ext.Manager3D_Bounding3D(p,x0,x1,y0,y1,z0,z1);
		Ext.Manager3D_MonoSelect();
		Ext.Manager3D_Select(p);
		name=baseName+"_HC_"+x0+"_"+y0+"_"+z0+"_"+z1;
		Ext.Manager3D_Rename(name);
	}
	Ext.Manager3D_DeselectAll();

	
	Ext.Manager3D_DeselectAll();
	Ext.Manager3D_SelectAll();
	Ext.Manager3D_Measure();
    Ext.Manager3D_SaveMeasure(Dir_res_geo_hc+baseName+"_geometric.tsv");
    Ext.Manager3D_CloseResult("M");
	run("Close All");

	path=newDir42+filenames[f];
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_step=1");
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveQuantif(Dir_res_ch1_int+baseName+".tsv");
	Ext.Manager3D_CloseResult("Q");
	run("Close All");
	Ext.Manager3D_Close();

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
