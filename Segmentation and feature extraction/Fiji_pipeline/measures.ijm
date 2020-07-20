
args = split(getArgument()," "); 
dirsa=args[0];
nucleus = args[1];
TF = args[2];
tfname= args[3];

//dirsa=getDirectory("Please choose the source directory");
setBatchMode(true);
dirsa1=dirsa+"rawimages"+File.separator;
dirb= dirsa + "3d objects"+ File.separator;
run("3D Manager Options", "volume surface compactness fit_ellipse 3d_moments integrated_density mean_grey_value std_dev_grey_value mode_grey_value feret minimum_grey_value maximum_grey_value centroid_(pix) centroid_(unit) distance_to_surface centre_of_mass_(pix) centre_of_mass_(unit) bounding_box radial_distance surface_contact closest use distance_between_centers=10 distance_max_contact=1.80");
Dir_res_ellp=dirsa + "3D ellipsoid"+ File.separator;
Dir_res_ellp1=Dir_res_ellp + "ellipsoid images"+ File.separator;  
Dir_res_ellp2=Dir_res_ellp + "results"+ File.separator; 
Dir_res_ellp3=Dir_res_ellp + "log"+ File.separator;   
Dir_res_geo=dirsa + "3D geometerical_simple"+ File.separator;
Dir_res_shape=dirsa + "3D shape measure"+ File.separator;
Dir_res_geo_m=dirsa + "3D geometrical data"+ File.separator;  
Dir_res_int=dirsa + "3D int_data"+ File.separator;  
Dir_res_ch1_int=Dir_res_int + "DNA" + File.separator;  
Dir_res_ch2_int=Dir_res_int + tfname +File.separator;  

File.makeDirectory(Dir_res_ellp); 
File.makeDirectory(Dir_res_ellp1); 
File.makeDirectory(Dir_res_ellp2); 
File.makeDirectory(Dir_res_ellp3); 
File.makeDirectory(Dir_res_geo); 
File.makeDirectory(Dir_res_shape); 
File.makeDirectory(Dir_res_geo_m); 
File.makeDirectory(Dir_res_int); 
File.makeDirectory(Dir_res_ch1_int); 
File.makeDirectory(Dir_res_ch2_int); 


filenames=getFileList(dirb);
filenamesa=getFileList(dirsa1);

for(f=0;f<filenames.length;f++){
	open(dirb+filenames[f]);
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
		name=baseName+"_nucleus_"+x0+"_"+y0+"_"+z0+"_"+z1;
		Ext.Manager3D_Rename(name);
	}
	Ext.Manager3D_DeselectAll();

	
	Ext.Manager3D_DeselectAll();
	Ext.Manager3D_SelectAll();
	Ext.Manager3D_Measure();
    Ext.Manager3D_SaveMeasure(Dir_res_geo_m+baseName+"_geometric.tsv");
    Ext.Manager3D_CloseResult("M");
	run("Close All");

	path=dirsa1+filenamesa[f];
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=nucleus c_end=nucleus c_step=1");
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveQuantif(Dir_res_ch1_int+baseName+".tsv");
	Ext.Manager3D_CloseResult("Q");
	run("Close All");

	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=TF c_end=TF c_step=1");
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveQuantif(Dir_res_ch2_int+baseName+".tsv");
	Ext.Manager3D_CloseResult("Q");
	run("Close All");
	
	Ext.Manager3D_Close();
}
//3D geometric and shape measures

for(f=0;f<filenames.length;f++){
	open(dirb+filenames[f]);
	obj_image=getTitle();
	baseNameEnd=indexOf(obj_image, ".tiff"); 
	baseName=substring(obj_image, 0, baseNameEnd); 

	run("3D Ellipsoid Fitting");
	selectWindow("Ellipsoids");
	name=baseName+"ellipsoid";
	saveAs("Tiff", Dir_res_ellp1+baseName+"ell"); 	
	selectWindow("Log");  //select Log-window 
	saveAs("Text", Dir_res_ellp3+baseName+"Log.txt"); 	
	selectWindow("Results");  //select Log-window 
	saveAs("Results", Dir_res_ellp2+baseName+"Results.csv");
	run("Clear Results");
	
	selectWindow(obj_image);   
	run("3D Geometrical Measure");
	selectWindow("Results");  //select Log-window 
	saveAs("Results", Dir_res_geo+baseName+".csv");
	run("Clear Results");
	
	selectWindow(obj_image);  
	run("3D Shape Measure");
	selectWindow("Results");  //select Log-window 
	saveAs("Results", Dir_res_shape+baseName+".csv");
	run("Clear Results");

	print("\\Clear");
	run("Close All");
}




