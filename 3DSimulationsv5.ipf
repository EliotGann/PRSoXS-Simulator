#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=1		// Use modern global access method.
#include "opticalconstantsDB"
//#include "fftstuff2"
//#include "MCAlignment"
structure ThreeDSystem
// parameter inputs
	variable rot // boolean saying if we are rotating the system or not (90 degrees is always included) 
	variable num // resolution of simulation (in each of the long axes of the film)
	variable thickness // voxels of thickness along the beam direction
	variable voxelsize // the size of each pixel in the simulation (in nanometers)
	variable size // the dominant size of features (ie the radius of the sphere, if we're simulating spheres)
	variable materialnum // number of materials
	string paramstring //specific parameters for creating density matricies according to the model
	string modelname //store the model name somewhere
	string path // the path to put important things for this calculation
	string name // name of simulation
	variable movie // are we doing a movie

// morphology outputs, alignment inputs

	//technically not needed by simulations, but nice to have around for visualizations
	wave density1 // the density of material 1 output from the morphology module, this will be changed into
	wave density2 // optional density of material 2 (implies 3 materials)
	wave density3 // implies 4 materials
	wave density4 // implies 5 materials
	wave density5 // cannot be created by morphology, but only as the remnant of the other materials (to make total density everywhere 1)

// Inputs for alignment calculations
	string materialalignmentstring // list of materials alignment (for creating alignment)
	
	
// Morphology or alignment Outputs, Scattering Inputs

	wave m1 //4 dimensional density of material 1 in (x direction, y direction, z direction, unaligned)  sum = relative concentration sum(between -1 - 1 for x, y and z)^2 + unaligned = total material density
	wave m2 //4 dimensional density of material 2 in (x direction, y direction, z direction, unaligned)  sum = relative concentration sum(between -1 - 1 for x, y and z)^2 + unaligned = total material density
	wave m3 //4 dimensional density of material 3 in (x direction, y direction, z direction, unaligned)  sum = relative concentration sum(between -1 - 1 for x, y and z)^2 + unaligned = total material density
	wave m4 //4 dimensional density of material 4 in (x direction, y direction, z direction, unaligned)  sum = relative concentration sum(between -1 - 1 for x, y and z)^2 + unaligned = total material density
	wave m5 //4 dimensional density of material 5 in (x direction, y direction, z direction, unaligned)  sum = relative concentration sum(between -1 - 1 for x, y and z)^2 + unaligned = total material density

// parameter inputs for scattering
	string materials // list of materials
	string efield // two components which add to 1 which indicate the unit direction of the efield in the (ydirection, z direction)
				// x-ray propogation is always in the x direction



// Internal Variables
	wave /t logbook // keep a log of progress in the simulation
	variable progressbar
	variable timer
	variable timesofar

// Scattering internal calculations

	wave /d/c n //a single complex number for each material which will be set at each energy step, and used in calculations
	wave /d /c px // induced polarization in the x direction for real space, x,y,z
	wave /d /c pxFFT // fourier transform of above, in qx, qy, qz space
	wave /d /c py // same in y direction
	wave /d /c pyFFT// fourier transform of above, in qx, qy, qz space
	wave /d /c pz // same in z direction
	wave /d /c pzFFT// fourier transform of above, in qx, qy, qz space
	wave /d /c Escatx, Escaty, Escatz // 
	wave /d EscatSqr // the sum of the magsquared values as described above, before interpolation
	
	

	
// scattering outputs for each step
	wave/d scatter3D // 2d scattering simulation from 3D data
	wave/d int3D // radially summed intensity
	wave/d int3Dp0 // radially averaged intensity
	wave/d int3Dpara // radially summed intensity in direction of electric field
	wave/d int3Dp0para // radially averaged intensity in direction of electric field
	wave/d int3Dperp // radially summed intensity in direction normal to the electric field
	wave/d int3Dp0perp // radially averaged intensity in direction normal to the electric field
	
	
// scattering outputs saved (if chosen) and output for all energies
	wave/d Scatter3DSave // Full 2D Scattering Scattering Pattern Saved vs Energy
	wave /d int3DvsEn // same but for the 3D calculated data
	wave /d ratio3DvsEn // same but for the 3D calculated data
	wave /d Para3DvsEn // same but for the 3D calculated data
	wave /d Perp3DvsEn // same but for the 3D calculated data
	
// internal scattering variables
	variable step // current energy step (integer for storing data)
	variable en // current energy = 1239.84197 / wavelength =  ( 1239.84197/(2*pi)) * wavevector
	variable wavelength // current wavelength = 1239.84197 / energy = 2*pi / wavevector (in nanometers)
	variable k // current wavevector magnitude (2*pi / wavelength) = (2 * pi / 1239.84197 ) * Energy (in inverse nanometers)
					// current code has the X-rays propotaging in the x direction
	
endstructure

function addtologbook(s3d,lognote)
	struct ThreeDSystem &s3d
	string lognote
	s3d.timesofar +=stopmstimer(s3d.timer)/10^6
	s3d.timer = startmstimer
	variable loc = dimsize(s3d.logbook,0)
	insertpoints /M=0 loc,1, s3d.logbook
	s3d.logbook[loc][0] = time2str(s3d.timesofar)
	s3d.logbook[loc][1] = lognote
	make /t /free timew = {time2str(s3d.timesofar)}
	make /t /free notew = {lognote}
	ListBox progresslistbox win=SIM_STATUS,row=max(0,loc-8)
	updateprogress(s3d)
	pathinfo save3dscattering
	if(v_flag)
		Save/A=2/J/M="\n"/O/F/p=Save3DScattering timew,notew as s3d.name+".log"
	endif
end

function updateprogress(s3d)
	struct ThreeDSystem &s3d
	nvar  pbar = $(s3d.path+"progressbar")
	pbar = s3d.progressbar
	doupdate
end

function updatephase(s3d,phase)
	struct ThreeDSystem &s3d
	string phase
	SetDrawLayer/w=SIM_Status /k UserBack
	SetDrawEnv/w=SIM_Status fsize= 16,textxjust= 1,textyjust= 1
	DrawText/w=SIM_Status 425,15,"Simulation: "+s3d.name+" using "+ s3d.modelname + " morphology.  Progress displayed : " + phase
end
	

function model3D(modelname,voxelsize,sizescale,resolution,thickness,paramstring,materiallist,materialalignmentstring,efield,energymin,energymax,energynum[movie, save3d,moviepath,enwave,rotatesys])
	// electric field is always in the y z plane
	// modelname - a string with the name of the model  accepted models are specificed below
// morphology parameters	
	// sizsescale - basic parameter of all models, incidating general size scale
	// resolution - basic parameter of all models, indicating the grid size of the system (the pixel is the basic size)
	// paramstring - a model specific string containing all parameters needed to create the system
//alignment parameters	
	// materiallist - a list of all the materials in the system.  the number of materials needs to match the model requirements seperated by ","
	// materialalignment - a list the length of material list (material1alignment;material2alignment ...)
		//each materialalignment item is a comma seperated list of alignment parameters
			//( alignment (0 faceon, 1 edgeon, other = none)  if none, then whatever is produced by morphology is passed through (this may include alignment) , 
			//  volfrac (negative means as calculated by morphology) , 
			//  ratio of entropy cost of alignment vs interface alignment)
//scattering parameters
	// efield - a string containing a vector (y component, z component)
	// energymin - the minimum x-ray energy to calculate
	// energymax - the maximum x-ray energy to calculate
	// energynum - the number of energysteps
//simulation options (optional)
	//movie = make a movie of the process?
	// save3d - save the full scattering pattern at each energy
	// moviepath - avoid a popup if you supply a valid path for the movie
	// enwave - overwrites energy min, max and num and just calculated the energy values in this wave, if it exists
	// boolean to rotate the system - greatly increases the quality o he calculation, can be done in he fuure by roating the efield, rather than the elaborate way it is done now
	string modelname,paramstring,materiallist,materialalignmentstring,efield, moviepath
	variable sizescale, resolution, energymin,energymax,energynum, movie, save3d, thickness,voxelsize, rotatesys
	wave/z enwave // this will override the energy start, step etc, with arbitrary energy steps

		
	struct ThreeDSystem s3D
	variable result, timer1
	for(timer1=0;timer1<11;timer1+=1)
		s3d.timesofar = stopmstimer(timer1)
	endfor
	s3d.voxelsize = voxelsize
	s3d.timesofar = 0
	s3d.timer = startmstimer
	s3d.thickness = thickness
	s3d.num = resolution
	s3d.size = sizescale
	s3d.materials = materiallist
	s3d.efield = efield
	s3d.rot = paramisdefault(rotatesys)? 0 : rotatesys
	s3d.modelname = modelname
	s3d.paramstring = paramstring
	s3d.movie = movie
	s3d.materialalignmentstring = materialalignmentstring
	s3d.materialnum = itemsinlist(materiallist,",")
	s3d.progressbar = 0
	s3d.name = getdatafolder(0)
	string dfolder = getdatafolder(1)
	make /n=(0,2)/o /t Progresslist
	wave /t s3d.logbook = Progresslist
	variable /g progressbar
	s3d.path = getdatafolder(1)
	killwindow/z SIM_Status
	NewPanel /k=1/n=SIM_Status/W=(800,767,1650,1014)  as getdatafolder(0) + " Progress"
	updatephase(s3d,"Starting Up")
	ValDisplay valdispProgress,win=SIM_Status,limits={0,100,0},barmisc={10,1},highColor= (16385,16388,65535),pos={37.00,39.00},size={747.00,41.00},bodyWidth=747
	ValDisplay valdispProgress,win=SIM_Status,frame=4,value=#(dfolder+"progressbar")
	ListBox progresslistbox,win=SIM_STATUS,pos={18.00,98.00},size={796.00,134.00},listWave=$(dfolder+"Progresslist"),widths={57,688}
	if(paramisdefault(moviepath) || strlen(moviepath)<3)
			open /D/F=".mov" fileref as getdatafolder(0)+".mov"
			//close fileref
			moviepath = S_filename
	endif
	newpath /o/z/q Save3DScattering, parseFilePath(1,moviepath,":",1,0)
	if(!v_flag)
		Save/J/M="\n"/O/F/p=Save3DScattering s3d.logbook as s3d.name+".log"
	endif
	addtologbook(s3d,"Starting Up - Using " + s3d.modelname + " Morphology")
	addtologbook(s3d,"parameters: " + s3d.paramstring)
	addtologbook(s3d,"Voxelsize: " + num2str(s3d.voxelsize))
	addtologbook(s3d,"thickness: " + num2str(thickness))
	addtologbook(s3d,"resolution: " + num2str(resolution))
	addtologbook(s3d,"sizescale: " + num2str(sizescale))
	addtologbook(s3d,"materials: " + materiallist)
	addtologbook(s3d,"efield: along Z")
	addtologbook(s3d,"modelname: " + modelname)
	addtologbook(s3d,"movie: " + num2str(movie))
	addtologbook(s3d,"Material Slignment String: " + materialalignmentstring)
	
	string cmdstr = "Model3D("
	cmdstr += "\"" + modelname + "\","
	cmdstr += num2str(voxelsize) + ","
	cmdstr += num2str(sizescale) + ","
	cmdstr += num2str(resolution) + ","
	cmdstr += num2str(thickness) + ","
	cmdstr += "\"" + paramstring + "\","
	cmdstr += "\"" + materiallist + "\","
	cmdstr += "\"" + materialalignmentstring + "\","
	cmdstr += "\"" + efield + "\","
	cmdstr += num2str(energymin) + ","
	cmdstr += num2str(energymax) + ","
	cmdstr += num2str(energynum) + ","
	if(!paramisdefault(movie))
		cmdstr += "movie=" + num2str(movie) + ","
	endif
	if(!paramisdefault(save3d))
		cmdstr += "save3d=" + num2str(save3d) + ","
	endif
	if(!paramisdefault(moviepath))
		cmdstr += "\r\tmoviepath=\"" + moviepath + "\",\r\t	"
	endif
	if(!paramisdefault(enwave))
		cmdstr += "enwave=" + nameofwave(enwave) + ","
	endif
	if(!paramisdefault(rotatesys))
		cmdstr += "rotatesys=" + num2str(rotatesys) + ","
	endif
	cmdstr = removeending(cmdstr,",")
	cmdstr += ")"
	
	
	addtologbook(s3d,"Command to run again: " + replacestring("\r\t",cmdstr,""))
	print/len=1000 cmdstr
	
	if(paramisdefault(enwave))
		make /n=(energynum)/o enwave=energymin + p * (energymax-energymin)/(energynum-1)
	else
		if(waveexists(enwave))
			energynum = numpnts(enwave)
			energymin = wavemin(enwave)
			energymax = wavemax(enwave)
		else
			make /n=(energynum)/o enwave=energymin + p * (energymax-energymin)/(energynum-1)
		endif
	endif
	
	duplicate/o enwave, enwavedisp
	insertpoints 0,1,enwavedisp
	enwavedisp[0]=enwave[0]-.1
	
	
	 // set up the display and movie output
	 
	dowindow /k simulation_layout
	execute("Simulation_Layout()")
	if(s3d.movie)
		Variable DebugEnab
		DebuggerOptions
		DebugEnab = V_debugOnError //check for debug on error
		if (DebugEnab) //if it is on,
			DebuggerOptions debugOnError=0 //turn it off
		endif
		try
			closemovie;AbortOnRTE
		catch
			Variable err = GetRTError(1)	
			addtologbook(s3d," - Creating new movie")
		endtry
		if(DebugEnab)
			DebuggerOptions debugOnError=1
		endif
		SavePict/O/WIN=Simulation_Layout /E=-5/P=_PictGallery_/w=(0,0,800,800) as "SimPict"
		variable fileref
		
		//print moviepath
		newmovie /z /PICT=SimPict /O as moviepath
		if(v_flag)
			s3d.movie=0
			addtologbook(s3d,"movie could not be started")
		else
			addtologbook(s3d,"movie started")
		endif
	endif
	
	// try to call the chosen morphology creator
	
	if(exists("model3d_"+modelname)==6)
		updatephase(s3d,"Building Morphology")
		funcref model3d_existing creatfunc=$("model3d_"+modelname)
		creatfunc(s3d)
	else
		print "no recognized model could be created"
		return -1
	endif
	variable alignmentincluded
	if(exists("special_"+modelname)==6) // this is a unique function which if exists, will give some extra information about the model, in this case wether alignment is created
												// in the model, or it needs to be calculated - with modern models it is almos always generated within the model creation process
		funcref special_existing specfunc=$("special_"+modelname)
		alignmentincluded = stringmatch(specfunc(),"*IncludesAlignment*")
	endif
	
	//calculate alignment of system if needed - some morphologies already do this (code above checkes for this)
	
	if( stringmatch(materialalignmentstring,"none") )
		//Print "Time : "+time2str(s3d.timesofar) +"  -   Loading Existing material Alignment"
		
		addtologbook(s3d,"Loading Existing Material Alignment")
		updatephase(s3d,"Loading Alignment")
		loadexistingmaterialalignment(s3d)
	elseif(  !alignmentincluded   )
		//Print "Time : "+time2str(s3d.timesofar) +"  -   Calculating alignment in system"
		updatephase(s3d,"Calculating Alignment")
		
		addtologbook(s3d,"Calculating alignment in system")
		Align3Dsystem(s3d)
	endif
	addtologbook(s3d,"Alignment Loaded - Checking System Density for consistancy")
	
	// check the m1-m5 matrices to make sure they make sense
	
	updatephase(s3d,"Checking Morphology and Alignment")
	if(sum3dsystem(s3d)<0)
		addtologbook(s3d,"warning : System Creation failed to conserve density throughout system")
		//Print "Time : "+time2str(s3d.timesofar) +"  -   warning : System Creation failed to conserve density throughout system"
	elseif(sum3dsystem(s3d)>0)
		addtologbook(s3d,"Error : System Creation failed no density matrices are found")
		//print  "Time : "+time2str(s3d.timesofar) +"  -   Error : System Creation failed no density matrices are found"
		return -1
	else
		addtologbook(s3d,"Check of system density checks out")
		//Print "Time : "+time2str(s3d.timesofar) +"  -  Check of system density checks out"
	endif
	
	
	
	
	
	// At this point we have all the values that we need.  we could save the model now to disk, or load a model from disk at this point
	// data to save to model file
	// m1-m5 (if they exist)
		// first three dimensions are the aligned component - should go into mat_(1-5)_alignment vectors field
		// fourth dimension is the unaligned component should go into mat_(1-5)_unaligned scalar field

	variable h5file, groupid
	HDF5CreateFile /O h5file as replacestring("mov",replacestring("mp4",moviepath,"hd5"),"hd5")
	HDF5CreateGroup h5file , "vector_morphology" , groupID
	
	newdatafolder /o/s hd5output
	make /o /n=(dimsize(s3d.m1,0),dimsize(s3d.m1,1),dimsize(s3d.m1,2),3) /d Mat_1_alignment = s3d.m1[p][q][r][t]
	hdf5saveData /GZIP = {9,1} /LAYO={2,1,64,64,1} /MAXD={256,2048,2048,3} Mat_1_alignment, groupiD
	make /o /n=(dimsize(s3d.m1,0),dimsize(s3d.m1,1),dimsize(s3d.m1,2)) /d Mat_1_unaligned = s3d.m1[p][q][r][3]
	hdf5saveData /GZIP = {9,1} /LAYO={2,64,64,64} /MAXD={256,2048,2048} Mat_1_unaligned, groupiD
	if(waveexists(s3d.m2))
		make /o /n=(dimsize(s3d.m2,0),dimsize(s3d.m2,1),dimsize(s3d.m2,2),3) /d Mat_2_alignment = s3d.m2[p][q][r][t]
		make /o /n=(dimsize(s3d.m2,0),dimsize(s3d.m2,1),dimsize(s3d.m2,2)) /d Mat_2_unaligned = s3d.m2[p][q][r][3]
		hdf5saveData /GZIP = {9,1} /LAYO={2,1,64,64,1} /MAXD={256,2048,2048,3} Mat_2_alignment, groupiD
		hdf5saveData /GZIP = {9,1} /LAYO={2,64,64,64} /MAXD={256,2048,2048} Mat_2_unaligned, groupiD
		if(waveexists(s3d.m3))
			make /o /n=(dimsize(s3d.m3,0),dimsize(s3d.m3,1),dimsize(s3d.m3,2),3) /d Mat_3_alignment = s3d.m3[p][q][r][t]
			make /o /n=(dimsize(s3d.m3,0),dimsize(s3d.m3,1),dimsize(s3d.m3,2)) /d Mat_3_unaligned = s3d.m3[p][q][r][3]
			hdf5saveData /GZIP = {9,1} /LAYO={2,1,64,64,1} /MAXD={256,2048,2048,3} Mat_3_alignment, groupiD
			hdf5saveData /GZIP = {9,1} /LAYO={2,64,64,64} /MAXD={256,2048,2048} Mat_3_unaligned, groupiD
			if(waveexists(s3d.m4))
				make /o /n=(dimsize(s3d.m4,0),dimsize(s3d.m4,1),dimsize(s3d.m4,2),3) /d Mat_4_alignment = s3d.m4[p][q][r][t]
				make /o /n=(dimsize(s3d.m4,0),dimsize(s3d.m4,1),dimsize(s3d.m4,2)) /d Mat_4_unaligned = s3d.m4[p][q][r][3]
				hdf5saveData /GZIP = {9,1} /LAYO={2,1,64,64,1} /MAXD={256,2048,2048,3} Mat_4_alignment, groupiD
				hdf5saveData /GZIP = {9,1} /LAYO={2,64,64,64} /MAXD={256,2048,2048} Mat_4_unaligned, groupiD
				if(waveexists(s3d.m5))
					make /o /n=(dimsize(s3d.m5,0),dimsize(s3d.m5,1),dimsize(s3d.m5,2),3) /d Mat_5_alignment = s3d.m5[p][q][r][t]
					make /o /n=(dimsize(s3d.m5,0),dimsize(s3d.m5,1),dimsize(s3d.m5,2)) /d Mat_5_unaligned = s3d.m5[p][q][r][3]
					hdf5saveData /GZIP = {9,1} /LAYO={2,1,64,64,1} /MAXD={256,2048,2048,3} Mat_5_alignment, groupiD
					hdf5saveData /GZIP = {9,1} /LAYO={2,1,64,64} /MAXD={256,2048,2048} Mat_5_unaligned, groupiD
				endif
			endif
		endif
	endif
	
	newdatafolder /o/s igor_paramaters
	
	// parameters to save to parameter file 																				(*unique to Igor)
	variable /g igorrotation =  s3d.rot // boolean saying if we are rotating the system or not (90 degrees is always included) 		* - this should be universal
	variable /g igornum =  s3d.num // resolution of simulation (in each of the long axes of the film) 								* - the data itself should have this
	variable /g igorthickness =  s3d.thickness // voxels of thickness along the beam direction 												* - the data itself should have this
	variable /g igorvoxelsize =  s3d.voxelsize // the size of each pixel in the simulation (in nanometers) 									* - can be changed trivially, it is just a scaling parameter
	variable /g igormaterialnum =  s3d.materialnum // number of materials																				* - can be read from materials list
	string /g igorparamstring =  s3d.paramstring //specific parameters for creating density matricies according to the model				* - used for morphology creation
	string /g igormodelname =  s3d.modelname //store the model name somewhere																	* - used for morphology creation
	string /g igorpath =  s3d.path // the path to put important things for this calculation	
	string /g igorname =  s3d.name // name of simulation
	variable /g igormovie =  s3d.movie // are we doing a movie																					* - for visualization and some feedback
	string /g igormaterials =  s3d.materials // list of materials
	string /g igorefield =  s3d.efield // direction of efield - in the future, if combined with rotation, this doesn't need to be defined, maybe - all orientations will be calculated
	setdatafolder ::
	newdatafolder /o/s morphology_variables
	make /o /n=3 film_normal = {1,0,0}
	variable /g voxel_size_nm = s3d.voxelsize
	string /g morphology_creator = "EG_"+s3d.modelname
	string /g version = secs2date(datetime,-2)
	string /g creation_date = secs2date(datetime,-2) +" " + secs2time(datetime,1)
	string /g name = s3d.name
	setdatafolder ::
	
	
	HDF5CreateGroup h5file , "igor_parameters" , groupID
	HDF5SaveGroup /L=7 /O /R igor_paramaters , h5file , "igor_parameters" 
	HDF5CreateGroup h5file , "morphology_variables" , groupID
	HDF5SaveGroup /L=7 /O /R morphology_variables , h5file , "morphology_variables" 
	HDF5CloseFile h5file
	
	setdatafolder ::
	
	killdatafolder /z hd5output
	
	
	
	// make waves to store integrations and ratios
	
	make/d/o/n=(floor(s3d.num/sqrt(2)),energynum) int3DvsEn=0,ratio3DvsEn=0, para3dvsen, perp3dvsen
	make/d/o/n=(s3d.num,s3d.num,energynum) scatter3DSave
	wave s3d.scatter3DSave=scatter3DSave, s3d.int3DvsEn=int3DvsEn, s3d.ratio3DvsEn=ratio3DvsEn, s3d.para3dvsen=para3dvsen, s3d.perp3dvsen=perp3dvsen
	setscale /p x,pi/s3d.num,pi/s3d.num, s3d.int3DvsEn,s3d.ratio3DvsEn, s3d.para3dvsen, s3d.perp3dvsen // in physics units of q
	setscale /i x, -pi/s3d.voxelsize,pi/s3d.voxelsize, s3d.scatter3DSave
	setscale /i y, -pi/s3d.voxelsize,pi/s3d.voxelsize, s3d.scatter3DSave
	setscale /i y, energymin, energymax, s3d.int3DvsEn,s3d.ratio3DvsEn, s3d.para3dvsen, s3d.perp3dvsen // only valid if not using a energy wave // best to not used it, and loot at the wave
	setscale /i z, energymin, energymax, s3d.scatter3DSave // only valid if not using a wave, best to just use the wave
	
	
	variable en, j
	variable enstep = (energymax-energymin)/(energynum-1)
	variable angle=0, anglestep=2 //(1 degree steps from 0 to 180) - 180 steps
	if(!s3d.rot)
		anglestep = 90
	endif
	variable offsetx, offsety // offset for rotation to get rid of speckle
	s3d.step=0
	updatephase(s3d,"Calculating Scattering")
	
	
	// expand the M matrices (this will allow rotation, although it assumes periodic type boundary conditions, otherwise this will add scattering artifacts.)
	
	
	make /n=(dimsize(s3d.m1,0),dimsize(s3d.m1,1)*1.8,dimsize(s3d.m1,2)*1.8,dimsize(s3d.m1,3))/o m1ext = s3d.m1[p][mod(q+.9*dimsize(s3d.m1,1),dimsize(s3d.m1,1))][mod(r+.9*dimsize(s3d.m1,2),dimsize(s3d.m1,2))][s]
	duplicate/o s3d.m1, m1save
	if(s3d.materialnum>1)
		make /n=(dimsize(s3d.m2,0),dimsize(s3d.m2,1)*1.8,dimsize(s3d.m2,2)*1.8,dimsize(s3d.m2,3))/o m2ext = s3d.m2[p][mod(q+.9*dimsize(s3d.m2,1),dimsize(s3d.m2,1))][mod(r+.9*dimsize(s3d.m2,2),dimsize(s3d.m2,2))][s]
		duplicate/o s3d.m2, m2save
	endif
	if(s3d.materialnum>2)
		make /n=(dimsize(s3d.m3,0),dimsize(s3d.m3,1)*1.8,dimsize(s3d.m3,2)*1.8,dimsize(s3d.m3,3))/o m3ext = s3d.m3[p][mod(q+.9*dimsize(s3d.m3,1),dimsize(s3d.m3,1))][mod(r+.9*dimsize(s3d.m3,2),dimsize(s3d.m3,2))][s]
		duplicate/o s3d.m3, m3save
	endif
	if(s3d.materialnum>3)
		make /n=(dimsize(s3d.m4,0),dimsize(s3d.m4,1)*1.8,dimsize(s3d.m4,2)*1.8,dimsize(s3d.m4,3))/o m4ext = s3d.m4[p][mod(q+.9*dimsize(s3d.m4,1),dimsize(s3d.m4,1))][mod(r+.9*dimsize(s3d.m4,2),dimsize(s3d.m4,2))][s]
		duplicate/o s3d.m4, m4save
	endif
	if(s3d.materialnum>4)
		make /n=(dimsize(s3d.m5,0),dimsize(s3d.m5,1)*1.8,dimsize(s3d.m5,2)*1.8,dimsize(s3d.m5,3))/o m5ext = s3d.m5[p][mod(q+.9*dimsize(s3d.m5,1),dimsize(s3d.m5,1))][mod(r+.9*dimsize(s3d.m5,2),dimsize(s3d.m5,2))][s]
		duplicate/o s3d.m5, m5save
	endif
	
	
	createsimlayout() // open the live display for the scattering simulation
	
	
	for( j=0; j <dimsize(enwave,0) ; j+=1) // main loop through energies
		en=enwave[j]
		//s3d.progressbar = 100*(en-energymin) / (energymax-energymin)
		s3d.progressbar = 100*(j) / (dimsize(enwave,0))
		s3d.en = en
		s3d.step = j
		s3d.wavelength = 1239.84197/en
		s3d.k = 2*pi/s3d.wavelength
		addtologbook(s3d,"calculating polarization for " + num2str(s3d.en) + " eV")
		//Print "Time : "+time2str(s3d.timesofar) +"  -  Creating Indexwave for energy : " + num2str(s3d.en)
		
		
		// loop through polarization and scattering calculations as the system is rotated
		//  this assumes that the enlarged systems we have created were from periodic boundary conditions.  Unfortunately once rotated, the boundaries cannot maintain periodicity, and will introduce some artifacts\
		//  this is far outweighed by the benefit particulary when looking at anisotropy
		// For Future development****  instead of rotating the system, rotate the efield, and then rotate the scattering pattern back.  it is very important that sums are only made with efield pointing in the same direction
		//   otherwise all anisotropy will be removed (simulating circular polarization)
		for(angle=0;angle<180;angle+=anglestep)
			offsetx = 0//floor(enoise(dimsize(s3d.m1,1)/10)) // this was rotating and offsetting the center randomly, I found that the offsetting didn't matter
			offsety = 0//floor(enoise(dimsize(s3d.m1,2)/10))
			
			sysrotate(m1ext,s3d.m1,angle,offsetx,offsety)
			if(s3d.materialnum>1)
				sysrotate(m2ext,s3d.m2,angle,offsetx,offsety)
			endif
			if(s3d.materialnum>2)
				sysrotate(m3ext,s3d.m3,angle,offsetx,offsety)
			endif
			if(s3d.materialnum>3)
				sysrotate(m4ext,s3d.m4,angle,offsetx,offsety)
			endif
			if(s3d.materialnum>4)
				sysrotate(m5ext,s3d.m5,angle,offsetx,offsety)
			endif
			
			result = CalculatePolarizationWave(s3d) // This is the main calculation, taking efield and m matrices with npara and nperp and producing the induced polarization within each voxel
			if(result<0)
				addtologbook(s3d,"quitting because of critical errors in calculating polarization")
				//print "quitting because of critical errors in calculating polarization"
				return -1
			elseif(result>0 && s3d.en == energymin && angle==0)
				addtologbook(s3d,"warning - no aligned optical constants for "+stringfromlist(result-1, s3d.materials,",")+" were found, using normal Optical Constants")
				//print "warning - no aligned optical constants for "+stringfromlist(result-1, s3d.materials,",")+" were found, using normal Optical Constants"
			endif
	
		//	addtologbook(s3d,"Angle="+num2str(angle))
			//Print "Time : "+time2str(s3d.timesofar) +"  -  Creating 3D Scattering Pattern"
			Calculate3DScatter(s3d) // this projects the induced polarization to the far field, simulating the scattering which emits from this polarization field
			if(angle==0)
				duplicate/free s3d.scatter3d, scatter3dsum
			else
				scatter3dsum+= s3d.scatter3D
			endif
			if(angle/10 == round(angle/10) && angle + anglestep < 180)
				if(angle>0)
					s3d.scatter3D = scatter3Dsum/(angle/anglestep) // add up the scattering from different rotations
					
				else
					RadialIntegratesystem(s3d)
					killwindow /z simulation_layout
					createsimlayout()
					TextBox/w=Simulation_Layout/C/N=text5/F=0/A=RT/X=51.85/Y=37.15/G=(16386,65535,16385) "\\Z24"+num2str(s3d.en)+" eV"
					SetDrawLayer /w=RatioIntensity /k userfront
					SetDrawEnv /w=RatioIntensity ycoord= left,linefgc= (65280,0,0),dash= 1,linethick= 3.00;DelayUpdate
					DrawLine/w=RatioIntensity 0,s3d.en,1,s3d.en
					SetDrawLayer /w=ScatteringIntensity /k userfront
					SetDrawEnv /w=ScatteringIntensity ycoord= left,linefgc= (65280,0,0),dash= 1,linethick= 3.00;DelayUpdate
					DrawLine/w=ScatteringIntensity 0,s3d.en,1,s3d.en
				endif
				RadialIntegratesystem(s3d)
				ModifyGraph/w=ParaPerpInt hbFill(int3Dpara)=5
				doupdate
				dowindow /f simulation_layout
				if(s3d.movie)
					SavePict/O/WIN=Simulation_Layout /E=-5/P=_PictGallery_/w=(0,0,800,800) as "SimPict"
					AddMovieFrame/PICT=SimPict
				endif
			endif
		endfor
		s3d.scatter3D = scatter3Dsum/(180/anglestep)
		//imagerotate /o/H scatter3Dsum
		//s3d.scatter3D += scatter3Dsum/(180/anglestep)
		//imagerotate /o/V scatter3Dsum
		//s3d.scatter3D += scatter3Dsum/(180/anglestep)
		//imagerotate /o/H scatter3Dsum
		//s3d.scatter3D += scatter3Dsum/(180/anglestep)
		
		RadialIntegratesystem(s3d)

		
		
		StoreIntegrations(s3d)
		// visualizescatter(s3d) // updates displays of scattering profiles
		//(re)create layout of scattering for an energy
		//CreateSimLayout()
		// update the layout (projections, 2D scattering pattern, label, and cuts)
		//doupdate
		
		//doupdate
		if(s3d.movie)
			killwindow /z simulation_layout
		
		// make the display
			createsimlayout()
			TextBox/w=Simulation_Layout/C/N=text5/F=0/A=RT/X=51.85/Y=37.15/G=(16386,65535,16385) "\\Z24"+num2str(s3d.en)+" eV"
			SetDrawLayer /w=RatioIntensity /k userfront
			SetDrawEnv /w=RatioIntensity ycoord= left,linefgc= (65280,0,0),dash= 1,linethick= 3.00;DelayUpdate
			DrawLine/w=RatioIntensity 0,s3d.en,1,s3d.en
			SetDrawLayer /w=ScatteringIntensity /k userfront
			SetDrawEnv /w=ScatteringIntensity ycoord= left,linefgc= (65280,0,0),dash= 1,linethick= 3.00;DelayUpdate
			DrawLine/w=ScatteringIntensity 0,s3d.en,1,s3d.en
			SavePict/O/WIN=Simulation_Layout /E=-5/P=_PictGallery_/w=(0,0,800,800) as "SimPict"
			dowindow /f simulationlayout
			ModifyGraph/w=ParaPerpInt hbFill(int3Dpara)=2
			doupdate
			AddMovieFrame/PICT=SimPict
			AddMovieFrame/PICT=SimPict
			AddMovieFrame/PICT=SimPict
			AddMovieFrame/PICT=SimPict
			AddMovieFrame/PICT=SimPict
			AddMovieFrame/PICT=SimPict
			AddMovieFrame/PICT=SimPict
			AddMovieFrame/PICT=SimPict
		endif
		s3d.step +=1
	endfor
	if(s3d.movie)
		addtologbook(s3d,"Finishing Movie, and closing file")
		//Print "Time : "+time2str(s3d.timesofar) +"  -  Finishing Movie, and closing file"
		closemovie
	endif
	updatephase(s3d,"Complete")
	s3d.progressbar=100
	addtologbook(s3d,"Complete.")
	//killwindow /z SIM_Status
	//Print "Time : "+time2str(s3d.timesofar) +"  -  Complete."
	//display 2D rato and integrations vs q and energy
end
function createsimlayout()
	//dowindow /k Simulation_Layout
	debuggeroptions debugonerror=0
	ScatteringIntensitydisp()
	RatioIntensitydisp()
	ScatteringPatterndisp()
	ParaPerpIntdisp()
	DeltaProjectiondisp()
	BetaProjectiondisp()
	dowindow Simulation_Layout
	if(!v_flag)
		execute("Simulation_Layout()")
	endif
	dowindow /f Simulation_Layout
	doupdate
end

function align3dsystem(s3d) // given a morphology of only scalar density, produce an alignment field
	struct ThreeDSystem &s3d
	string matstr = s3d.materialalignmentstring
	string mat,alignment,alignmentwidth, alignmentstrength
	variable nummaterials = itemsinlist(matstr) , i
	variable align, width, strength
	variable volfrac
	if(nummaterials==2)
		mat = stringfromlist(0,matstr)
		if(stringmatch(mat,"*Face-on*") || stringmatch(mat, "1*"))
			align = 1
		elseif(stringmatch(mat, "*Edge-on*") || stringmatch(mat, "0*"))
			align=0
		else
			//don't do the next calculation
			align=-1
		endif
		if(align==-1)
			//don't do the calculation for m1
			make /o/n=(dimsize(s3d.density1,0),dimsize(s3d.density1,1),dimsize(s3d.density1,2),4) s3d.m1=0
			s3d.m1[][][][3] = s3d.density1[p][q][r]
		else
			splitstring/e="^[^,]*,([^,]*),([^,]*)" mat, alignmentwidth, alignmentstrength
			strength = str2num(alignmentstrength)
			width = str2num(alignmentwidth)
			s3d.timesofar += stopmstimer(s3d.timer)/1e6
			s3d.timer = startmstimer
			wave alignmentw = createalignmentdensity(s3d.density1,width,strength, align, s3d.movie, s3d.timesofar)
			//createbinarysystem(s3d.density1,width,strength, 1, 0, 1- volfrac, 0, 1,s=s)
			duplicate/o alignmentw, s3d.m1
			//s3d.m1 *=s3d.density1[p][q][r]
			s3d.m1 *=sqrt(s3d.density1[p][q][r])
			s3d.m1[][][][3] *=sqrt(s3d.density1[p][q][r])
		endif
		mat = stringfromlist(1,matstr)
		if(stringmatch( mat,"*Face-on*")  || stringmatch(mat, "1*") )
			align = 1
		elseif(stringmatch( mat,"*Edge-on*") || stringmatch(mat, "0*") )
			align=0
		else
			duplicate/o s3d.m1, s3d.m2
			s3d.m2[][][][0,2]=0
			s3d.m2[][][][3] = 1 - s3d.m1[p][q][r][0]^2 - s3d.m1[p][q][r][1]^2 - s3d.m1[p][q][r][2]^2 -s3d.m1[p][q][r][3] // all density that isn't accounted for in m1 is unaligned m2
			align=-1
		endif
		if(align==-1)
			// we've already created both m1 and m2, so we're done!
		else
			splitstring/e="^[^,]*,([^,]*),([^,]*)" mat, alignmentwidth, alignmentstrength
			strength = str2num(alignmentstrength)
			width = str2num(alignmentwidth)
			volfrac = mean(s3d.density1)
			s3d.timesofar += stopmstimer(s3d.timer)/1e6
			s3d.timer = startmstimer
			wave alignmentw = createalignmentdensity(s3d.density1,width,strength, align, s3d.movie, s3d.timesofar)
		//	wave alignmentw = root:packages:ScatterSim3D:test:sn
			//createbinarysystem(s3d.density1,width,strength, 1, 0, 1- volfrac, 0, 1,s=s)
			duplicate/o alignmentw, s3d.m2
			s3d.m2 *=sqrt(1-s3d.density1[p][q][r])
			s3d.m2[][][][3] *=sqrt(1-s3d.density1[p][q][r])
		endif
	else
		for(i=0;i<nummaterials;i+=1)
			mat = stringfromlist(i,matstr)
			
			// make density and m matricies for all dimensions
			if(i==0)
				make /o/n=(dimsize(s3d.density1,0),dimsize(s3d.density1,1),dimsize(s3d.density1,2),4) s3d.m1=0
				s3d.m1[][][][3] = s3d.density1[p][q][r]
			elseif(i==1)
				duplicate /o s3d.m1,s3d.m2
				s3d.m2 = 0
				if(i< nummaterials - 1)
					s3d.m2[][][][3] = s3d.density2[p][q][r]
				else
					//this is the last material, so density wave doesn't exist, need to calculate it
					s3d.m2[][][][3] = 1 - s3d.m1[p][q][r][0]^2 - s3d.m1[p][q][r][1]^2 - s3d.m1[p][q][r][2]^2 -s3d.m1[p][q][r][3]
					duplicate /o s3d.density1,s3d.density2
					s3d.density2[][][] = s3d.m2[p][q][r][3]
				endif
			elseif(i==2)
				duplicate /o s3d.m1,s3d.m3
				s3d.m3 = 0
				if(i< nummaterials - 1)
					s3d.m3[][][][3] = s3d.density3[p][q][r]
				else
					//this is the last material, so density wave doesn't exist, need to calculate it
					s3d.m3[][][][3] = 1
					s3d.m3[][][][3] -= s3d.m1[p][q][r][0]^2 + s3d.m1[p][q][r][1]^2 + s3d.m1[p][q][r][2]^2 + s3d.m1[p][q][r][3]
					s3d.m3[][][][3] -= s3d.m2[p][q][r][0]^2 + s3d.m2[p][q][r][1]^2 + s3d.m2[p][q][r][2]^2 + s3d.m2[p][q][r][3]
					duplicate /o s3d.density1,s3d.density3
					s3d.density3[][][] = s3d.m3[p][q][r][3]
				endif
			elseif(i==3)
				duplicate /o s3d.m1,s3d.m4
				s3d.m4 = 0
				if(i< nummaterials - 1)
					s3d.m4[][][][3] = s3d.density4[p][q][r]
				else
					//this is the last material, so density wave doesn't exist, need to calculate it
					s3d.m4[][][][3] = 1
					s3d.m4[][][][3] -= s3d.m1[p][q][r][0]^2 + s3d.m1[p][q][r][1]^2 + s3d.m1[p][q][r][2]^2 + s3d.m1[p][q][r][3]
					s3d.m4[][][][3] -= s3d.m2[p][q][r][0]^2 + s3d.m2[p][q][r][1]^2 + s3d.m2[p][q][r][2]^2 + s3d.m2[p][q][r][3]
					s3d.m4[][][][3] -= s3d.m3[p][q][r][0]^2 + s3d.m3[p][q][r][1]^2 + s3d.m3[p][q][r][2]^2 + s3d.m3[p][q][r][3]
					duplicate /o s3d.density1,s3d.density4
					s3d.density4[][][] = s3d.m4[p][q][r][3]
				endif
			elseif(i==4)
				duplicate /o s3d.m1,s3d.m4
				s3d.m4 = 0
				if(i< nummaterials - 1)
					s3d.m5[][][][3] = s3d.density5[p][q][r] // this shouldn't exist, we are limited to 5 materials so 4 densities in the current code
				else
					//this is the last material, so density wave doesn't exist, need to calculate it
					s3d.m5[][][][3] = 1
					s3d.m5[][][][3] -= s3d.m1[p][q][r][0]^2 + s3d.m1[p][q][r][1]^2 + s3d.m1[p][q][r][2]^2 + s3d.m1[p][q][r][3]
					s3d.m5[][][][3] -= s3d.m2[p][q][r][0]^2 + s3d.m2[p][q][r][1]^2 + s3d.m2[p][q][r][2]^2 + s3d.m2[p][q][r][3]
					s3d.m5[][][][3] -= s3d.m3[p][q][r][0]^2 + s3d.m3[p][q][r][1]^2 + s3d.m3[p][q][r][2]^2 + s3d.m3[p][q][r][3]
					s3d.m5[][][][3] -= s3d.m4[p][q][r][0]^2 + s3d.m4[p][q][r][1]^2 + s3d.m4[p][q][r][2]^2 + s3d.m4[p][q][r][3]
					duplicate /o s3d.density1,s3d.density5
					s3d.density5[][][] = s3d.m5[p][q][r][3]
				endif
			endif
			
			
			
			if(stringmatch("*Face-on*", mat)  || stringmatch(mat, "1*"))
				align = 1
			elseif(stringmatch("*Edge-on*", mat) || stringmatch(mat, "0*"))
				align=0
			else
				//don't do the calculation for this material, the previously made density and m matrix is sufficient
				continue // move to next material
			endif
			splitstring/e="^[^,]*,([^,]*),([^,]*)" mat, alignmentwidth, alignmentstrength
			strength = str2num(alignmentstrength)
			width = str2num(alignmentwidth)
			volfrac = mean(s3d.density1)
			s3d.timesofar += stopmstimer(s3d.timer)/1e6
			s3d.timer = startmstimer
			if(i==0)
				wave alignmentw = createalignmentdensity(s3d.density1,width,strength, align, s3d.movie, s3d.timesofar)
				duplicate/o alignmentw, s3d.m1
				
				s3d.m1 *=sqrt(s3d.density1[p][q][r])
				s3d.m1[][][][3] *=sqrt(s3d.density1[p][q][r]) // last dimension is density not sqrt density, so it needs an extra multiplication
				
			elseif(i==1)
				wave alignmentw = createalignmentdensity(s3d.density2,width,strength, align, s3d.movie, s3d.timesofar)
				duplicate/o alignmentw, s3d.m2
				
				s3d.m2 *=sqrt(s3d.density2[p][q][r])
				s3d.m2[][][][3] *=sqrt(s3d.density2[p][q][r])
			elseif(i==2)
				wave alignmentw = createalignmentdensity(s3d.density3,width,strength, align, s3d.movie, s3d.timesofar)
				duplicate/o alignmentw, s3d.m3
				
				s3d.m3 *=sqrt(s3d.density3[p][q][r])
				s3d.m3[][][][3] *=sqrt(s3d.density3[p][q][r])
			elseif(i==3)
				wave alignmentw = createalignmentdensity(s3d.density4,width,strength, align, s3d.movie, s3d.timesofar)
				duplicate/o alignmentw, s3d.m4
				
				s3d.m4 *=sqrt(s3d.density4[p][q][r])
				s3d.m4[][][][3] *=sqrt(s3d.density4[p][q][r])
			elseif(i==4)
				wave alignmentw = createalignmentdensity(s3d.density5,width,strength, align, s3d.movie, s3d.timesofar)
				duplicate/o alignmentw, s3d.m5
				
				s3d.m5 *=sqrt(s3d.density5[p][q][r])
				s3d.m5[][][][3] *=sqrt(s3d.density5[p][q][r])
			endif
		
		endfor
	endif
end
function storeintegrations(s3d)
	struct ThreeDSystem &s3d
	s3d.int3DvsEn[][s3d.step] =  s3d.int3d[p] * x^2
	s3d.ratio3DvsEn[][s3d.step] = (s3d.int3dpara(x) - s3d.int3dperp(x)) / (s3d.int3dpara(x) + s3d.int3dperp(x))
	s3d.scatter3dSave[][][s3d.step] = s3d.scatter3d(x)(y)
	s3d.perp3Dvsen[][s3d.step] =  s3d.int3dperp(x)
	s3d.para3Dvsen[][s3d.step] = s3d.int3dpara(x)
	
end
function radialintegratesystem(s3d)
	//radially integrate the 2D and 3D transforms and populate the corresponding waves in s (see structure declaration for details)
	struct ThreeDSystem &s3d
	variable ey =0// str2num(stringfromlist(0,s3d.efield,","))
	variable ez =1// str2num(stringfromlist(1,s3d.efield,","))
	variable pe = atan(ey/ez)*180/pi
	variable pa = atan(-ez/ey)*180/pi
	wave s3d.int3dp0 = radialintegratew(s3d.scatter3d,0,90,"int3dp0")
	wave s3d.int3dp0para = radialintegratew(s3d.scatter3d,pa-20,pa+20,"int3dp0para")
	wave s3d.int3dp0perp = radialintegratew(s3d.scatter3d,pe-20,pe+20,"int3dp0perp")
	duplicate/o s3d.int3dp0, s3d.int3d
	duplicate/o s3d.int3dp0para, s3d.int3dpara
	duplicate/o s3d.int3dp0perp, s3d.int3dperp
//	s3d.int3d *= x^2
//	s3d.int3dpara *= x^2
//	s3d.int3dperp *= x^2
end



function Calculate3DScatter(s3d) // turn the polarization fields into scattering fields (real space to fourier space)
	struct ThreeDSystem &s3d
	variable thicknum = dimsize(s3d.px,0) , num = s3d.num
	if(num/2 != round(num/2) )
		num +=1
	endif
	if(thicknum/2 != round(thicknum/2) )
		thicknum +=1
	endif
	wave px = s3d.px, py=s3d.py, pz=s3d.pz
	
	// each element is fourier transformed seperatly, the different components only add incoherently
	
	FFT /WINF=Hanning /pad={1*thicknum,1*num,1*num} /Dest=pxfft px
	FFT /WINF=Hanning /pad={1*thicknum,1*num,1*num} /Dest=pyfft py
	FFT /WINF=Hanning /pad={1*thicknum,1*num,1*num} /Dest=pzfft pz
	wave/c s3d.pxfft = pxfft, s3d.pyfft=pyfft, s3d.pzfft=pzfft
	//make /n=(dimsize(s3d.pxfft,0),dimsize(s3d.pxfft,1),dimsize(s3d.pxfft,2),9) /d/o s3d.qtensor
	make /n=(dimsize(s3d.pxfft,0),dimsize(s3d.pxfft,1),dimsize(s3d.pxfft,2)) /d/o s3d.EscatSqr =0
	wave EscatSqr=s3d.EscatSqr
	setscale /i x, -pi/s3d.voxelsize, pi/s3d.voxelsize, pxfft, pyfft, pzfft, EscatSqr // changing from math units (1/voxelsize) to physics units (2pi/voxelsize)
	setscale /i y, -pi/s3d.voxelsize, pi/s3d.voxelsize, pxfft, pyfft, pzfft, EscatSqr
	setscale /i z, -pi/s3d.voxelsize, pi/s3d.voxelsize, pxfft, pyfft, pzfft, EscatSqr // all three dimensions have the same range.  The x dimension may have a lower number of pixels, but they cover the same range

	variable k = s3d.k
	multithread EscatSqr = magsqr(x * (2 * k + x) * pxfft[p][q][r] + y * (k + x) * pyfft[p][q][r] + z * (k + x) * pzfft[p][q][r]) // Equation 8 from paper (qx -> -qx ends up equivalent somehow - checked with mathematica)
	multithread EscatSqr += magsqr(y * (k + x)* pxfft[p][q][r]+ (y^2 - k^2) * pyfft[p][q][r] + y * z * pzfft[p][q][r])
	multithread EscatSqr += magsqr(z * (k + x)* pxfft[p][q][r]+ y * z* pyfft[p][q][r]+ (z^2 - k^2 ) * pzfft[p][q][r])
	make /d/n=(s3d.num,s3d.num)/o s3d.scatter3D // this will hold the final scattering result
	wave scatter3d = s3d.scatter3d
	setscale /i x, -pi/s3d.voxelsize, pi/s3d.voxelsize, s3d.scatter3D // changing from math units (1/voxelsize) to physics units (2pi/voxelsize)
	setscale /i y, -pi/s3d.voxelsize, pi/s3d.voxelsize, s3d.scatter3D
	multithread scatter3D = k^2 > x^2 + y^2 ? interp3d(EscatSqr,-k + sqrt(k^2 - x^2 -y^2),x,y) : nan // projection to the EWalds sphere
	
	duplicate /o s3d.EscatSqr, pxreal
	//pxreal = magsqr(s3d.pxfft)
	imagetransform/meth=2 xprojection, pxreal
	duplicate/o m_xprojection, pxfftproj
	
	return 0
end

function CalculatePolarizationWave(s3d)
	struct ThreeDSystem &s3d
	variable result
	//create polarization waves
	make /o/c/d/n=(dimsize(s3d.m1,0),s3d.num,s3d.num) s3d.px = 0, s3d.py = 0, s3d.pz = 0
	setscale /p x, dimoffset(s3d.m1,0), dimdelta(s3d.m1,0), s3d.px, s3d.py, s3d.pz
	setscale /p y, dimoffset(s3d.m1,1), dimdelta(s3d.m1,1), s3d.px, s3d.py, s3d.pz
	setscale /p z, dimoffset(s3d.m1,2), dimdelta(s3d.m1,2), s3d.px, s3d.py, s3d.pz
	wave/c px= s3d.px, py=s3d.py, pz=s3d.pz
	//for each materials in list
	string material
	variable i, matnum = itemsinlist(s3d.materials,",")
	variable delpara, betapara, delperp, betaperp,ey,ez
	variable /c nperp, npara
	ey = 0 //str2num(stringfromlist(0,s3d.efield,","))
	ez = 1 //str2num(stringfromlist(1,s3d.efield,","))  // the complication of adding depolarization, makes the equations only work with incident polarization in the z direction
	if(ez+ey != 1)
		print "error - Electric field components are invalid : "+ s3d.efield
		return -1
	endif
	
	for(i=0;i<matnum;i+=1) // loop through materials
		material = stringfromlist(i,s3d.materials,",")
		//if there exist both para and perp optical constants load the corresponding values for that energy, otherwise throw warning and load same OC for both
		if(ocsexist(material)==2)
			delpara = getdelta(material,1239.84197/s3d.en,aligned=1)
			betapara = getbeta(material,1239.84197/s3d.en,aligned=1)
			delperp = getdelta(material,1239.84197/s3d.en,aligned=0)
			betaperp = getbeta(material,1239.84197/s3d.en,aligned=0)
			result = 0
		elseif(ocsexist(material)==1)
			delpara = getdelta(material,1239.84197/s3d.en)
			betapara = getbeta(material,1239.84197/s3d.en)
			delperp = getdelta(material,1239.84197/s3d.en)
			betaperp = getbeta(material,1239.84197/s3d.en)
			result = i+1
		else
			print "error - optical constants for " + material+ " were not found"
			return -1
		endif
		nperp = 1-delperp + sqrt(-1) * betaperp
		npara = 1-delpara + sqrt(-1) * betapara
		// select the right material polarized density wave
		if(i==0)
			wave m = s3d.m1
		elseif(i==1)
			wave m = s3d.m2
		elseif(i==2)
			wave m = s3d.m3
		elseif(i==3)
			wave m = s3d.m4
		elseif(i==4)
			wave m = s3d.m5
		else
			print " Error - invalid number of materials specified"
			return -1
		endif
		multithread px[][][] +=   ( (npara^2 - nperp^2)/(4*pi) ) * m[p][q][r][2]* m[p][q][r][0] // equation 5 from paper
		multithread py[][][] +=   ( (npara^2 - nperp^2)/(4*pi) ) * m[p][q][r][2]* m[p][q][r][1] // equation 5 from paper
		
				// following lines are equation 4 from paper
		multithread pz[][][] -=   (m[p][q][r][0]^2 + m[p][q][r][1]^2 + m[p][q][r][2]^2 + m[p][q][r][3])/(4*pi) // the overall material (phi) fourth term
		multithread pz[][][] +=   m[p][q][r][3]/(36*pi) * (npara + 2*nperp)^2 // first term (unaligned portion)
		multithread pz[][][] +=   npara^2 * m[p][q][r][2]^2 /(4*pi) // term 2
		multithread pz[][][] +=   nperp^2 * (m[p][q][r][0]^2 + m[p][q][r][1]^2) /(4*pi) //term 3
		//multithread pz[][][] +=   (m[p][q][r][0]^2 + m[p][q][r][1]^2 + m[p][q][r][2]^2 + m[p][q][r][3])/(4*pi) +  (m[p][q][r][3]/(12*pi)) * (npara^2 + 2*nperp^2) + (npara^2 * m[p][q][r][2]^2 + nperp^2 * (m[p][q][r][0]^2 +m[p][q][r][1]^2) )/(4*pi)
	endfor
	
	
	// make a bunch of projections to try to visualize what this looks like
	duplicate /o s3d.density1, preal, pimag
	preal = real(s3d.px)
	pimag = imag(s3d.px)
	imagetransform/meth=2 xprojection, preal
	duplicate/o m_xprojection, pxrealproj
	imagetransform/meth=2 xprojection, pimag
	duplicate/o m_xprojection, pximagproj
	preal = real(s3d.py)
	pimag = imag(s3d.py)
	imagetransform/meth=2 xprojection, preal
	duplicate/o m_xprojection, pyrealproj
	imagetransform/meth=2 xprojection, pimag
	duplicate/o m_xprojection, pyimagproj
	preal = real(s3d.pz)
	pimag = imag(s3d.pz)
	imagetransform/meth=2 xprojection, preal
	duplicate/o m_xprojection, pzrealproj
	imagetransform/meth=2 xprojection, pimag
	duplicate/o m_xprojection, pzimagproj
	return result
end

function sum3dsystem(s3d) // this sums all of the relative densities of all the materials and determines how much each voxel holds.
	struct ThreeDsystem &s3d
	make /o/n=(dimsize(s3d.m1,0),s3d.num,s3d.num) sumwave
	if(waveexists(s3d.m5))
		multithread sumwave = s3d.m1[p][q][r][0]^2 + s3d.m1[p][q][r][1]^2 + s3d.m1[p][q][r][2]^2 
		multithread sumwave += s3d.m2[p][q][r][0]^2 + s3d.m2[p][q][r][1]^2 + s3d.m2[p][q][r][2]^2 
		multithread sumwave += s3d.m3[p][q][r][0]^2 + s3d.m3[p][q][r][1]^2 +s3d.m3[p][q][r][2]^2  
		multithread sumwave += s3d.m4[p][q][r][0]^2 + s3d.m4[p][q][r][1]^2 +s3d.m4[p][q][r][2]^2  
		multithread sumwave += s3d.m5[p][q][r][0]^2 + s3d.m5[p][q][r][1]^2 +s3d.m5[p][q][r][2]^2 
		multithread sumwave += s3d.m1[p][q][r][3]+s3d.m2[p][q][r][3]+s3d.m3[p][q][r][3]+s3d.m4[p][q][r][3]+s3d.m5[p][q][r][3]
	elseif(waveexists(s3d.m4))
		multithread sumwave = s3d.m1[p][q][r][0]^2 + s3d.m1[p][q][r][1]^2 + s3d.m1[p][q][r][2]^2 
		multithread sumwave += s3d.m2[p][q][r][0]^2 + s3d.m2[p][q][r][1]^2 + s3d.m2[p][q][r][2]^2 
		multithread sumwave += s3d.m3[p][q][r][0]^2 + s3d.m3[p][q][r][1]^2 +s3d.m3[p][q][r][2]^2  
		multithread sumwave += s3d.m4[p][q][r][0]^2 + s3d.m4[p][q][r][1]^2 +s3d.m4[p][q][r][2]^2  
		multithread sumwave +=s3d.m1[p][q][r][3]+s3d.m2[p][q][r][3]+s3d.m3[p][q][r][3]+s3d.m4[p][q][r][3]
	elseif(waveexists(s3d.m3))
		multithread sumwave = s3d.m1[p][q][r][0]^2 + s3d.m1[p][q][r][1]^2 + s3d.m1[p][q][r][2]^2 + s3d.m2[p][q][r][0]^2 + s3d.m2[p][q][r][1]^2 + s3d.m2[p][q][r][2]^2 + s3d.m3[p][q][r][0]^2 + s3d.m3[p][q][r][1]^2 +s3d.m3[p][q][r][2]^2 +s3d.m1[p][q][r][3]+s3d.m2[p][q][r][3]+s3d.m3[p][q][r][3]
	elseif(waveexists(s3d.m2))
		multithread sumwave = s3d.m1[p][q][r][0]^2 + s3d.m1[p][q][r][1]^2 + s3d.m1[p][q][r][2]^2 + s3d.m2[p][q][r][0]^2 + s3d.m2[p][q][r][1]^2 + s3d.m2[p][q][r][2]^2 +s3d.m1[p][q][r][3]+s3d.m2[p][q][r][3]
	elseif(waveexists(s3d.m1))
		multithread sumwave = s3d.m1[p][q][r][0]^2 + s3d.m1[p][q][r][1]^2 + s3d.m1[p][q][r][2]^2 +s3d.m1[p][q][r][3]
	else
		return 1
	endif
	
	
	if(!(wavemax(sumwave)<1.01) || !(wavemin(sumwave) >.99))
		killwaves /z sumwave
		return -1
	endif
	killwaves /z sumwave
	return 0
end

function /s variables_spheres2()
	return "Interpenetration [pixels],SetVariable,2;Minimum Seperation,SetVariable,.1;Number of Particles (Max),SetVariable,500;PolyDispursity (sigma of radiuses),setvariable,5;Noise,SetVariable,0;"
end
function model3D_Spheres2(s3d)
	//Creates a spherical system, with two components, aligned 
	// utilizes extra parameters, number of particles, polydispursity (a sigma of radiuses) and populates the 3D space with 
		//aligned particles
	// material 1 is insize the spheres, material 2 is outside (although the concentrations are determined by m2volfrac
	//paramstring = 
		//MORPHOLOGY parameters
		//0 interpenetration of spheres (roughness smaller than wavelength)
		//1 minimum seperation (when placing a new sphere into the system, it needs to be this far away from any others)
			//	this can be a fraction of the radius (if less than 1) or a flat number of pixels (if greater than 1)
		//2 number of particles // the maximum number of particles . the program will stop when it cannot fit anymore
		//3 polydispursity (sigma of radiuses)
		//4 thickness (this can be different than size)
		//5 noise (percentage of vacuum - density variations)
	//num = resolution (in each of 3 dimensions)
	//size = average size of sphere (radius)
	
	struct ThreeDSystem &s3d
	if(itemsinlist(s3d.paramstring,",")<5)
		return -1
	endif
	newdatafolder /o/s CreatingSpheres
	variable interpenetration = 	str2num(stringfromlist( 0 ,s3d.paramstring,","))
	variable minsep = 			str2num(stringfromlist( 1 ,s3d.paramstring,","))
	variable particlesnum = 	str2num(stringfromlist( 2 ,s3d.paramstring,","))
	variable pd = 				str2num(stringfromlist( 3 ,s3d.paramstring,","))
	variable thickness = 		s3d.thickness
	variable noise = 			str2num(stringfromlist( 4 ,s3d.paramstring,","))
	

	
	make /o /n=(thickness,s3d.num,s3d.num) mat=1,xwave, ywave, zwave
	
	make/B/U /o /n=(thickness,s3d.num,s3d.num,30) exmat= (p <= t) || (q <= t) || (r <= t) || (p >= thickness-t) || (q >= s3d.num-t) || (r >= s3d.num-t) ? 0 : 1
	make/B/U /o /n=(thickness,s3d.num,s3d.num) tempwave
	if(s3d.movie)
		Execute("Spheres3Ddisp(" +num2str(s3d.num)+", \""+getwavesdatafolder(mat,2)+"\")")
		Execute("exportgizmo wave=\"testimage\"  ;Spinoidal3DLayout();Spinoidal3DImage(\""+getdatafolder(1)+"testimage\")")
	endif
	setscale /i x, -thickness/2, thickness/2, mat, exmat,xwave, ywave, zwave
	setscale /i y, -s3d.num/2, s3d.num/2, mat, exmat,xwave, ywave, zwave
	setscale /i z, -s3d.num/2, s3d.num/2, mat, exmat,xwave, ywave, zwave
	xwave = x
	ywave = y
	zwave = z
	redimension /n=(thickness*s3d.num*s3d.num) xwave, ywave, zwave
	variable i,radius, orad, cx, cy, cz, failed, fnum =0, xmn,xmx,ymn,ymx,zmn,zmx, loc
	
	for(i=0;i<particlesnum;i+=1)
		fnum=0
		do
			failed = 0
			radius = abs(gnoise(pd)+s3d.size)
			radius =radius < 1 ? 1 : radius
			if(minsep<1)
				orad = radius*(1+minsep/2)
			else
				orad = radius + minsep/2
			endif
			//duplicate/o /r=()()()(ceil(2*orad)) exmat,tempwave
			redimension /n=(thickness,s3d.num,s3d.num) tempwave
			multithread tempwave[][][] = exmat[p][q][r][ceil(orad)]
//			imagefilter /n=(ceil(2*orad)) /o min3d tempwave
//			multithread tempwave = tempwave<1? 0 : 1 
//			multithread tempwave = (p <= orad) || (q <= orad) || (r <= orad) || (p >= thickness-orad) || (q >= s3d.num-orad) || (r >= s3d.num-orad) ? 0 : tempwave[p][q][r]
  			if(wavemax(tempwave)<1)
				//there are no possible locations for this radius, find another
				failed=1
			else
				// randomly pick a pixel that is good for the center
				redimension /n=(thickness*s3d.num*s3d.num) tempwave
				integrate tempwave /D=intwave
				loc = binarysearch(intwave, enoise(wavemax(intwave)/2)+wavemax(intwave)/2)
				cx = xwave[loc]
				cy = ywave[loc]
				cz = zwave[loc]
			endif
			
			if(failed && (fnum>10))
				print "warning : failed to find unoccupied location - only " +num2str(i) +" particles created"		
				break // get out of the loop
				//we have failed 30 times to create a radius or find a location
			endif
			fnum +=1
		while(failed==1)
		if(failed)
			break // do not add this sphere to the system
					// we are done adding spheres
		endif
		// subtract out this sphere from the matrix  // matrix starts at 1s, within this sphere, multiply this by 0, outside multiply by 1
		multithread mat*= (x-cx)^2 + (y-cy)^2 + (z-cz)^2 < radius^2 ? 0 : 1 
		multithread exmat*= (x-cx)^2 + (y-cy)^2 + (z-cz)^2 <= (orad+t)^2 ? 0 : 1 
		if(s3d.movie)
			execute("ModifyGizmo /n=Spheres3D update=2")
			doupdate
			Execute "exportgizmo wave=\"testimage\"   "
			//TextBox/w=Spinoidal3DLayout/C/N=text0/A=LT/X=0.00/Y=0.00 "\Z32" + time2str2(ttot)
			doupdate
			savepict /p=_PictGallery_ /E=-5 /N=Spinoidal3DLayout /o as "Frame3D"
			addmovieframe /pict=Frame3D
		endif
	endfor
	setdatafolder ::
	imagefilter /n=(interpenetration)/o gauss3d mat
	duplicate /o mat,s3d.density1 // this returns the density matrix of material 1 (the matrix) for alignment etc later on
end

function loadexistingmaterialalignment(s3d)
	// loads existing model into memory for the rest of the calculation
	// no parameter wave is required
	struct ThreeDSystem &s3d
	wave/z s3d.m1 = m1
	wave/z s3d.m2 = m2
	wave/z s3d.m3 = m3
	wave/z s3d.m4 = m4
	wave/z s3d.m5 = m5
end

function/s variables_existing()
	return ""
end
function/s special_existing()
	return ""
end
function model3D_existing(s3d)
	// loads existing model into memory for the rest of the calculation
	// no parameter wave is required
	struct ThreeDSystem &s3d
	wave/z s3d.density1 = density1
	if(!waveexists(s3d.density1))
		print "Density1 wave was not found"
		return -1
	endif
	wave/z s3d.density2 = density2
	if(!waveexists(s3d.density2))
		print "Warning : Density2 wave was not found - If the number of materials is less than 3, this is fine :)"
		return 1
	endif
	wave/z s3d.density3 = density3
	if(!waveexists(s3d.density3))
		print "Warning : Density3 wave was not found - If the number of materials is less than 4, this is fine :)"
		return 1
	endif
	wave/z s3d.density4 = density4
	if(!waveexists(s3d.density4))
		print "Warning : Density4 wave was not found - If the number of materials is less than 5, this is fine :)"
		return 1
	endif
end

function/s variables_existingal()
	return ""
end
function/s special_existingal()
	return "IncludesAlignment"
end
function model3D_existingal(s3d)
	// loads existing model into memory for the rest of the calculation
	// no parameter wave is required
	struct ThreeDSystem &s3d
	wave/z s3d.density1 = density1
	
	wave/z s3d.m1 = m1
	wave/z s3d.m2 = m2
	wave/z s3d.m3 = m3
	wave/z s3d.m4 = m4
	wave/z s3d.m5 = m5
	
	
	
	if(!waveexists(s3d.density1))
		print "Density1 wave was not found"
		return -1
	endif
	wave/z s3d.density2 = density2
	if(!waveexists(s3d.density2))
		print "Warning : Density2 wave was not found - If the number of materials is less than 3, this is fine :)"
		return 1
	endif
	wave/z s3d.density3 = density3
	if(!waveexists(s3d.density3))
		print "Warning : Density3 wave was not found - If the number of materials is less than 4, this is fine :)"
		return 1
	endif
	wave/z s3d.density4 = density4
	if(!waveexists(s3d.density4))
		print "Warning : Density4 wave was not found - If the number of materials is less than 5, this is fine :)"
		return 1
	endif
	
end

function /wave radialintegratew(wave1,minangle,maxangle,outputwavename)
	wave wave1 // input wave to integrate (typically a 2d scattering pattern)
	variable maxangle,minangle // between 0 and 90  ie 80 and 90 will select vertical 20 degrees
							// and likewise, 0 and 10 will collect horizontal 20 degrees
	string outputwavename
	variable centx=0, centy=0
	if(numpnts(wave1)<10)
		make /o /n=0 outputwave
		return outputwave
	endif
	duplicate /o wave1,radialdistance,mask, maskeddata
	//mask = abs(atan((centy-y)/(centx-x))) <maxangle*pi/180 && abs(atan((centy-y)/(centx-x))) >minangle*pi/180 ? 1 : nan
	mask = wave1[p][q]*0==0 && ((atan((centy-y)/(centx-x)) <maxangle*pi/180 && atan((centy-y)/(centx-x)) >minangle*pi/180) || ( (atan((centy-y)/(centx-x)) < (maxangle-180)*pi/180) && (atan((centy-y)/(centx-x)) > (minangle-180)*pi/180) )  ) ? 1 : nan
	radialdistance = sqrt((centx-x)^2 + (centy-y)^2)
	maskeddata *= mask
	radialdistance *= mask
	redimension/n=(numpnts(wave1)) maskeddata
	redimension/n=(numpnts(wave1)) radialdistance
	wavetransform zapnans maskeddata
	wavetransform zapnans radialdistance
	variable mn = wavemin(radialdistance) // the minimum distance for the integration
	variable mx= wavemax(radialdistance) // the maximum distance value for the integration
	variable dmin = 2*min(dimdelta(wave1,0),dimdelta(wave1,1)) // the minimum resolution for the integration
	variable range = round((mx-mn)/(dmin)) // calculate the number of points in the range, with that minimum spacing
	make /d/o $outputwavename , npoints // make the output wave (the reference to which will be returned)
	wave outputwave = $outputwavename
	outputwave = 0
	npoints = 0
	if(mn*dmin*range * 0 != 0)
		make /o /n=0 outputwave
		return outputwave
	endif
	Histogram /B={mn,dmin,range} /c radialdistance, npoints
	Histogram  /B={mn,dmin,range} /c /w=maskeddata radialdistance, outputwave
	outputwave /=npoints
	setscale /p x,mn,dmin, outputwave
	return outputwave
end

function /s time2str(secs)
	variable secs
	string timestr
	variable hours, minutes, seconds
	hours = floor(secs / 3600)
	minutes = floor(secs / 60) - 60*hours
	seconds = 0.01 * floor(secs*100) -60*minutes - 3600*hours
	if(hours>0)
		sprintf timestr "%d:%02d:%02d", hours, minutes, seconds
	elseif(minutes>0)
		sprintf timestr "%d:%02d", minutes, seconds
	else
		sprintf timestr "0:%02d", seconds
	endif
	return timestr
end

window Simulation_Layout() : Layout
	PauseUpdate; Silent 1		// building window...//w=(0,0,100*11,100*8)
	Layout/P=Landscape/C=1/K=1/W=(323,96,1323,846) ScatteringPattern(18,216,396,594)/O=1/F=0
	ModifyLayout mag=1
	PrintSettings /i /w=Simulation_Layout margins={.1,.1,.1,.1}
	Append ScatteringIntensity(400,59,598,599)/O=1/F=0
	Append RatioIntensity(590,56,788,596)/O=1/F=0,BetaProjection(17,17,207,208)/O=1/F=0
	Append DeltaProjection(208,17,398,208)/O=1/F=0
//	Append ParaPerpInt(11.25,339.75,377.25,584.25)/O=1/F=0
//	TextBox/C/N=text0/F=0/A=LB/X=28.63/Y=66.06 "Real part of Polarized Field\rProjection"
//	TextBox/C/N=text1/F=0/A=LB/X=4.77/Y=66.32 "Imaginary part of Polarized Field\rProjection"
//	TextBox/C/N=text2/F=0/A=LB/X=76.74/Y=93.21 "\\Z16Anisotropic Ratio \rvs Q and Energy"
//	TextBox/C/N=text3/F=0/A=LB/X=51.69/Y=93.21 "\\Z16Total Scatter (*q^2) \rvs Q and Energy"
//	TextBox/C/N=text4/F=0/A=LB/X=1.39/Y=97.13 "\\Z16Simulated Scattering Pattern"
//	TextBox/C/N=text5/F=0/A=LB/X=33.80/Y=91.12 "\\F'Arial Black'\\Z24286eV"
//	SetDrawLayer UserFront
//	DrawLine 0,291,1,297
	
	TextBox/C/N=text0/F=0/B=1/A=LT/X=26.5/Y=2/G=(16386,65535,16385) "\Z18Real part of\rPolarized Field"
//	ModifyLayout width(text0)=117,height(text0)=24
	TextBox/C/N=text1/F=0/B=1/A=LT/X=2/Y=2/G=(16386,65535,16385) "\Z18Imaginary part of\rPolarized Field"
//	ModifyLayout width(text1)=139,height(text1)=24
	TextBox/C/N=text2/F=0/B=1/A=LB/X=79.50/Y=93.23 "\\Z16Anisotropic Ratio \rvs Q and Energy"
	TextBox/C/N=text3/F=0/B=1/A=LB/X=54.76/Y=93.40 "\\Z16Total Scatter (*q^2) \rvs Q and Energy"
	TextBox/C/N=text4/F=0/B=1/A=LB/X=3.21/Y=59.36/G=(16386,65535,16385) "\\Z18Simulated Scattering Pattern"
//	ModifyLayout width(text4)=216,height(text4)=25
	TextBox/C/N=text5/F=0/B=1/A=RT/X=51.85/Y=37.15/G=(16386,65535,16385) "\\Z24289 eV"
	//modifylayout trans(text0)=1,trans(text1)=1,trans(text2)=1,trans(text3)=1,trans(text4)=1,trans(text5)=1
	//modifylayout frame=0
	SetDrawLayer UserFront
	DrawLine 0,291,1,297
EndMacro

function ScatteringIntensitydisp()
	dowindow /k ScatteringIntensity
	Display/n=ScatteringIntensity /W=(819,176,991.5,665.75)/K=1 
	AppendImage/w=ScatteringIntensity /T int3DvsEn vs {*,enwavedisp}
	ModifyImage/w=ScatteringIntensity  int3DvsEn ctab= {1e-08,100,BlueHot,0}
	ModifyImage/w=ScatteringIntensity  int3DvsEn log= 1
	ModifyGraph/w=ScatteringIntensity  margin(left)=18,margin(bottom)=14,margin(top)=18,margin(right)=36
	ModifyGraph/w=ScatteringIntensity  log(top)=1
	ModifyGraph/w=ScatteringIntensity  mirror=2
	ModifyGraph/w=ScatteringIntensity  nticks(left)=8
	ModifyGraph/w=ScatteringIntensity  minor=1
	ModifyGraph/w=ScatteringIntensity  fSize(left)=14,fSize(top)=8
	ModifyGraph/w=ScatteringIntensity  standoff=0
	ModifyGraph/w=ScatteringIntensity  tkLblRot(left)=90
	ModifyGraph/w=ScatteringIntensity  btLen=3
	ModifyGraph/w=ScatteringIntensity  tlOffset=-2
	SetAxis/w=ScatteringIntensity /A left
	SetAxis/w=ScatteringIntensity /A  top 
	ColorScale/C/N=text0/F=0/B=1/A=RC/X=0.00/Y=0.00/E=2 image=int3DvsEn, width=10
	ColorScale/C/N=text0 nticks=2, log=1, lblMargin=10
	SetDrawLayer/w=ScatteringIntensity  UserFront
	SetDrawEnv/w=ScatteringIntensity  xcoord= prel,ycoord= left,linethick= 3,linefgc= (65280,0,0),dash= 1
	DrawLine/w=ScatteringIntensity  0,319.999999999984,1,319.999999999984
	//ModifyGraph axRGB=(16385,28398,65535),tlblRGB=(16385,28398,65535),alblRGB=(16385,28398,65535)
End

function RatioIntensitydisp()
	dowindow /k RatioIntensity
	Display/n=RatioIntensity /W=(1002.75,176.75,1175.25,664.25)/K=1 
	AppendImage/w=RatioIntensity /T ratio3DvsEn vs {*,enwavedisp}
	ModifyImage/w=RatioIntensity  ratio3DvsEn ctab= {-1,1,RedWhiteBlue,0}
	ModifyGraph/w=RatioIntensity  margin(left)=18,margin(bottom)=14,margin(top)=18,margin(right)=36
	ModifyGraph/w=RatioIntensity  log(top)=1
	ModifyGraph/w=RatioIntensity  mirror=2
	ModifyGraph/w=RatioIntensity  nticks(left)=8
	ModifyGraph/w=RatioIntensity  minor=1
	ModifyGraph/w=RatioIntensity  fSize(left)=14,fSize(top)=8
	ModifyGraph/w=RatioIntensity  standoff=0
	ModifyGraph/w=RatioIntensity  tkLblRot(left)=90
	ModifyGraph/w=RatioIntensity  btLen=3
	ModifyGraph/w=RatioIntensity  tlOffset=-2
	SetAxis/w=RatioIntensity /A left
	SetAxis/w=RatioIntensity /A top
	ColorScale/C/N=text0/F=0/B=1/A=RC/E=2 image=ratio3DvsEn, width=10, nticks=2
	ColorScale/C/N=text0 lblMargin=5
	SetDrawLayer/w=RatioIntensity  UserFront
	SetDrawEnv/w=RatioIntensity  xcoord= prel,ycoord= left,linethick= 3,linefgc= (65280,0,0),dash= 1
	DrawLine/w=RatioIntensity  0,319.999999999984,1,319.999999999984
End

function ScatteringPatterndisp()
	dowindow /k ScatteringPattern
	Display/n=ScatteringPattern /W=(1151.25,48.5,1361.25,209)/K=1 
	AppendImage/w=ScatteringPattern /T scatter3D
	ModifyImage/w=ScatteringPattern  scatter3D ctab= {1e-05,30000,Grays,0}
	ModifyImage/w=ScatteringPattern  scatter3D ctabAutoscale=1
	ModifyImage/w=ScatteringPattern  scatter3D log= 1
	ModifyGraph/w=ScatteringPattern  margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=58
	ModifyGraph/w=ScatteringPattern  mirror=2
	ModifyGraph/w=ScatteringPattern  nticks=4
	ModifyGraph/w=ScatteringPattern  minor=1
	ModifyGraph/w=ScatteringPattern  fSize=8
	ModifyGraph/w=ScatteringPattern  standoff=0
	ModifyGraph/w=ScatteringPattern  tkLblRot(left)=90
	ModifyGraph/w=ScatteringPattern  btLen=3,mirror=0,nticks=0,minor=0,margin=1
	ModifyGraph/w=ScatteringPattern  tlOffset=-2,height={Plan,1,left,top}
	SetAxis/w=ScatteringPattern /R left 0.5,-0.5
	SetAxis/w=ScatteringPattern  top -0.5,0.5
	ColorScale/w=ScatteringPattern /C/N=text0/F=0/A=MC/X=69.19/Y=2.81 image=scatter3D, log=1
	
	
	appendtograph/w=ScatteringPattern /l=scatl /b=scatb int3Dpara,int3Dperp
	ModifyGraph/w=ScatteringPattern lSize=2
	ModifyGraph/w=ScatteringPattern rgb(int3Dpara)=(65280,16384,16384),rgb(int3Dperp)=(5808,0,5808)
	ModifyGraph/w=ScatteringPattern log(scatl)=1,log(scatb)=1
	ModifyGraph/w=ScatteringPattern mode=7,useNegRGB=1,usePlusRGB=1,toMode=1,plusRGB=(16385,28398,65535)
	ModifyGraph/w=ScatteringPattern hbFill(int3Dpara)=2,rgb=(3,52428,1),standoff=0,lstyle(int3Dpara)=4,lstyle(int3Dperp)=3
	ModifyGraph/w=ScatteringPattern axisEnab(scatl)={0.10,0.5},freePos(scatl)={0.15,kwFraction},axisEnab(scatb)={0.15,1},freePos(scatb)={0.10,kwFraction}
	ModifyGraph axRGB=(16386,65535,16385),tlblRGB=(16386,65535,16385),alblRGB=(16386,65535,16385)
End

function BetaProjectiondisp()
	dowindow/k BetaProjection
	Display/n=BetaProjection /W=(1130.25,433.25,1301.25,594.5)/K=1 
	AppendImage/w=BetaProjection /T pzimagproj
	ModifyImage/w=BetaProjection  pzimagproj ctab= {*,*,Mocha,0}
	ModifyGraph/w=BetaProjection  margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=22
	ModifyGraph/w=BetaProjection  mirror=2
	ModifyGraph/w=BetaProjection  nticks=3
	ModifyGraph/w=BetaProjection  minor=1
	ModifyGraph/w=BetaProjection  fSize=8
	ModifyGraph/w=BetaProjection  standoff=0
	ModifyGraph/w=BetaProjection  tkLblRot(left)=90
	ModifyGraph/w=BetaProjection  btLen=3,mirror=0,nticks=0,minor=0,margin=1
	ModifyGraph/w=BetaProjection  tlOffset=-2,height={Plan,1,left,top}
	setaxis /w=BetaProjection left 0,min(200,dimsize(pzimagproj,1))
	setaxis /w=BetaProjection top 0,min(200,dimsize(pzimagproj,0))
	ColorScale/C/N=text0/F=0/A=RC/X=0.00/Y=0.00 image=pzimagproj, tickLen=1, nticks=2
	ColorScale/C/N=text0 lblMargin=0, lowTrip=0.01, highTrip=100, notation=1
	
	
	
End

function ParaPerpIntdisp()
	dowindow/k ParaPerpInt
	Display/n=ParaPerpInt /W=(798,43.25,1173,213.5)/K=1  int3Dpara,int3Dperp
	ModifyGraph/w=ParaPerpInt lSize=2
	ModifyGraph/w=ParaPerpInt rgb(int3Dpara)=(65280,16384,16384),rgb(int3Dperp)=(5808,0,5808)
	ModifyGraph/w=ParaPerpInt log=1
	ModifyGraph/w=ParaPerpInt mode(int3Dperp)=7,useNegRGB(int3Dperp)=1,usePlusRGB(int3Dperp)=1,toMode(int3Dperp)=1,plusRGB(int3Dperp)=(0,65535,0)
	ModifyGraph/w=ParaPerpInt mode=7,useNegRGB=1,usePlusRGB=1,toMode=1,plusRGB=(0,65535,0)
	ModifyGraph/w=ParaPerpInt hbFill(int3Dpara)=2,rgb(int3Dpara)=(65280,16384,16384,0)
	ModifyGraph/w=ParaPerpInt rgb(int3Dperp)=(5808,0,5808,0)
	//SetAxis/w=ParaPerpInt /A  left
	//Legend/w=ParaPerpInt /C/N=text0/J/F=0/A=LB/X=0.00/Y=0.00 "\\Z14\\s(int3Dpara) Parallel with E-Field\r\\s(int3Dperp) Perpendicular to E-Field"
//	SetDrawLayer/w=ParaPerpInt  UserFront
//	SetDrawEnv/w=ParaPerpInt  xcoord= bottom,ycoord= left,linethick= 3,linefgc= (0,0,65280)
//	DrawLine/w=ParaPerpInt  0.1,100,0.1,1e-05
//	SetDrawEnv/w=ParaPerpInt  xcoord= bottom,ycoord= left,linethick= 3,linefgc= (0,52224,0)
//	DrawLine/w=ParaPerpInt  0.2,100,0.2,1e-05
//	SetDrawEnv/w=ParaPerpInt  linefgc= (0,0,65280),fname= "Arial Black",fsize= 10,textrgb= (0,0,65280)
//	DrawText/w=ParaPerpInt  0.394230769230769,0.180790960451978,"Shell Diameter"
//	SetDrawEnv/w=ParaPerpInt  fname= "Arial Black",fsize= 10,textrgb= (0,52224,0)
//	DrawText/w=ParaPerpInt  0.819711538461539,0.288135593220339,"Core Diameter"
End

function DeltaProjectiondisp()
	dowindow /k deltaprojection	
	Display /n=deltaprojection /W=(891,437.75,1083.75,599)/K=1 
	AppendImage/w=deltaprojection/T pzrealproj
	ModifyImage/w=deltaprojection pzrealproj ctab= {*,*,Mocha,0}
	ModifyGraph/w=deltaprojection margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=22
	ModifyGraph/w=deltaprojection mirror=2
	ModifyGraph/w=deltaprojection nticks=3
	ModifyGraph/w=deltaprojection minor=1
	ModifyGraph/w=deltaprojection fSize=8
	ModifyGraph/w=deltaprojection standoff=0
	ModifyGraph/w=deltaprojection tkLblRot(left)=90
	ModifyGraph/w=deltaprojection btLen=3,mirror=0,nticks=0,minor=0,margin=1
	ModifyGraph/w=deltaprojection tlOffset=-2,height={Plan,1,left,top}
	SetAxis/A/R/w=deltaprojection left
	setaxis /w=deltaprojection left 0,min(200,dimsize(pzrealproj,1))
	setaxis /w=deltaprojection top 0,min(200,dimsize(pzrealproj,0))
	ColorScale/C/N=text0/F=0/A=RC/X=0.00/Y=0.00 image=pzrealproj, tickLen=1, nticks=2
	ColorScale/C/N=text0 lblMargin=0, lowTrip=0.01, highTrip=100, notation=1
End

function PyrealProjectiondisp()
	dowindow /k Pyrealprojection	
	Display /n=Pyrealprojection /W=(891,437.75,1083.75,599)/K=1 
	AppendImage/w=Pyrealprojection /T pyrealproj
	ModifyImage/w=Pyrealprojection pyrealproj ctab= {*,*,Grays256,0}
	ModifyGraph/w=Pyrealprojection margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=22
	ModifyGraph/w=Pyrealprojection mirror=2
	ModifyGraph/w=Pyrealprojection nticks=3
	ModifyGraph/w=Pyrealprojection minor=1
	ModifyGraph/w=Pyrealprojection fSize=8
	ModifyGraph/w=Pyrealprojection standoff=0
	ModifyGraph/w=Pyrealprojection tkLblRot(left)=90
	ModifyGraph/w=Pyrealprojection btLen=3
	ModifyGraph/w=Pyrealprojection tlOffset=-2
	SetAxis/A/R/w=Pyrealprojection left
	setaxis /w=Pyrealprojection left 0,*
	setaxis /w=Pyrealprojection top 0,*
	//ColorScale/w=deltaprojection/C/N=text0/F=0/A=MC/X=41.90/Y=5.59 image=pyrealproj, nticks=2, minor=1
End
function PyimagProjectiondisp()
	dowindow /k Pyimagprojection	
	Display /n=Pyimagprojection /W=(891,437.75,1083.75,599)/K=1 
	AppendImage/w=pyimagprojection/T pyimagproj
	ModifyImage/w=pyimagprojection pyimagproj ctab= {*,*,Grays256,0}
	ModifyGraph/w=pyimagprojection margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=22
	ModifyGraph/w=pyimagprojection mirror=2
	ModifyGraph/w=pyimagprojection nticks=3
	ModifyGraph/w=pyimagprojection minor=1
	ModifyGraph/w=pyimagprojection fSize=8
	ModifyGraph/w=pyimagprojection standoff=0
	ModifyGraph/w=pyimagprojection tkLblRot(left)=90
	ModifyGraph/w=pyimagprojection btLen=3
	ModifyGraph/w=pyimagprojection tlOffset=-2
	SetAxis/A/R/w=pyimagprojection left
	setaxis /w=pyimagprojection left 0,*
	setaxis /w=pyimagprojection top 0,*
	//ColorScale/w=deltaprojection/C/N=text0/F=0/A=MC/X=41.90/Y=5.59 image=pyrealproj, nticks=2, minor=1
End
function pxrealProjectiondisp()
	dowindow /k pxrealprojection	
	Display /n=pxrealprojection /W=(891,437.75,1083.75,599)/K=1 
	AppendImage/w=pxrealprojection/T pxrealproj
	ModifyImage/w=pxrealprojection pxrealproj ctab= {*,*,Grays256,0}
	ModifyGraph/w=pxrealprojection margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=22
	ModifyGraph/w=pxrealprojection mirror=2
	ModifyGraph/w=pxrealprojection nticks=3
	ModifyGraph/w=pxrealprojection minor=1
	ModifyGraph/w=pxrealprojection fSize=8
	ModifyGraph/w=pxrealprojection standoff=0
	ModifyGraph/w=pxrealprojection tkLblRot(left)=90
	ModifyGraph/w=pxrealprojection btLen=3
	ModifyGraph/w=pxrealprojection tlOffset=-2
	SetAxis/A/R/w=pxrealprojection left
	setaxis /w=pxrealprojection left 0,*
	setaxis /w=pxrealprojection top 0,*
	//ColorScale/w=deltaprojection/C/N=text0/F=0/A=MC/X=41.90/Y=5.59 image=pyrealproj, nticks=2, minor=1
End
function pximagProjectiondisp()
	dowindow /k pximagprojection	
	Display /n=pximagprojection /W=(891,437.75,1083.75,599)/K=1 
	AppendImage/w=pximagprojection/T pximagproj
	ModifyImage/w=pximagprojection pximagproj ctab= {*,*,Grays256,0}
	ModifyGraph/w=pximagprojection margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=22
	ModifyGraph/w=pximagprojection mirror=2
	ModifyGraph/w=pximagprojection nticks=3
	ModifyGraph/w=pximagprojection minor=1
	ModifyGraph/w=pximagprojection fSize=8
	ModifyGraph/w=pximagprojection standoff=0
	ModifyGraph/w=pximagprojection tkLblRot(left)=90
	ModifyGraph/w=pximagprojection btLen=3
	ModifyGraph/w=pximagprojection tlOffset=-2
	SetAxis/A/R/w=pximagprojection left
	setaxis /w=pximagprojection left 0,*
	setaxis /w=pximagprojection top 0,*
	//ColorScale/w=deltaprojection/C/N=text0/F=0/A=MC/X=41.90/Y=5.59 image=pyrealproj, nticks=2, minor=1
End


function /wave createalignmentdensity(wavein,sigma,intensity, alignment, movie, timeoffset)
	wave wavein
	variable sigma, intensity // intensity is between 0 and 1
	variable alignment // face on is 1 - the uniaxial vector to be pointing parallel to the interface normal
					// edge on is 0 - the uniaxial vector to be pointing along the interface , perpindicular to normal
	variable movie // wether to add the alignment step to a movie
	variable timeoffset
	variable timer = startmstimer
	duplicate/o wavein, asm, asum, aadd, dx, dy, dz
	//imagefilter /o/n=(3) gauss3d asm
	differentiate /dim=0 asm /D=dx
	differentiate /dim=1 asm /D=dy
	differentiate /dim=2 asm /D=dz
	variable magscale =max(max( max(wavemax(dx), -wavemin(dx)), max(wavemax(dy), -wavemin(dy))) ,max(wavemax(dz), -wavemin(dz)) )
	dx/=magscale
	dy/=magscale
	dz/=magscale
	//create s(x,y,z) vector defining vector of every point in the system
	// sb is the current boundary condition (strength is proportional to magnitude)
	// the first two terms resist a change in s, the last encourages sb to be constant
	// s changes each step, we are trying to minimize E so if a change locally minimizes E for that location, we are good
	make /o/n=(dimsize(wavein,0),dimsize(wavein,1),dimsize(wavein,2),4)  sn
	sn[][][][0] = dx[p][q][r]
	sn[][][][1] = dy[p][q][r]
	sn[][][][2] = dz[p][q][r]
	duplicate/o sn, sb
	//createvectorslice(sn)
	make /o/n=(dimsize(sn,1),dimsize(sn,2)) slice = wavein[10][p][q][0]
	make /o/n=(dimsize(sn,1)*dimsize(sn,2)) yloc, xloc
	make /o/n=(dimsize(sn,1)*dimsize(sn,2),4) arrowsyay
	Alignmentmapdisp()
	if(alignment==1)
		timeoffset += stopmstimer(timer)/1e6
		timer = startmstimer
		//evolvesystem(sn,sb,-500,0,100,-1000,.02,.0001,100, movie, timeoffset) 
		evolvesystem(sn,sb,-10,10,10,-100,.2,.02,200, movie, timeoffset)
	//	timeoffset += stopmstimer(timer)/1e6
	//	timer = startmstimer
		//evolvesystem(sn,sb,-100,100,1000,-1000,.01,.001,100, movie, timeoffset)
	else
		sn=enoise(.1)
		timeoffset += stopmstimer(timer)/1e6
		timer = startmstimer
//		evolvesystem(sn,sb,-1,1,1000,30000,.1,.001,100, movie, timeoffset)
//		timeoffset += stopmstimer(timer)/1e6
//		timer = startmstimer
//		evolvesystem(sn,sb,-100,100,-100,1000,.1,.02,100, movie, timeoffset)
//		timeoffset += stopmstimer(timer)/1e6
//		timer = startmstimer
//		evolvesystem(sn,sb,-100,100,700,1000,.01,.001,100, movie, timeoffset) 
		evolvesystem(sn,sb,-10,10,10,100,.2,.02,200, movie, timeoffset)
	endif
	// enforce alignment strength (intensity) and alignment width (sigma)
	duplicate/o dx, interfaceloc
	multithread interfaceloc = dx^2 + dy^2 + dz^2
	imagefilter /o/n=(sigma) gauss3d interfaceloc
	make/o/n=50 temphist;histogram /B=1 interfaceloc, temphist
	findpeak/Q temphist
	magscale = V_PeakLoc
	multithread interfaceloc = interfaceloc > magscale ? 1 : interfaceloc/magscale
	sn[][][][0,2] *= interfaceloc[p][q][r]
	sn[][][][3] = 1 - sn[p][q][r][0]^2- sn[p][q][r][1]^2- sn[p][q][r][2]^2
	return sn
end	
function evolvesystem(sn,sb,Ec,Ea,Eb,Es,jump,temp,n,movie, timeoffset)
	wave sn, sb
	variable Ec // lower energy when aligned with nearest neighbors	- should be negative
	variable Ea // cost of being different amplitudes of neighbors 	- should be positive
	variable Eb // entropic cost of being non zero  	- should be positive
	variable Es // alignment with barrier 			- should be negative
	variable jump, temp
	variable n, movie, timeoffset
	variable timer = startmstimer, timesofar=timeoffset
	variable numx=dimsize(sn,0), numy=dimsize(sn,1), numz=dimsize(sn,2)
	make/o /n=(dimsize(sn,1)*dimsize(sn,2)) xloc = mod(p,dimsize(sn,1)),yloc = floor(p/dimsize(sn,1))
	make/o /n=(dimsize(sn,1)*dimsize(sn,2)) xcomp,ycomp,zcomp
	make/o /n=(dimsize(sn,1),dimsize(sn,2)) slice
	make/o /n=(dimsize(sn,1)*dimsize(sn,2),2) arrowsyay
	
	
	make /o/n=(dimsize(sn,0),dimsize(sn,1),dimsize(sn,2))  En
	newdatafolder /o/s alignment
	make /o/n=(dimsize(sn,0),dimsize(sn,1),dimsize(sn,2))  Ensbmag, snmag,sbmag, stestmag, Etest, magsn, magstest, magsn1, magsn2, magsn3, magsn4, magsn5, magsn6,magsn7, magsn8, magsn9, magsn10, magsn11, magsn12,magsn13, magsn14, magsn15, magsn16,magsn17, magsn18, tempw
	duplicate/o sn,stest 
	
	sbmag = sqrt(sb[p][q][r][0]^2+sb[p][q][r][1]^2+sb[p][q][r][2]^2)
	//sn = enoise(sqrt(3))
	variable i
	for(i=0;i<n;i+=1)
		multithread snmag = sqrt(sn[p][q][r][0]^2+sn[p][q][r][1]^2+sn[p][q][r][2]^2)
		multithread sn /= snmag[p][q][r]>1 ? snmag[p][q][r] : 1 // make sure magnitute of n can't get bigger than 1
		multithread snmag = sqrt(sn[p][q][r][0]^2+sn[p][q][r][1]^2+sn[p][q][r][2]^2)
		
		multithread stest = max(-1,min(1,sn+gnoise(jump)))
		
		multithread stestmag = sqrt(stest[p][q][r][0]^2+stest[p][q][r][1]^2+stest[p][q][r][2]^2)
		multithread stest /= stestmag[p][q][r]>1 ? stestmag[p][q][r] : 1 // make sure magnitute of test can't get bigger than 1
		multithread stestmag = sqrt(stest[p][q][r][0]^2+stest[p][q][r][1]^2+stest[p][q][r][2]^2)
		
		multithread magsn1[][][] = snmag[mod(p-1+numx,numx)][q][r]
		multithread magsn2[][][] = snmag[mod(p+1+numx,numx)][q][r]
		multithread magsn3[][][] = snmag[p][mod(q-1+numy,numy)][r]
		multithread magsn4[][][] = snmag[p][mod(q+1+numy,numy)][r]
		multithread magsn5[][][] = snmag[p][q][mod(r-1+numz,numz)]
		multithread magsn6[][][] = snmag[p][q][mod(r+1+numz,numz)]
		
		multithread magsn7[][][] = snmag[mod(p+1+numx,numx)][mod(q+1+numy,numy)][r]
		multithread magsn8[][][] = snmag[mod(p+1+numx,numx)][mod(q-1+numy,numy)][r]
		multithread magsn9[][][] = snmag[mod(p-1+numx,numx)][mod(q+1+numy,numy)][r]
		multithread magsn10[][][] = snmag[mod(p-1+numx,numx)][mod(q-1+numy,numy)][r]
		multithread magsn11[][][] = snmag[mod(p+1+numx,numx)][q][mod(r+1+numz,numz)]
		multithread magsn12[][][] = snmag[mod(p+1+numx,numx)][q][mod(r-1+numz,numz)]
		multithread magsn13[][][] = snmag[mod(p-1+numx,numx)][q][mod(r+1+numz,numz)]
		multithread magsn14[][][] = snmag[mod(p-1+numx,numx)][q][mod(r-1+numz,numz)]
		multithread magsn15[][][] = snmag[p][mod(q+1+numy,numy)][mod(r+1+numz,numz)]
		multithread magsn16[][][] = snmag[p][mod(q+1+numy,numy)][mod(r-1+numz,numz)]
		multithread magsn17[][][] = snmag[p][mod(q-1+numy,numy)][mod(r+1+numz,numz)]
		multithread magsn18[][][] = snmag[p][mod(q-1+numy,numy)][mod(r-1+numz,numz)]

		// make sure the deominators aren't zero, just very small
		multithread magsn = magsn<.001 ?.001 : magsn
		multithread magsn1 = magsn1<.001 ?.001 : magsn1
		multithread magsn2 = magsn2<.001 ?.001 : magsn2
		multithread magsn3 = magsn3<.001 ?.001 : magsn3
		multithread magsn4 = magsn4<.001 ?.001 : magsn4
		multithread magsn5 = magsn5<.001 ?.001 : magsn5
		multithread magsn6 = magsn6<.001 ?.001 : magsn6
		multithread magstest = magstest<.001 ?.001 : magstest
		
		multithread Etest = ((stest[p][q][r][0] * sn[mod(p-1+numx,numx)][q][r][0]+stest[p][q][r][1] * sn[mod(p-1+numx,numx)][q][r][1]+stest[p][q][r][2] * sn[mod(p-1+numx,numx)][q][r][2])   )^2//  /(magsn2[p][q][r]*magstest[p][q][r]))^2
		multithread Etest += ((stest[p][q][r][0] * sn[mod(p+1+numx,numx)][q][r][0]+stest[p][q][r][1] * sn[mod(p+1+numx,numx)][q][r][1]+stest[p][q][r][2] * sn[mod(p+1+numx,numx)][q][r][2])   )^2//  /(magsn2[p][q][r]*magstest[p][q][r]))^2
		multithread Etest += ((stest[p][q][r][0] * sn[p][mod(q-1+numy,numy)][r][0]+stest[p][q][r][1] * sn[p][mod(q-1+numy,numy)][r][1]+stest[p][q][r][2] * sn[p][mod(q-1+numy,numy)][r][2])   )^2//  /(magsn3[p][q][r]*magstest[p][q][r]))^2
		multithread Etest += ((stest[p][q][r][0] * sn[p][mod(q+1+numy,numy)][r][0]+stest[p][q][r][1] * sn[p][mod(q+1+numy,numy)][r][1]+stest[p][q][r][2] * sn[p][mod(q+1+numy,numy)][r][2])   )^2//  /(magsn4[p][q][r]*magstest[p][q][r]))^2
		multithread Etest += ((stest[p][q][r][0] * sn[p][q][mod(r-1+numz,numz)][0]+stest[p][q][r][1] * sn[p][q][mod(r-1+numz,numz)][1]+stest[p][q][r][2] * sn[p][q][mod(r-1+numz,numz)][2])  )^2//  /(magsn5[p][q][r]*magstest[p][q][r]))^2
		multithread Etest += ((stest[p][q][r][0] * sn[p][q][mod(r+1+numz,numz)][0]+stest[p][q][r][1] * sn[p][q][mod(r+1+numz,numz)][1]+stest[p][q][r][2] * sn[p][q][mod(r+1+numz,numz)][2])  )^2//  /(magsn6[p][q][r]*magstest[p][q][r]))^2
		multithread Etest *= Ec // lower energy if aligned with neighbors// between 0 and 6, 
		multithread Etest += Ea * 2*((magsn1[p][q][r] + magsn2[p][q][r] + magsn3[p][q][r] + magsn4[p][q][r] + magsn5[p][q][r] + magsn6[p][q][r])/6 - magstest[p][q][r])^2 // magnitudes, between 0 and 6,  0 - magnitude is average of neighboring
		multithread Etest += Ea * ((magsn7[p][q][r] + magsn8[p][q][r] + magsn9[p][q][r] + magsn10[p][q][r] + magsn11[p][q][r] + magsn12[p][q][r] +magsn13[p][q][r] + magsn14[p][q][r] + magsn15[p][q][r] + magsn16[p][q][r] + magsn17[p][q][r] + magsn18[p][q][r])/12 - magstest[p][q][r])^2 // magnitudes, between 0 and 6,  0 - magnitude is average of neighboring
		multithread Etest += Eb * ( stest[p][q][r][0]^2 + stest[p][q][r][1]^2+ stest[p][q][r][2]^2 ) // entropic cost of alignment
		multithread Etest += Es * abs(stest[p][q][r][0] *sb[p][q][r][0] + stest[p][q][r][1] * sb[p][q][r][1] + stest[p][q][r][2] * sb[p][q][r][2] )
		multithread Etest -= Es * sqrt((stest[p][q][r][1] *sb[p][q][r][0] - stest[p][q][r][0] * sb[p][q][r][1])^2  + (stest[p][q][r][2] *sb[p][q][r][0] - stest[p][q][r][0] * sb[p][q][r][2])^2  + (stest[p][q][r][1] *sb[p][q][r][2] - stest[p][q][r][2] * sb[p][q][r][1])^2 ) // lower energy, the better aligned with boundary

		
		multithread En = ((sn[p][q][r][0] * sn[mod(p-1+numx,numx)][q][r][0]+sn[p][q][r][1] * sn[mod(p-1+numx,numx)][q][r][1]+sn[p][q][r][2] * sn[mod(p-1+numx,numx)][q][r][2])   )^2//  /(magsn2[p][q][r]*magsn[p][q][r]))^2
		multithread En += ((sn[p][q][r][0] * sn[mod(p+1+numx,numx)][q][r][0]+sn[p][q][r][1] * sn[mod(p+1+numx,numx)][q][r][1]+sn[p][q][r][2] * sn[mod(p+1+numx,numx)][q][r][2])   )^2//  /(magsn2[p][q][r]*magsn[p][q][r]))^2
		multithread En += ((sn[p][q][r][0] * sn[p][mod(q-1+numy,numy)][r][0]+sn[p][q][r][1] * sn[p][mod(q-1+numy,numy)][r][1]+sn[p][q][r][2] * sn[p][mod(q-1+numy,numy)][r][2])   )^2//  /(magsn3[p][q][r]*magsn[p][q][r]))^2
		multithread En += ((sn[p][q][r][0] * sn[p][mod(q+1+numy,numy)][r][0]+sn[p][q][r][1] * sn[p][mod(q+1+numy,numy)][r][1]+sn[p][q][r][2] * sn[p][mod(q+1+numy,numy)][r][2])   )^2//  /(magsn4[p][q][r]*magsn[p][q][r]))^2
		multithread En += ((sn[p][q][r][0] * sn[p][q][mod(r-1+numz,numz)][0]+sn[p][q][r][1] * sn[p][q][mod(r-1+numz,numz)][1]+sn[p][q][r][2] * sn[p][q][mod(r-1+numz,numz)][2])  )^2//  /(magsn5[p][q][r]*magsn[p][q][r]))^2
		multithread En += ((sn[p][q][r][0] * sn[p][q][mod(r+1+numz,numz)][0]+sn[p][q][r][1] * sn[p][q][mod(r+1+numz,numz)][1]+sn[p][q][r][2] * sn[p][q][mod(r+1+numz,numz)][2])  )^2//  /(magsn6[p][q][r]*magsn[p][q][r]))^2
		multithread En *= Ec // lower energy if aligned with neighbors// between 0 and 6, 
		multithread En += Ea * 2*((magsn1[p][q][r] + magsn2[p][q][r] + magsn3[p][q][r] + magsn4[p][q][r] + magsn5[p][q][r] + magsn6[p][q][r])/6 - magsn[p][q][r])^2 // magnitudes, between 0 and 6,  0 - magnitude is average of neighboring
		multithread En += Ea * ((magsn7[p][q][r] + magsn8[p][q][r] + magsn9[p][q][r] + magsn10[p][q][r] + magsn11[p][q][r] + magsn12[p][q][r] +magsn13[p][q][r] + magsn14[p][q][r] + magsn15[p][q][r] + magsn16[p][q][r] + magsn17[p][q][r] + magsn18[p][q][r])/12 - magsn[p][q][r])^2 // magnitudes, between 0 and 6,  0 - magnitude is average of neighboring
		multithread En += Eb * ( sn[p][q][r][0]^2 + sn[p][q][r][1]^2+ sn[p][q][r][2]^2 ) // entropic cost of alignment
		multithread En += Es * abs(sn[p][q][r][0] *sb[p][q][r][0] + sn[p][q][r][1] * sb[p][q][r][1] + sn[p][q][r][2] * sb[p][q][r][2] )
		multithread En -= Es * sqrt((sn[p][q][r][1] *sb[p][q][r][0] - sn[p][q][r][0] * sb[p][q][r][1])^2  + (sn[p][q][r][2] *sb[p][q][r][0] - sn[p][q][r][0] * sb[p][q][r][2])^2  + (sn[p][q][r][1] *sb[p][q][r][2] - sn[p][q][r][2] * sb[p][q][r][1])^2 ) // lower energy, the better aligned with boundary
		
		timesofar += stopmstimer(timer)/1e6
		timer = startmstimer
		print "Alignment Time : "+Time2str(timesofar) +" - Alignment Step " + num2str(i) + " of "+num2str(n)+ "  -  Total Alignment Energy : " + num2str(mean(En))
//metropolis algorithm 
		multithread tempw[][][] = Exp(-(Etest[p][q][r]- En[p][q][r])/temp) > enoise(.5)+.5 ? 1 : 0 // one test for all three components at each location ie, don't just keep the x component and throw the y and z out
		multithread sn[][][][] = tempw[p][q][r] ? stest[p][q][r][t] : sn[p][q][r][t]
		//multithread sn[][][][1] *= sn[p][q][r][0] < 0 ? -1 : 1
		//multithread sn[][][][2] *= sn[p][q][r][0] < 0 ? -1 : 1
		//multithread sn[][][][0] *= sn[p][q][r][0] < 0 ? -1 : 1
		
		multithread xcomp = sn[10][xloc][yloc][0]
		multithread ycomp = sn[10][xloc][yloc][1]
		multithread zcomp = sn[10][xloc][yloc][2]
		multithread arrowsyay[][1] = atan(zcomp[p]/ycomp[p])
		multithread arrowsyay[][0] = sqrt(zcomp[p]^2 + ycomp[p]^2) * 20
		doupdate
		if(movie)
			addmovieframe
		endif
	endfor
	timesofar += stopmstimer(timer)/1e6
	killwaves /a/z
	setdatafolder ::
	killdatafolder/z alignment
end

function Alignmentmapdisp()
	dowindow /k alignmentmap
	wave yloc, xloc, slice
	Display /k=1/n=alignmentmap /W=(215.25,45.5,1011.75,731)/K=1  yloc vs xloc as "Slice of Alignment"
	AppendImage/w=alignmentmap slice
	ModifyImage/w=alignmentmap slice ctab= {-1,1,BlueRedGreen,0}
	ModifyGraph/w=alignmentmap mode=3
	ModifyGraph/w=alignmentmap rgb=(65535,65535,65535)
	ModifyGraph/w=alignmentmap msize=0.01
	ModifyGraph/w=alignmentmap mrkThick=0.1
	ModifyGraph/w=alignmentmap arrowMarker(yloc)={arrowsyay,1,0,0,1}
	ModifyGraph/w=alignmentmap mirror=0
	ModifyGraph/w=alignmentmap axOffset(left)=-4.46154,axOffset(bottom)=-0.962963
	setaxis /w=alignmentmap left 0,min(100,dimsize(slice,0))
	setaxis /w=alignmentmap bottom 0,min(100,dimsize(slice,1))
EndMacro

Function Model3DTabProc(tca) : TabControl
	STRUCT WMTabControlAction &tca

	switch( tca.eventCode )
		case 2: // mouse up
			Variable tab = tca.tab
			nvar usealignment = root:Packages:ScatterSim3D:useprecalcalignment
			nvar usemorphalignment =  root:Packages:ScatterSim3D:usemorphalignment
			svar morphology = root:Packages:ScatterSim3D:morphology
			groupbox MorphOptions,disable= (tab!=0)
			popupmenu MorphologyPop,disable= (tab!=0)
			setvariable setres,disable= (tab!=0)
			setvariable setsize,disable= (tab!=0)
			setvariable setsize1,disable= (tab!=0)
			setvariable SetThickness,disable= (tab!=0) 
			button CalcMorphBut,disable= (tab!=0)
			popupmenu MaterialPop,disable= (tab!=1 ||  usealignment)
			button AddMaterial,disable= (tab!=1 ||  usealignment)
			button CalcAlignBut,disable= (tab!=1 ||  usealignment)
			button removematerial,disable= (tab!=1 ||  usealignment)
			listbox materials ,disable= (tab!=1 ||  usealignment)
			groupbox MaterialsListText,disable= (tab!=1 ||  usealignment)
			checkbox usealign,disable= (tab!=1 || cmpstr(morphology,"existing"))
			setvariable setmaxen,disable= (tab!=2)
			setvariable setensteps,disable= (tab!=2)
			setvariable setminen,disable= (tab!=2)
			setvariable setefield,disable= (tab==2)+1
			CheckBox SaveMovieChk,disable= (tab!=2)
			CheckBox Save2DDataEn,disable= (tab!=2) 
			CheckBox rotatesystemck,disable= (tab!=2) 
			CheckBox SaveQEn,disable= (tab!=2)
			CheckBox SaveParaPerp,disable= (tab!=2)
			CheckBox SaveAnisotropy,disable= (tab!=2)
			svar controllist = root:Packages:ScatterSim3D:controllist
			variable i
			string controlname, controltype
			for(i=0;i<itemsinlist(controllist);i+=1)
				controlname = stringfromlist(0,stringfromlist(i,controllist),",")
				controltype = stringfromlist(2,stringfromlist(i,controllist),",")
				if(stringmatch("checkbox",controltype))
					checkbox $controlname win=ScatteringSimPanel, disable=(tab!=0)
				elseif(stringmatch("SetVariable",controltype))
					setvariable $controlname win=ScatteringSimPanel, disable=(tab!=0)
				elseif(stringmatch("String",controltype))
					setvariable $controlname win=ScatteringSimPanel, disable=(tab!=0)
				endif
			endfor
			svar material = root:Packages:ScatterSim3D:material
			if(OCsexist(material)!=2)
				variable val = (tab!=1 || usealignment || usemorphalignment) ? 1 : 2
				PopupMenu AlignmentPop disable=val, mode=1
				SetVariable AlignmentSize disable=val
				SetVariable AlignmentStrength disable=val
			else
				PopupMenu AlignmentPop disable=(tab!=1 ||  usealignment || usemorphalignment)
				SetVariable AlignmentSize disable=(tab!=1 ||  usealignment|| usemorphalignment)
				SetVariable AlignmentStrength disable=(tab!=1 ||  usealignment|| usemorphalignment)
			endif
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

function Model3DPanel()
	string foldersave = getdatafolder(1)
	if(!datafolderexists("root:opticalconstants"))
		loadoc()
	endif
	Newdatafolder /O/S root:packages
	newdatafolder /o/s ScatterSim3D
	PauseUpdate; Silent 1		// building window...
	dowindow /k ScatteringSimPanel
	NewPanel/n=ScatteringSimPanel /k=1/W=(1352,181,1829,688) as "Scattering Simulations"
	modifypanel /w=ScatteringSimPanel fixedsize=1
	SetDrawLayer UserBack
	string funcs = functionlist("model3d_*",";","")
	variable num = itemsinlist(funcs)
	variable i
	string funcname
	svar/z extralist
	if(svar_exists(extralist))
		svar funcnames, morphology, material, alignment, SimName, efield,controllist,extralist
		nvar sizescale, thickness, resolution, startxrayenergy, endxrayenergy, numensteps, movie, voxelsize, usealignment, rotatesys
		nvar alignmentsize, alignmentstrength, useprecalcalignment,SaveScattering1D, SaveParaPerp1D, SaveAnisotropic1D, SaveScattering2D
	else
		string/g funcnames="", morphology="existing", material="pcbm", alignment="None", SimName = "test", efield = " ( 0 , 1 )"
		string /g controllist="", extralist=""
		variable /g sizescale=5, thickness=16, resolution=128, startxrayenergy=260, endxrayenergy=320, numensteps=600, movie=0, voxelsize=5
		variable /g alignmentsize=10, alignmentstrength=1, useprecalcalignment=0, usemorphalignment=0
		variable /g  SaveScattering1D=1, SaveParaPerp1D=1, SaveAnisotropic1D=1, SaveScattering2D=1, RotateSys=1
	endif
	make /n=(0,2)/o/t ListOfMaterials, ListColumnNames = {"\JCMaterial","\JCAlignment"}
	make /n=(0,2)/o SelectedMaterials
	funcnames = ""
	for(i=0;i<num;i+=1)
		splitstring/e="^[^_]*_(.*)$" stringfromlist(i,funcs), funcname
		funcnames = addlistitem(funcname,funcnames,";") 
	endfor
	string funclistname = "root:packages:ScatterSim3D:funcnames"
	
	PopupMenu MorphologyPop,pos={25,57},size={131,21},proc=PopMorph,title="Morphology : "
	PopupMenu MorphologyPop,mode=2,popvalue="existing",value= #"root:packages:ScatterSim3D:funcnames"
	ListBox Materials,pos={219,55},size={225,175},disable=1,frame=4
	ListBox Materials,listWave=root:Packages:ScatterSim3D:ListOfMaterials
	ListBox Materials,selWave=root:Packages:ScatterSim3D:SelectedMaterials
	ListBox Materials,titleWave=root:Packages:ScatterSim3D:ListColumnNames,mode= 4
	ListBox Materials,editStyle= 1,widths={50,80},userColumnResize= 1
	PopupMenu AlignmentPop,pos={62,80},size={113,21},disable=1,proc=AlignmentPop,title="Alignment : "
	PopupMenu AlignmentPop,mode=3,popvalue="None",value= #"\"None;Face-On;Edge-On\""
	Button AddMaterial,pos={40,177},size={134,37},disable=1,proc=Model3DProc,title="Add Material ->"
	PopupMenu MaterialPop,pos={30,48},size={165,21},bodyWidth=116,disable=1,proc=MaterialPop,title="Material : "
	PopupMenu MaterialPop,mode=12,popvalue="PCBM",value= #"listocs()"
	SetVariable SetMinEn,pos={23,56},size={180,16},disable=1,title="Start X-ray Energy [eV]"
	SetVariable SetMinEn,limits={10,30000,0.05},value= root:Packages:ScatterSim3D:startxrayenergy,live= 1
	SetVariable SetMaxEn,pos={24,77},size={180,16},disable=1,title="End X-ray Energy [eV]"
	SetVariable SetMaxEn,limits={10,30000,0.05},value= root:Packages:ScatterSim3D:endxrayenergy,live= 1
	SetVariable SetEnSteps,pos={24,97},size={180,16},disable=1,title="Number of Steps"
	SetVariable SetEnSteps,limits={1,30000,0.05},value= root:Packages:ScatterSim3D:numensteps,live= 1
	GroupBox MorphOptions,pos={20,87},size={431,314},title="Variables for Morphology"
	Button RemoveMaterial,pos={40,216},size={134,39},disable=1,proc=RemoveMaterialButton,title="Remove Material"
	TabControl tab0,pos={8,7},size={458,443},proc=Model3dTabProc,tabLabel(0)="Morphology"
	TabControl tab0,tabLabel(1)="Materials",tabLabel(2)="Scattering",value= 0
	GroupBox MaterialsListText,pos={206,37},size={249,206},disable=1,title="List of Materials"
	SetVariable SetRes,pos={20,31},size={275,16},title="Resolution (Film Width and Length) [Voxels]"
	SetVariable SetRes,limits={10,5000,1},value= root:Packages:ScatterSim3D:resolution,live= 1
	SetVariable SetSize,pos={180,60},size={145,16},title="Length Scale [Voxels]"
	SetVariable SetSize,limits={1,5000,1},value= root:packages:ScatterSim3D:sizescale,live= 1
	Button CalcMorphBut,pos={262,409},size={192,31},title="Calculate Morphology Now"
	Button CalcAlignBut,pos={206,250},size={243,38},disable=1,title="Calculate Material Alignment Now"
	CheckBox SaveMovieChk,pos={234,54},size={75,14},disable=1,title="Save Movie"
	CheckBox SaveMovieChk,variable= root:Packages:ScatterSim3D:movie
	CheckBox SaveQEn,pos={234,72},size={161,14},disable=1,title="Save 1D Scattering vs Energy"
	CheckBox SaveQEn,variable= root:Packages:ScatterSim3D:SaveScattering1D
	CheckBox Save2DDataEn,pos={234,121},size={187,14},disable=1,title="Save 2D Scattering Data vs Energy"
	CheckBox Save2DDataEn,variable= root:Packages:ScatterSim3D:SaveScattering2D
	CheckBox rotatesystemck,pos={234,137},size={187,14},disable=1,title="Integrate System rotated (best with periodic boundaries)"
	CheckBox rotatesystemck,variable= root:Packages:ScatterSim3D:RotateSys
	CheckBox SaveParaPerp,pos={234,88},size={206,14},disable=1,title="Save Parallel and Perpindicular 1D data"
	CheckBox SaveParaPerp,variable= root:Packages:ScatterSim3D:SaveParaPerp1D
	CheckBox SaveAnisotropy,pos={234,104},size={177,14},disable=1,title="Save Anisotropic Ratio Vs Energy"
	CheckBox SaveAnisotropy,variable= root:Packages:ScatterSim3D:SaveAnisotropic1D
	SetVariable SetThickness,pos={300,33},size={146,16},title="Thickness (Pixels)"
	SetVariable SetThickness,limits={1,5000,1},value= root:Packages:ScatterSim3D:thickness,live= 1
	SetVariable AlignmentSize,pos={43,108},size={153,30},disable=1,title="\\JCApproximate Width\r of Alignment [pixels]"
	SetVariable AlignmentSize,limits={0,1000,0.01},value= root:Packages:ScatterSim3D:AlignmentSize
	SetVariable AlignmentStrength,pos={21,149},size={173,16},disable=1,title="Strength of Alignment [0-1]"
	SetVariable AlignmentStrength,limits={0,1,0.01},value= root:Packages:ScatterSim3D:AlignmentStrength
	Button SimulateSystem,pos={253,455},size={210,47},proc=SimulateSystem,title="\\Z20Simulate System!"
	SetVariable NameOfSim,pos={8,472},size={240,16},title="Name for Simulation :"
	SetVariable NameOfSim,value= root:Packages:ScatterSim3D:SimName
	SetVariable setefield,pos={30,163},size={300,16},disable=1,proc=SetEField,title="E Field ( y component, z component )"
	SetVariable setefield,value= root:Packages:ScatterSim3D:efield
	CheckBox usealign,pos={266,410},size={158,14},title="Use Pre Calculated Alignment"
	CheckBox usealign,variable=useprecalcalignment,disable=1,proc=UserPreCalcChk
	SetVariable SetSize1,pos={335,61},size={115,16},title="Voxel Size [nm]"
	SetVariable SetSize1,limits={1,5000,1},value= root:packages:ScatterSim3D:voxelsize,live= 1
	setdatafolder foldersave
End
function /s time2str2(secs)
	variable secs
	string timestr
	variable hours, minutes, seconds
	hours = floor(secs / 3600)
	minutes = floor(secs / 60) - 60*hours
	seconds = secs -60*minutes - 3600*hours
	if(hours>0)
		sprintf timestr "%d:%02d:%02d", hours, minutes, seconds
	elseif(minutes>0)
		sprintf timestr "%d:%02d", minutes, seconds
	else
		if(seconds<1e-9)
			sprintf timestr "%03.2f ps", (seconds*1e12)
		elseif(seconds<1e-6)
			sprintf timestr "%03.2f ns", (seconds*1e9)
		elseif(seconds<1e-3)
			sprintf timestr "%03.2f µs", (seconds*1e6)
		elseif(seconds<1)
			sprintf timestr "%03.2f ms", (seconds*1e3)
		else
			sprintf timestr "0:%02d s", seconds
		endif
	endif
	return timestr
end
function /s variables_spinoidal()
//"Minimum Seperation,SetVariable,.1;Number of Particles (Max),SetVariable,500;PolyDispursity (sigma of radiuses),setvariable,5;"
	return "Number of Steps,SetVariable,100;Size Parameter (.001-.01),SetVariable,.005;Exponential component of Time Steps,SetVariable,.7;Linear component of Time Steps,SetVariable,5"
end
function model3d_spinoidal(s3d)
	struct ThreeDSystem &s3d
	variable ndim = s3d.num
	variable thickness = s3d.thickness
	variable num = str2num(stringfromlist(0,s3d.paramstring,","))
	variable eps = str2num(stringfromlist(1,s3d.paramstring,","))
	variable texp = str2num(stringfromlist(2,s3d.paramstring,","))
	variable t0 = str2num(stringfromlist(3,s3d.paramstring,","))
//	setdatafolder root:packages:ScatterSim3D
	make /o/n=(s3d.thickness,s3d.num,s3d.num) /o s3d.density1
	newdatafolder /o/s modelCreation
	spinoidalLB3Dn(ndim, num,eps, thickness,texp,t0,movie = s3d.movie)
	wave realu
	s3d.density1 = realu[q][r][p] // density's smaller dimension is p, whereas realu's smallest dimension is layers
	killwaves /A/Z
	setdatafolder ::
	s3d.density1 = s3d.density1 > 0 ? 1 : 0
end

function spinoidalLB3Dn(ndim, num,epsilon, thickness,texp,t0,[movie])
	variable ndim, num, thickness, texp, t0, movie,epsilon
	variable M=ndim,N=ndim
	variable delx = 1/(M-1)
	variable delx2 = delx^2
	variable visualupdate =10
	variable typeupdate = 10
	//variable delt=.00005
	//variable epsilon = .004
	variable eps2 = epsilon^2
	variable a = 2
	movie = paramisdefault(movie) ? 1 : movie
	variable comptime = 0
	variable timer = startmstimer
	make /o/C/n=(M,N, thickness) Leig, seig, Cheig,  hatU, hatrhs, U , fU
	make /o /n=(M,N, thickness) realu, imagu
	make /o/n=(num) delt = .000005*exp(p^texp/t0)
	
	variable lam1 = delt[0]/delx2
	variable lam2 = eps2*lam1/delx2
	
	multithread Leig = (2*Exp(sqrt(-1)*pi*(p)/((M-1)/2))-2) + (2*Exp(sqrt(-1)*pi*(q)/((N-1)/2))-2)+ (2*Exp(sqrt(-1)*pi*(r)/((thickness-1)/2))-2)
	multithread imagu = real(leig)
	multithread leig = cmplx(imagu,0)
	
	multithread cheig = 1-a*lam1*leig[p][q][r] + lam2*leig[p][q][r]^2
	multithread seig = lam1*Leig

	multithread realu= enoise(.05)
	if(movie)
		variable sizeofdisplay
		if(ndim>128)
			sizeofdisplay = 128
			make/o /n=(128,128,thickness) realudisp
		else
			sizeofdisplay = ndim
			duplicate/o realu, realudisp
		endif
		Execute("Spin3Ddisp(" +num2str(sizeofdisplay)+", \""+getwavesdatafolder(realudisp,2)+"\")")
		Execute("exportgizmo wave=\"testimage\"   ;Spinoidal3DLayout_2();Spinoidal3DImage(\""+getdatafolder(1)+"testimage\")")
		makespinoidalslices(realu)
	endif
	
	multithread U = cmplx(realu,0)
	//newimage realu
	makespinoidalslices(realu)
	fft /dest=hatu U
	variable i, ttot=0
	for(i=0;i<num;i+=1)
		comptime +=stopmstimer(timer)/10^6
		timer = startmstimer
		ttot +=delt[i]
		print "Frame " + num2str(i) + "  -  " + time2str2(ttot) + " in simulation  -  " +time2str2(comptime) + " computation time elapsed"
		
		lam1 = delt[i]/delx2
		lam2 = eps2*lam1/delx2
		multithread cheig = 1-a*lam1*leig[p][q] + lam2*leig[p][q]^2
		multithread seig = lam1*Leig

		multithread fU = U[p][q][r]^3 - ((1+a)*U[p][q][r]);
//		multithread fu[0][][] *= fu*(-1)^(p+q+r) < 0 ? -1 : 1  // these two lines are my attempt to add wetting at the surfaces, although it doesn't seem to work that well.
//		multithread fu[inf][][] *= fu*(-1)^(p+q+r) > 0 ? -1 : 1
		
		fft /dest=tempfu fU
		
		multithread hatrhs = hatu + seig[p][q][r] * tempfu[p][q][r]
		multithread hatu = hatrhs[p][q][r] / cheig[p][q][r]
		ifft /c/dest=u hatU
		multithread u = cmplx( real(u[p][q][r]),0)
		multithread realu =real(u[p][q][r])* (-1)^(p+q+r)
		if(movie)
			multithread realudisp = realu[p][q][r]
			execute("ModifyGizmo /n=Spheres3D update=2")
			doupdate
			Execute "exportgizmo wave=\"testimage\"   "
			//TextBox/w=Spinoidal3DLayout/C/N=text0/A=LT/X=0.00/Y=0.00 "\Z32" + time2str2(ttot)
			dowindow /F Spinoidal3DLayout_2
			doupdate
			savepict /p=_PictGallery_ /E=-5 /N=Spinoidal3DLayout_2 /o as "Frame3D"
			addmovieframe /pict=Frame3D
		endif
		doupdate
	endfor
	if(movie)
		dowindow /k Spheres3D // the generaic 3D model
		dowindow /k Spinoidal3Dlayout_2
		dowindow /k Spinoidal3Dimage
	endif
	dowindow /k unscaledslice
	dowindow /k scaledslice
end
	
window Spheres3Ddisp(size, densitywaveloc) : GizmoPlot
	variable size
	string densitywaveloc
	PauseUpdate; Silent 1	// Building Gizmo 6 window...
	dowindow /k Spheres3D
	// Do nothing if the Gizmo XOP is not available.
	if(exists("NewGizmo")!=4)
		DoAlert 0, "Gizmo XOP must be installed"
		return 0
	endif

	NewGizmo/N=Spheres3D/T="3D Spheres" /W=(30,30,1310,1054) /k=1
	ModifyGizmo startRecMacro
	ModifyGizmo scalingMode=8
	ModifyGizmo setOuterBox={-size/2,size/2,-size/2,size/2,-size/2,size/2}
	ModifyGizmo scalingOption=60
	AppendToGizmo isoSurface=$densitywaveloc,name=isoSurface0
	ModifyGizmo ModifyObject=isoSurface0 property={ surfaceColorType,1}
	ModifyGizmo ModifyObject=isoSurface0 property={ lineColorType,0}
	ModifyGizmo ModifyObject=isoSurface0 property={ lineWidthType,0}
	ModifyGizmo ModifyObject=isoSurface0 property={ fillMode,2}
	ModifyGizmo ModifyObject=isoSurface0 property={ lineWidth,1}
	ModifyGizmo ModifyObject=isoSurface0 property={ isoValue,0.5}
	ModifyGizmo ModifyObject=isoSurface0 property={ frontColor,1,0,0,1}
	ModifyGizmo ModifyObject=isoSurface0 property={ backColor,0,0,1,1}
	ModifyGizmo modifyObject=isoSurface0 property={calcNormals,1}
	AppendToGizmo light=Directional,name=light0
	ModifyGizmo light=light0 property={ liveUpdate,1}
	ModifyGizmo light=light0 property={ position,-0.272735,-0.570263,-0.774865,0.000000}
	ModifyGizmo light=light0 property={ direction,-0.272735,-0.570263,-0.774865}
	ModifyGizmo light=light0 property={ specular,1.000000,1.000000,1.000000,1.000000}
	ModifyGizmo setDisplayList=0, object=light0
	ModifyGizmo setDisplayList=1, opName=ortho0, operation=ortho, data={-1.61826,1.61826,-1.61826,1.61826,-2,2}
	ModifyGizmo setDisplayList=2, object=isoSurface0
	ModifyGizmo SETQUATERNION={0.696323,0.580486,0.415407,-0.074868}
	ModifyGizmo currentGroupObject=""
	ModifyGizmo compile

	ModifyGizmo userString={wmgizmo_df,"Gizmo0"}
	ModifyGizmo endRecMacro
End
window FIbrals3Ddisp(size, densitywaveloc) : GizmoPlot
	variable size
	string densitywaveloc
	PauseUpdate; Silent 1	// Building Gizmo 6 window...
	dowindow /k Fibrals3D
	// Do nothing if the Gizmo XOP is not available.
	if(exists("NewGizmo")!=4)
		DoAlert 0, "Gizmo XOP must be installed"
		return 0
	endif

	NewGizmo/N=Fibrals3D/T="3D Fibrals" /W=(30,30,1310,1054) /k=1
	ModifyGizmo startRecMacro
	ModifyGizmo scalingMode=8
	ModifyGizmo setOuterBox={0,size,0,size,0,size}
	ModifyGizmo scalingOption=60
	AppendToGizmo isoSurface=$densitywaveloc,name=isoSurface0
	ModifyGizmo ModifyObject=isoSurface0 property={ surfaceColorType,1}
	ModifyGizmo ModifyObject=isoSurface0 property={ lineColorType,0}
	ModifyGizmo ModifyObject=isoSurface0 property={ lineWidthType,0}
	ModifyGizmo ModifyObject=isoSurface0 property={ fillMode,2}
	ModifyGizmo ModifyObject=isoSurface0 property={ lineWidth,1}
	ModifyGizmo ModifyObject=isoSurface0 property={ isoValue,0.5}
	ModifyGizmo ModifyObject=isoSurface0 property={ frontColor,1,0,0,1}
	ModifyGizmo ModifyObject=isoSurface0 property={ backColor,0,0,1,1}
	ModifyGizmo modifyObject=isoSurface0 property={calcNormals,1}
	AppendToGizmo light=Directional,name=light0
	ModifyGizmo light=light0 property={ liveUpdate,1}
	ModifyGizmo light=light0 property={ position,-0.272735,-0.570263,-0.774865,0.000000}
	ModifyGizmo light=light0 property={ direction,-0.272735,-0.570263,-0.774865}
	ModifyGizmo light=light0 property={ specular,1.000000,1.000000,1.000000,1.000000}
	ModifyGizmo setDisplayList=0, object=light0
	ModifyGizmo setDisplayList=1, opName=ortho0, operation=ortho, data={-1.61826,1.61826,-1.61826,1.61826,-2,2}
	ModifyGizmo setDisplayList=2, object=isoSurface0
	ModifyGizmo SETQUATERNION={0.696323,0.580486,0.415407,-0.074868}
	ModifyGizmo currentGroupObject=""
	ModifyGizmo compile

	ModifyGizmo userString={wmgizmo_df,"Gizmo0"}
	ModifyGizmo endRecMacro
End

Window Spinoidal3Dview() : GizmoPlot
	PauseUpdate; Silent 1	// Building Gizmo 6 window...

	// Do nothing if the Gizmo XOP is not available.
	if(exists("NewGizmo")!=4)
		DoAlert 0, "Gizmo XOP must be installed"
		return
	endif
	
	NewGizmo/k=1/N=Spinoidal3Dview/T="Spinoidal View 3D"
	ModifyGizmo startRecMacro
	MoveWindow 44.25,126.5,462,444.5 
	ModifyGizmo scalingMode=8
	ModifyGizmo setOuterBox={0,127,0,127,0,100}
	ModifyGizmo scalingOption=15
	AppendToGizmo isoSurface=root:realu,name=isoSurface0
	ModifyGizmo ModifyObject=isoSurface0 property={ surfaceColorType,1}
	ModifyGizmo ModifyObject=isoSurface0 property={ lineColorType,0}
	ModifyGizmo ModifyObject=isoSurface0 property={ lineWidthType,0}
	ModifyGizmo ModifyObject=isoSurface0 property={ fillMode,2}
	ModifyGizmo ModifyObject=isoSurface0 property={ lineWidth,1}
	ModifyGizmo ModifyObject=isoSurface0 property={ isoValue,0.5}
	ModifyGizmo ModifyObject=isoSurface0 property={ frontColor,0,0,0.996109,1}
	ModifyGizmo ModifyObject=isoSurface0 property={ backColor,0.996109,0,0,1}
	ModifyGizmo modifyObject=isoSurface0 property={calcNormals,1}
	AppendToGizmo Axes=boxAxes,name=axes0
	ModifyGizmo ModifyObject=axes0,property={0,axisRange,-1,-1,-1,1,-1,-1}
	ModifyGizmo ModifyObject=axes0,property={1,axisRange,-1,-1,-1,-1,1,-1}
	ModifyGizmo ModifyObject=axes0,property={2,axisRange,-1,-1,-1,-1,-1,1}
	ModifyGizmo ModifyObject=axes0,property={3,axisRange,-1,1,-1,-1,1,1}
	ModifyGizmo ModifyObject=axes0,property={4,axisRange,1,1,-1,1,1,1}
	ModifyGizmo ModifyObject=axes0,property={5,axisRange,1,-1,-1,1,-1,1}
	ModifyGizmo ModifyObject=axes0,property={6,axisRange,-1,-1,1,-1,1,1}
	ModifyGizmo ModifyObject=axes0,property={7,axisRange,1,-1,1,1,1,1}
	ModifyGizmo ModifyObject=axes0,property={8,axisRange,1,-1,-1,1,1,-1}
	ModifyGizmo ModifyObject=axes0,property={9,axisRange,-1,1,-1,1,1,-1}
	ModifyGizmo ModifyObject=axes0,property={10,axisRange,-1,1,1,1,1,1}
	ModifyGizmo ModifyObject=axes0,property={11,axisRange,-1,-1,1,1,-1,1}
	ModifyGizmo ModifyObject=axes0,property={-1,axisScalingMode,1}
	ModifyGizmo ModifyObject=axes0,property={-1,axisColor,0,0,0,1}
	ModifyGizmo ModifyObject=axes0,property={0,ticks,3}
	ModifyGizmo ModifyObject=axes0,property={1,ticks,3}
	ModifyGizmo ModifyObject=axes0,property={2,ticks,3}
	ModifyGizmo modifyObject=axes0 property={Clipped,0}
	AppendToGizmo light=Directional,name=light0
	ModifyGizmo light=light0 property={ liveUpdate,1}
	ModifyGizmo light=light0 property={ position,-0.467891,-0.748782,-0.469471,0.000000}
	ModifyGizmo light=light0 property={ direction,-0.467891,-0.748782,-0.469471}
	ModifyGizmo light=light0 property={ specular,1.000000,1.000000,1.000000,1.000000}
	AppendToGizmo isoSurface=root:realu,name=isoSurface1
	ModifyGizmo ModifyObject=isoSurface1 property={ surfaceColorType,1}
	ModifyGizmo ModifyObject=isoSurface1 property={ lineColorType,0}
	ModifyGizmo ModifyObject=isoSurface1 property={ lineWidthType,0}
	ModifyGizmo ModifyObject=isoSurface1 property={ fillMode,2}
	ModifyGizmo ModifyObject=isoSurface1 property={ lineWidth,1}
	ModifyGizmo ModifyObject=isoSurface1 property={ isoValue,-0.5}
	ModifyGizmo ModifyObject=isoSurface1 property={ frontColor,0.996109,0,0,1}
	ModifyGizmo ModifyObject=isoSurface1 property={ backColor,0,0.996109,0,1}
	ModifyGizmo modifyObject=isoSurface1 property={calcNormals,1}
	ModifyGizmo setDisplayList=0, object=light0
	ModifyGizmo setDisplayList=1, object=isoSurface0
	ModifyGizmo setDisplayList=2, object=isoSurface1
	ModifyGizmo setDisplayList=3, object=axes0
	ModifyGizmo SETQUATERNION={0.251634,0.518248,0.722122,0.382934}
	ModifyGizmo currentGroupObject=""
	ModifyGizmo compile

	ModifyGizmo showInfo
	ModifyGizmo infoWindow={1481,523,1898,835}
	ModifyGizmo bringToFront
	ModifyGizmo userString={wmgizmo_df,"Gizmo0"}
	ModifyGizmo endRecMacro
End

window Spin3Ddisp(size, densitywaveloc)  : GizmoPlot
	variable size
	string densitywaveloc
	PauseUpdate; Silent 1	// Building Gizmo 6 window...
	dowindow /k Spheres3D
	// Do nothing if the Gizmo XOP is not available.
	if(exists("NewGizmo")!=4)
		DoAlert 0, "Gizmo XOP must be installed"
		return 0
	endif

	NewGizmo/K=1/N=Spheres3D/T="3D Spheres" /W=(29,54,1366,1024)
	ModifyGizmo startRecMacro
	ModifyGizmo scalingMode=8
	ModifyGizmo setOuterBox={0,size-1,0,size-1,0,size-1}
	ModifyGizmo scalingOption=15
	AppendToGizmo isoSurface=$densitywaveloc,name=isoSurface0
	ModifyGizmo ModifyObject=isoSurface0 property={ surfaceColorType,1}
	ModifyGizmo ModifyObject=isoSurface0 property={ lineColorType,0}
	ModifyGizmo ModifyObject=isoSurface0 property={ lineWidthType,0}
	ModifyGizmo ModifyObject=isoSurface0 property={ fillMode,2}
	ModifyGizmo ModifyObject=isoSurface0 property={ lineWidth,1}
	ModifyGizmo ModifyObject=isoSurface0 property={ isoValue,0.5}
	ModifyGizmo ModifyObject=isoSurface0 property={ frontColor,0,0,0.996109,1}
	ModifyGizmo ModifyObject=isoSurface0 property={ backColor,0.996109,0,0,1}
	ModifyGizmo modifyObject=isoSurface0 property={calcNormals,1}
	AppendToGizmo isoSurface=$densitywaveloc,name=isoSurface1
	ModifyGizmo ModifyObject=isoSurface1 property={ surfaceColorType,1}
	ModifyGizmo ModifyObject=isoSurface1 property={ lineColorType,0}
	ModifyGizmo ModifyObject=isoSurface1 property={ lineWidthType,0}
	ModifyGizmo ModifyObject=isoSurface1 property={ fillMode,2}
	ModifyGizmo ModifyObject=isoSurface1 property={ lineWidth,1}
	ModifyGizmo ModifyObject=isoSurface1 property={ isoValue,-0.5}
	ModifyGizmo ModifyObject=isoSurface1 property={ frontColor,0.996109,0,0,1}
	ModifyGizmo ModifyObject=isoSurface1 property={ backColor,0,0.996109,0,1}
	ModifyGizmo modifyObject=isoSurface1 property={calcNormals,1}
	AppendToGizmo light=Directional,name=light0
	ModifyGizmo light=light0 property={ liveUpdate,1}
	ModifyGizmo light=light0 property={ position,-0.467891,-0.748782,-0.469471,0.000000}
	ModifyGizmo light=light0 property={ direction,-0.467891,-0.748782,-0.469471}
	ModifyGizmo light=light0 property={ specular,1.000000,1.000000,1.000000,1.000000}
	ModifyGizmo setDisplayList=0, opName=ortho0, operation=ortho, data={-1.48235,1.48235,-2.0278,0.936895,-2,2}
	ModifyGizmo setDisplayList=1, object=light0
	ModifyGizmo setDisplayList=2, object=isoSurface0
	ModifyGizmo setDisplayList=3, object=isoSurface1
	ModifyGizmo SETQUATERNION={0.139448,0.513616,0.801278,0.273315}
	ModifyGizmo currentGroupObject=""
	ModifyGizmo compile

	ModifyGizmo userString={wmgizmo_df,"Gizmo0"}
	ModifyGizmo endRecMacro
End
function makespinoidalslices(wavein)
	wave wavein
	dowindow /k scaledSlice
	Display/k=1/n=ScaledSlice /W=(395,44,551,200)
	AppendImage /w=ScaledSlice/T wavein
	ModifyImage /w=ScaledSlice $nameofwave(wavein) ctab= {-1,1,BlueRedGreen,0}
	ModifyGraph /w=ScaledSlice margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=14
	ModifyGraph /w=ScaledSlice mirror=2
	ModifyGraph /w=ScaledSlice nticks=4
	ModifyGraph /w=ScaledSlice minor=1
	ModifyGraph /w=ScaledSlice fSize=9
	ModifyGraph /w=ScaledSlice standoff=0
	ModifyGraph /w=ScaledSlice tkLblRot(left)=90
	ModifyGraph /w=ScaledSlice btLen=3
	ModifyGraph /w=ScaledSlice tlOffset=-2
	SetAxis/w=ScaledSlice/A/R left
	SetAxis/w=ScaledSlice left 0,100
	SetAxis/w=ScaledSlice top 0,100
	
	dowindow /k UnScaledSlice
	Display/k=1/n=UnScaledSlice /W=(40,44,196,200)
	AppendImage /w=UnScaledSlice/T wavein
	ModifyImage /w=UnScaledSlice $nameofwave(wavein) ctab= {*,*,Grays,0}
	ModifyGraph /w=UnScaledSlice margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=14
	ModifyGraph /w=UnScaledSlice mirror=2
	ModifyGraph /w=UnScaledSlice nticks=4
	ModifyGraph /w=UnScaledSlice minor=1
	ModifyGraph /w=UnScaledSlice fSize=9
	ModifyGraph /w=UnScaledSlice standoff=0
	ModifyGraph /w=UnScaledSlice tkLblRot(left)=90
	ModifyGraph /w=UnScaledSlice btLen=3
	ModifyGraph /w=UnScaledSlice tlOffset=-2
	SetAxis/w=UnScaledSlice/A/R left
	SetAxis/w=UnScaledSlice left 0,100
	SetAxis/w=UnScaledSlice top 0,100
end

Window Spinoidal3DLayout() : Layout
	dowindow /k Spinoidal3DLayout
	PauseUpdate; Silent 1		// building window...
	Layout/k=1/C=1/W=(100,100,782,600) Spinoidal3DImage(14,23,599,410)/O=1/F=0 as "Spinoidal3DLayout"
	ModifyLayout mag=1
EndMacro
Window Spinoidal3DLayout_2() : Layout
	dowindow /k Spinoidal3DLayout_2
	PauseUpdate; Silent 1		// building window...
	Layout/C=1/K=1/W=(111.75,116,762,692) Spinoidal3DImage(4.5,159,578.25,554.25)/O=1/F=0 as "Spinoidal3DLayout"
	printsettings /i /w=Spinoidal3DLayout_2 margins={.25 , .25 , .25 , .25 }
	Append UnscaledSlice(44.25,18.75,296.25,270.75)/O=1/F=0,ScaledSlice(324.75,18.75,576.75,270.75)/O=1/F=0
	TextBox/C/N=text0/D={1,8,-1}/F=0/B=1/A=LB/X=1.31/Y=64.71 "\\Z24\\F'Arial Bold'Relative Concentration"
	TextBox/C/N=text1/D={1,8,-1}/F=0/B=1/A=LB/X=63.84/Y=64.71 "\\Z24\\F'Arial Bold'Absolute"
	ModifyLayout mag=1
EndMacro
Window Spinoidal3DImage(imagename) : Graph
	string imagename
	dowindow /k Spinoidal3DImage
	PauseUpdate; Silent 1		// building window...
	Display /k=1/n=Spinoidal3DImage/W=(70,70,70+1281,70+1024) as "Spinoidal3DImage"
	AppendImage/T $imagename
	ModifyImage testimage ctab= {*,*,Grays,0}
	ModifyGraph margin(left)=1,margin(bottom)=1,margin(top)=1,margin(right)=1
	ModifyGraph mirror=2
	ModifyGraph nticks=0
	ModifyGraph minor=1
	ModifyGraph noLabel=2
	ModifyGraph fSize=8
	ModifyGraph standoff=0
	ModifyGraph axThick=0
	ModifyGraph tkLblRot(left)=90
	ModifyGraph btLen=3
	ModifyGraph tlOffset=-2
	SetAxis/A/R left
EndMacro





function /s addcontrols(controllist)
	string controllist
	variable i, height = 110, width = 30
	string controltype, titlestr, controlstring, val, varname, ctrlname, returnlist=""
	for(i=0;i<itemsinlist(controllist);i+=1)
		controlstring = stringfromlist(i,controllist)
		titlestr = stringfromlist(0,controlstring,",")
		controltype = stringfromlist(1,controlstring,",")
		val = stringfromlist(2,controlstring,",")
		varname=""
		ctrlname=""
		if(stringmatch("checkbox",controltype))
			varname = uniquename(cleanupname(titlestr,0),3,1,"ScatteringSimPanel")
			variable /g $varname = (str2num(val))
			ctrlname = uniquename(cleanupname(titlestr,0),15,1,"ScatteringSimPanel")
			checkbox $ctrlname win=ScatteringSimPanel, title=titlestr, pos={width,height}, size={200,25}, variable=$varname
		elseif(stringmatch("SetVariable",controltype))
			varname = uniquename(cleanupname(titlestr,0),3,1,"ScatteringSimPanel")
			variable /g $varname = (str2num(val))
			ctrlname = uniquename(cleanupname(titlestr,0),15,1,"ScatteringSimPanel")
			SetVariable $ctrlname win=ScatteringSimPanel, title=titlestr, pos={width,height}, size={200,25}, variable=$varname
		elseif(stringmatch("String",controltype))
			varname = uniquename(cleanupname(titlestr,0),3,1,"ScatteringSimPanel")
			string /g $varname = val
			ctrlname = uniquename(cleanupname(titlestr,0),15,1,"ScatteringSimPanel")
			SetVariable $ctrlname win=ScatteringSimPanel, title=titlestr, pos={width,height}, size={200,25}, variable=$varname
		endif
		returnlist = returnlist + ctrlname + "," + varname+"," + controltype + ";"
		height +=30
		if(height>380)
			width = 240
			height = 110
		endif
	endfor
	return returnlist
end
function removecontrols(controllist)
	string controllist
	string controltype, titlestr, controlstring, val, varname, ctrlname, returnlist=""
	variable i
	for(i=0;i<itemsinlist(controllist);i+=1)
		controlstring = stringfromlist(i,controllist)
		ctrlname = stringfromlist(0,controlstring,",")
		varname = stringfromlist(1,controlstring,",")
		controltype = stringfromlist(2,controlstring,",")
		if(stringmatch("checkbox",controltype))
			killcontrol /w=ScatteringSimPanel $ctrlname
			killvariables /z $varname
		elseif(stringmatch("SetVariable",controltype))
			killcontrol /w=ScatteringSimPanel $ctrlname
			killvariables /z $varname
		elseif(stringmatch("String",controltype))
			killcontrol /w=ScatteringSimPanel $ctrlname
			killstrings /z $varname
		endif
	endfor
end



Function Model3DProc(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// move material to list of materials
			string foldersave = getdatafolder(1)
			setdatafolder root:Packages:ScatterSim3D
			wave/t listofmaterials
			wave selectedmaterials
			svar material,alignment
			nvar alignmentstrength, alignmentsize
//			insertpoints dimsize(listofmaterials,0)+1 ,1, listofmaterials
			variable lastindex = dimsize(listofmaterials,0)
			redimension /n=(lastindex+1,2) listofmaterials, selectedmaterials
			listofmaterials[lastindex][0] = material
			if(stringmatch("none",alignment))
				listofmaterials[lastindex][1] = alignment
			else
				listofmaterials[lastindex][1] = alignment +","+ num2str(alignmentsize)+" pixels," +num2str(alignmentstrength) 
			endif
			selectedmaterials = p==lastindex ? 1 : 0
			setdatafolder foldersave
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function RemoveMaterialButton(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// Remove Selected Material from the list
			string foldersave = getdatafolder(1)
			setdatafolder root:Packages:ScatterSim3D
			variable i
			wave selectedmaterials
			wave /t listofmaterials
			for(i=dimsize(selectedmaterials,0);i>=0;i-=1)
				if(selectedmaterials[i][0])
					deletepoints/M=0 i,1, listofmaterials, selectedmaterials
				endif
			endfor
			if(numpnts(listofmaterials)<1)
				make/n=(0,2)/o/t listofmaterials
				make/n=(0,2)/o selectedmaterials
			endif
			setdatafolder foldersave
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function AlignmentPop(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			string foldersave = getdatafolder(1)
			setdatafolder root:Packages:ScatterSim3D
			svar alignment
			alignment = popstr
			setdatafolder foldersave
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function MaterialPop(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			string foldersave = getdatafolder(1)
			setdatafolder root:Packages:ScatterSim3D
			svar material, alignment
			material = popstr
			nvar usemorphalignment
			if(OCsexist(material)!=2  || usemorphalignment)
				PopupMenu AlignmentPop disable=2, mode=1
				alignment="None"
				SetVariable AlignmentSize disable=2
				SetVariable AlignmentStrength disable=2
			else
				PopupMenu AlignmentPop disable=0
				SetVariable AlignmentSize disable=0
				SetVariable AlignmentStrength disable=0
			endif
			
			
			
			setdatafolder foldersave
			break
			
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function SimulateSystem(ba) : ButtonControl
	STRUCT WMButtonAction &ba

	switch( ba.eventCode )
		case 2: // mouse up
			// SimulateSystem
			string foldersave = getdatafolder(1)
			setdatafolder root:Packages:ScatterSim3D
			svar controllist
			wave /t listofmaterials
			string paramstring = generateParamstring(controllist)
			string materialstring = generatematerialstring(listofmaterials)
			nvar usealignment = root:Packages:ScatterSim3D:useprecalcalignment
			string alignmentstring 
			if(usealignment)
				alignmentstring="None"
			else
				alignmentstring = generatealignmentstring(listofmaterials)
			endif
			svar SimName
			svar morphology
//			nvar alignment
			nvar sizescale
			nvar resolution
			nvar numensteps
			nvar movie
			nvar startxrayenergy
			nvar endxrayenergy
			nvar savescattering2D
			nvar voxelsize
			nvar thickness
			nvar rotatesys
//			nvar useprecalcalignment
			svar efield
			string ycomp, zcomp
			splitstring /e="^[( ]*([1234567890.]*)[^,]{0,5},[^,0.1]{0,5}([1234567890.]*)[) ]*" efield, ycomp, zcomp
			newdatafolder /o/s $cleanupname(simname,0)
			model3D(morphology,voxelsize,sizescale,resolution,thickness,paramstring,materialstring,alignmentstring,ycomp+","+zcomp,startxrayenergy,endxrayenergy,numensteps,movie = movie, save3d = savescattering2D, rotatesys = Rotatesys)
			setdatafolder foldersave
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End
function /s generatematerialstring(materialswave)
	wave /t materialswave
	string returnlist=""
	variable i
	for(i=0;i<dimsize(materialswave,0);i+=1)
		returnlist =returnlist + materialswave[i][0]+","
	endfor
	return returnlist
end
function /s generatealignmentstring(materialswave)
	wave /t materialswave
	string alignments, strengths, widths, returnlist=""
	variable i
	for(i=0;i<dimsize(materialswave,0);i+=1)
		if(stringmatch(materialswave[i][1] ,"*none*"))
			alignments = "-1"
			strengths = "0"
			widths = "0"
		elseif(stringmatch(materialswave[i][1] ,"*Face*"))
			alignments = "1"
			print materialswave[i][1]
			splitstring /e="^[^,]*,([^,]*) pixels,([^,]*)" materialswave[i][1], widths, strengths
		elseif(stringmatch(materialswave[i][1] ,"*Edge*"))
			alignments = "0"
			splitstring /e="^[^,]*,([^,]*) pixels,([^,]*)" materialswave[i][1], widths, strengths
		endif
		returnlist =returnlist + alignments+"," +widths+","+ strengths+";"
	endfor
	return returnlist
end

function /s generateparamstring(controllist)
	string controllist
	string controltype, val, varname, returnlist=""
	variable i
	for(i=0;i<itemsinlist(controllist);i+=1)
		varname = stringfromlist(1,stringfromlist(i,controllist),",")
		controltype = stringfromlist(2,stringfromlist(i,controllist),",")
		if(stringmatch("checkbox",controltype))
			nvar nvariable = $varname
			returnlist = returnlist + num2str(nvariable) + ","
		elseif(stringmatch("SetVariable",controltype))
			nvar nvariable = $varname
			returnlist = returnlist + num2str(nvariable) + ","
		elseif(stringmatch("String",controltype))
			svar svariable = $varname
			returnlist = returnlist + svariable + "," 
		endif
	endfor
	return returnlist
end

Function PopMorph(pa) : PopupMenuControl
	STRUCT WMPopupAction &pa

	switch( pa.eventCode )
		case 2: // mouse up
			Variable popNum = pa.popNum
			String popStr = pa.popStr
			string foldersave = getdatafolder(1)
			setdatafolder root:packages:ScatterSim3D
			svar morphology
			morphology = popstr
			svar controllist, extralist
			extralist = ""
			nvar useprecalcalignment, usemorphalignment
			removecontrols(controllist)
			execute/Q/Z ("controllist = variables_"+popstr+"()")
			execute/Q/Z ("extralist = special_"+popstr+"()")
			if(Stringmatch(extralist,"*IncludesAlignment*"))
				usemorphalignment=1
			else
				usemorphalignment=0
			endif
			useprecalcalignment=0
			controllist = addcontrols(controllist)
			setdatafolder foldersave
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End

Function SetEField(sva) : SetVariableControl
	STRUCT WMSetVariableAction &sva

	switch( sva.eventCode )
		case 1: // mouse up
		case 2: // Enter key
		case 3: // Live update
			Variable dval = sva.dval
			String sval = sva.sval
			string foldersave = getdatafolder(1)
			setdatafolder root:packages:ScatterSim3D
			svar efield
			string ycomps, zcomps
			splitstring /e="^[( ]*([1234567890.]*)[^,]{0,5},[^,0.1]{0,5}([1234567890.]*)[) ]*" efield, ycomps, zcomps
			variable ycomp = str2num(ycomps)
			variable zcomp = str2num(zcomps)
			variable mag
			if(ycomp*0==0) // y is ok
				if(zcomp*0==0) // y and z is ok
					mag = sqrt(ycomp^2 + zcomp^2)
					ycomp/=mag
					zcomp/=mag
				else // y is ok by z isn't
					ycomp=1
					zcomp=0
				endif
			else //y is bad
				if(zcomp*0==0) // y is bad and z is ok
					ycomp=0
					zcomp=1
				else // y and z aren't ok
					ycomp=1
					zcomp=0
				endif
			endif
			//for now we can only do z fields
			zcomp=1
			ycomp=0
			
			sprintf efield " ( %1.4f , %1.4f )" ycomp, zcomp 
			
			setdatafolder foldersave
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End


Function UserPreCalcChk(cba) : CheckBoxControl
	STRUCT WMCheckboxAction &cba

	switch( cba.eventCode )
		case 2: // mouse up
			Variable checked = cba.checked
			popupmenu MaterialPop,disable= (checked)
			popupmenu AlignmentPop,disable= (checked)
			Setvariable AlignmentSize, disable= (checked)
			Setvariable AlignmentStrength, disable= (checked)
			button AddMaterial,disable= (checked)
			button CalcAlignBut,disable= (checked)
			button removematerial,disable=(checked)
			listbox materials ,disable=(checked)
			groupbox MaterialsListText,disable=(checked)
			break
		case -1: // control being killed
			break
	endswitch

	return 0
End
function /s variables_fibrals()
	string variables
	variables = "Interpenetration [voxels],SetVariable,2;"
	variables = variables + "Number of Fibrals (Max),SetVariable,500;"
	variables = variables + "PolyDispursity (sigma of radiuses),setvariable,1;"
	variables = variables + "Minimum Radius,SetVariable,1;"
	variables = variables + "Maximum Radius,SetVariable,5;"
	variables = variables + "Sigma (around x-plane in degrees),SetVariable,5;"
	variables = variables + "Volume Fraction (<1),SetVariable,.5;"
	variables = variables + "Noise,SetVariable,0;"
	variables = variables + "Minimum Seperation,SetVariable,0.3;"
	variables = variables + "Minimum Fibral Length,SetVariable,5;"
	variables = variables + "Maximum Fibral Length,SetVariable,40;"
	variables = variables + "ShellAlignment,SetVariable,1;"
	return variables
end
function /s special_fibrals()
	return "IncludesAlignment"
end
function model3D_fibrals(s3d)

	struct ThreeDSystem &s3d
	if(itemsinlist(s3d.paramstring,",")<11)
		return -1
	endif
	newdatafolder /o/s CreatingFibrals
	variable interpenetration = 	str2num(stringfromlist( 0 ,s3d.paramstring,","))
	variable Fibralnum = 		str2num(stringfromlist( 1 ,s3d.paramstring,","))
	variable pd = 				str2num(stringfromlist( 2 ,s3d.paramstring,","))
	variable thickness = 		s3d.thickness
	variable minsize = 		str2num(stringfromlist( 3 ,s3d.paramstring,","))
	variable maxsize = 		str2num(stringfromlist( 4 ,s3d.paramstring,","))
	variable angsigma = 		str2num(stringfromlist( 5 ,s3d.paramstring,","))
	variable volfrac = 			str2num(stringfromlist( 6 ,s3d.paramstring,","))
	variable noise = 			str2num(stringfromlist( 7 ,s3d.paramstring,","))
	variable minsep = 		str2num(stringfromlist( 8 ,s3d.paramstring,","))
	variable minlength = 		str2num(stringfromlist( 9 ,s3d.paramstring,","))
	variable maxlength = 		str2num(stringfromlist( 10 ,s3d.paramstring,","))
	variable shellalignment = str2num(stringfromlist( 11 ,s3d.paramstring,","))
	
	make /n=(thickness*s3d.num^2,3)/o rmat
	make /o /n=(thickness,s3d.num,s3d.num) mat=0,xwave, ywave, zwave, ammat=0
	make /n=( thickness,s3d.num,s3d.num,3)/o vecmat=0
	if(minsep<1)
		make/B/U /o /n=(thickness,s3d.num,s3d.num,maxsize) exmat =  (p <= t/(1+minsep)) || (q <= t/(1+minsep)) || (r <= t/(1+minsep) ) || (p >= thickness-t/(1+minsep)) || (q >= s3d.num-t/(1+minsep)) || (r >= s3d.num-t/(1+minsep)) ? 0 : 1
	else
		make/B/U /o /n=(thickness,s3d.num,s3d.num,maxsize) exmat =  (p <= t-minsep) || (q <= t-minsep) || (r <= t-minsep) || (p >= thickness-t+minsep) || (q >= s3d.num-t+minsep) || (r >= s3d.num-t+minsep) ? 0 : 1
	endif
	make/B/U /o /n=(thickness,s3d.num,s3d.num) tempwave, tempx
	if(s3d.movie)
		Execute("Fibrals3Ddisp(" +num2str(s3d.num)+", \""+getwavesdatafolder(mat,2)+"\")")
		Execute("exportgizmo wave=\"testimage\"   ;Spinoidal3DLayout();Spinoidal3DImage(\""+getdatafolder(1)+"testimage\")")
	endif

	xwave = x
	ywave = y
	zwave = z
	redimension /n=(thickness*s3d.num*s3d.num) xwave, ywave, zwave
	rmat[][0] = xwave[p]
	rmat[][1] = ywave[p]
	rmat[][2] = zwave[p]
	variable testcx,testcy,testcz,i,radius, orad, cx, cy, cz, failed, fnum =0, f2num=0, xmn,xmx,ymn,ymx,zmn,zmx, loc, qfnum=0, theta, phi
	variable tx,ty,tz,hit, mx=thickness-1, my=s3d.num-1, mz=s3d.num-1, length
	fnum=0
	make /o /n=3 vec
	for(i=0;i<Fibralnum;i+=1)
		radius=max(1,gnoise(pd)+s3d.size)
		tempx = exmat[p][q][r][ceil(radius)]
		duplicate /o tempx, tempwave
		redimension /n=(numpnts(tempwave)) tempwave
		integrate tempwave /D=intwave
		if(wavemax(intwave) < 5)
			fnum+=1
			if(fnum<5)
				continue
			else
				print "warning:  can't fit in anymore fibrals, only " + num2str(i-fnum) + " fibrals were created"
				break
			endif
		endif
		loc = binarysearch(intwave, enoise(wavemax(intwave)/2)+wavemax(intwave)/2)
		
		cx = rmat[loc][0]
		cy = rmat[loc][1]
		cz = rmat[loc][2]
		
		theta = (90 + gnoise(angsigma)) * pi / 180 // polar angle of the fibral axis relative to the film normal (right now it is mostly in plane)
		phi = enoise(pi)  // (azimuthal angle of the axis of the fibral)
		vec[2] = sin(theta)*cos(phi)
		vec[1] = sin(theta)*sin(phi)
		vec[0] = cos(theta)
		
		// for Micelles, change the above angles, so that they are oriented how you want in the plane (right now the 
		
		//forward direction
		hit=0
		length=0
		tx=cx // center location
		ty=cy
		tz=cz
		do
			
			if(shellalignment>0)
				vecmat[max(tx-radius,0),min(tx+radius,mx)][max(ty-radius,0),min(ty+radius,my)][max(tz-radius,0),min(tz+radius,mz)][]= (p-tx)^2 + (q-ty)^2 + (r-tz)^2 < radius^2  && ammat[p][q][r]==0? vec[t] : vecmat
				vecmat[max(tx-radius,0),min(tx+radius,mx)][max(ty-radius,0),min(ty+radius,my)][max(tz-radius,0),min(tz+radius,mz)][]= (p-tx)^2 + (q-ty)^2 + (r-tz)^2 < (radius-1.5)^2 ? 0 : vecmat
				ammat[max(tx-radius,0),min(tx+radius,mx)][max(ty-radius,0),min(ty+radius,my)][max(tz-radius,0),min(tz+radius,mz)]= (p-tx)^2 + (q-ty)^2 + (r-tz)^2 < (radius-1.5)^2 ? 1 : ammat
			elseif(shellalignment < 0)
				vecmat[max(tx-radius,0),min(tx+radius,mx)][max(ty-radius,0),min(ty+radius,my)][max(tz-radius,0),min(tz+radius,mz)][]= (p-tx)^2 + (q-ty)^2 + (r-tz)^2 < (radius-1.5)^2 ? vec[t] : vecmat
				ammat[max(tx-radius,0),min(tx+radius,mx)][max(ty-radius,0),min(ty+radius,my)][max(tz-radius,0),min(tz+radius,mz)]   = (p-tx)^2 + (q-ty)^2 + (r-tz)^2 > (radius-1.5)^2 &&  (p-tx)^2 + (q-ty)^2 + (r-tz)^2 < (radius)^2 && abs(vecmat[p][q][r][0])+abs(vecmat[p][q][r][1])+abs(vecmat[p][q][r][2])==0? 1 : ammat
			endif
			// for micelles change the above line, (right now it's just making the dipole alignment in line with the fibral main axis)
			// radial, and two materials
			// need to add another material outer core, and add to exmat, probably with a core shell diameter parameters needed
			
			
			exmat[max(tx-radius-maxsize,0),min(tx+radius+maxsize,mx)][max(ty-radius-maxsize,0),min(ty+radius+maxsize,my)][max(tz-radius-maxsize,0),min(tz+radius+maxsize,mz)][]= (p-tx)^2 + (q-ty)^2 + (r-tz)^2 < (radius+t)^2 ? 0 : exmat
			tx += vec[0]
			ty += vec[1]
			tz += vec[2]
			hit = 1-tempx[Min(max(tx,0),mx)][Min(max(ty,0),my)][Min(max(tz,0),mz)]
			length +=1
		while( length<maxlength &&  hit==0 && tx <= mx && ty <= my && tz <= mz && tx>=0 && ty>=0 && tz>=0 )
		//other direction
		hit=0
		tx=cx
		ty=cy
		tz=cz
		length-=1
		do
			if(shellalignment>0)
				vecmat[max(tx-radius,0),min(tx+radius,mx)][max(ty-radius,0),min(ty+radius,my)][max(tz-radius,0),min(tz+radius,mz)][]= (p-tx)^2 + (q-ty)^2 + (r-tz)^2 < radius^2  && ammat[p][q][r]==0? vec[t] : vecmat
				vecmat[max(tx-radius,0),min(tx+radius,mx)][max(ty-radius,0),min(ty+radius,my)][max(tz-radius,0),min(tz+radius,mz)][]= (p-tx)^2 + (q-ty)^2 + (r-tz)^2 < (radius-1.5)^2 ? 0 : vecmat
				ammat[max(tx-radius,0),min(tx+radius,mx)][max(ty-radius,0),min(ty+radius,my)][max(tz-radius,0),min(tz+radius,mz)]= (p-tx)^2 + (q-ty)^2 + (r-tz)^2 < (radius-1.5)^2 ? 1 : ammat
			elseif(shellalignment < 0)
				vecmat[max(tx-radius,0),min(tx+radius,mx)][max(ty-radius,0),min(ty+radius,my)][max(tz-radius,0),min(tz+radius,mz)][]= (p-tx)^2 + (q-ty)^2 + (r-tz)^2 < (radius-1.5)^2 ? vec[t] : vecmat
				ammat[max(tx-radius,0),min(tx+radius,mx)][max(ty-radius,0),min(ty+radius,my)][max(tz-radius,0),min(tz+radius,mz)]   = (p-tx)^2 + (q-ty)^2 + (r-tz)^2 > (radius-1.5)^2 &&  (p-tx)^2 + (q-ty)^2 + (r-tz)^2 < (radius)^2 && abs(vecmat[p][q][r][0])+abs(vecmat[p][q][r][1])+abs(vecmat[p][q][r][2])==0? 1 : ammat
			endif
			
			exmat[max(tx-radius-maxsize,0),min(tx+radius+maxsize,mx)][max(ty-radius-maxsize,0),min(ty+radius+maxsize,my)][max(tz-radius-maxsize,0),min(tz+radius+maxsize,mz)][]= (p-tx)^2 + (q-ty)^2 + (r-tz)^2 < (radius+t)^2 ? 0 : exmat
			tx -= vec[0]
			ty -= vec[1]
			tz -= vec[2]
			hit = 1-tempx[Min(max(tx,0),mx)][Min(max(ty,0),my)][Min(max(tz,0),mz)] // if we hit another fibral or the edge, tempx will be 0 (hit will be 1)  -  until then, hit will equal 0
			length+=1
		while( length<maxlength && hit==0 && tx <= mx && ty <= my && tz <= mz && tx>=0 && ty>=0 && tz>=0 )

		mat = sqrt( vecmat[p][q][r][0]^2 + vecmat[p][q][r][1]^2 + vecmat[p][q][r][2]^2 + ammat[p][q][r])
		imagefilter /n=3 /o gauss3d, mat
		if(s3d.movie)
			execute("ModifyGizmo /n=Fibrals3D update=2")
			doupdate
			Execute "exportgizmo wave=\"testimage\"   "
			//TextBox/w=Spinoidal3DLayout/C/N=text0/A=LT/X=0.00/Y=0.00 "\Z32" + time2str2(ttot)
			doupdate
			savepict /p=_PictGallery_ /E=-5 /N=Spinoidal3DLayout /o as "Frame3D"
			addmovieframe /pict=Frame3D
		endif
		
		if(length < minlength)
			//print "Created a fibral too short"
			f2num+=1
		endif
		
		if(f2num>5 + i/10) // if ever more than 10% + 5 are short, then let's stop
			print "too many short fibrals have been created, ending fibral creation"
			break
		endif
		
	endfor
	if(interpenetration > 1)
		mat = vecmat[p][q][r][0]
		imagefilter /n=(interpenetration) /o gauss3d, mat
		vecmat[][][][0] = mat[p][q][r]
		mat = vecmat[p][q][r][1]
		imagefilter /n=(interpenetration) /o gauss3d, mat
		vecmat[][][][1] = mat[p][q][r]
		mat = vecmat[p][q][r][2]
		imagefilter /n=(interpenetration) /o gauss3d, mat
		vecmat[][][][2] = mat[p][q][r]
	endif
	mat = sqrt( vecmat[p][q][r][0]^2 + vecmat[p][q][r][1]^2 + vecmat[p][q][r][2]^2 + ammat[p][q][r])
	
	setdatafolder ::
	variable fibralvol = mean(mat)
	variable rhomatrix = (volfrac - fibralvol)/(1-fibralvol)
	if(rhomatrix < 0 )
		print "Volume fraction is too low, make less fibrals"
		return -1
	endif
	make /n=(thickness,s3d.num,s3d.num,4) /o m1=0, m2=0
	wave s3d.m1=m1, s3d.m2=m2
	s3d.m1[][][][0,2] = vecmat[p][q][r][t]
	s3d.m1[][][][3] = rhomatrix * (1-mat[p][q][r]- ammat[p][q][r]) + ammat[p][q][r]
	s3d.m2[][][][3] = (1-rhomatrix) * (1-mat[p][q][r])
	
	duplicate /o mat,s3d.density1 // this returns the density matrix of material 1 (the matrix) for alignment etc later on
end

function /s variables_lamella()
	string variables
	variables = "Nucleation Probability [%],SetVariable,0.002;" 							//		   0
	variables = variables + "Xtal thickness in lamella [%],SetVariable,67;" 			//		   1
	variables = variables + "Xtal variability sigma [px],SetVariable,1.5;" 			//		   2
	variables = variables + "Amorph Alignment [-1 1][L ||],SetVariable, -.5;" 		//		   3
	variables = variables + "Amorph Alignment width [px],SetVariable,3;" 				//		   4
	variables = variables + "Amorphous variability σ [px],SetVariable,2;" 				//		   5
 	variables = variables + "Inter-lamella angular σ [rad],SetVariable,0.4;" 			//		   6
	variables = variables + "Intra-lamella angular σ [rad],SetVariable,0.04;" 		//		   7
	variables = variables + "Relative crystalline density [delta%},SetVariable,10;" //		   8
	variables = variables + "Grain Boundary [px],SetVariable,1;" 							//		   9
	variables = variables + "Intra-to-Inter lamella growth factor,SetVariable,100;" //		  10
	variables = variables + "Lamella region volume [%},SetVariable,80;" 				//		  11
	variables = variables + "Amorphous region thickness [%},SetVariable,80;" 			//		  12
	variables = variables + "Material 2 in amorphous region [%},SetVariable,0;" 		//		  13
	variables = variables + "Test threshold for Edge finding [%],SetVariable,50;" 	//		  14
	variables = variables + "Lamella Edge Density [%],SetVariable,90;" 					//		  15
	variables = variables + "swap fringe/amorph alignment [0 1],SetVariable,0;" 	//		  16
	return variables
end
function /s special_lamella()
	return "IncludesAlignment"
end


function model3D_lamella(s3d)

	struct ThreeDSystem &s3d
	if(itemsinlist(s3d.paramstring,",")<14)
		return -1
	endif
	variable nucprob = 			str2num(stringfromlist( 0 ,s3d.paramstring,","))
	variable xtallampct = 			str2num(stringfromlist( 1 ,s3d.paramstring,","))
	variable lamthickness = s3d.size * xtallampct/100//str2num(stringfromlist( 1 ,s3d.paramstring,",")) 
	variable lamuncert = 			str2num(stringfromlist( 2 ,s3d.paramstring,","))
	variable nz = s3d.thickness
	variable amorphthickness = 	s3d.size * (1-xtallampct/100)//str2num(stringfromlist( 3 ,s3d.paramstring,","))
	variable amorphalignment =	str2num(stringfromlist( 3 ,s3d.paramstring,","))
	variable amorphAlignWid =	str2num(stringfromlist( 4 ,s3d.paramstring,","))
	variable amorphuncert = 		str2num(stringfromlist( 5 ,s3d.paramstring,","))
	variable angularuncert = 	str2num(stringfromlist( 6 ,s3d.paramstring,","))
	variable XtalAngUncert = 	str2num(stringfromlist( 7 ,s3d.paramstring,","))
	variable xtaldens = 			str2num(stringfromlist( 8 ,s3d.paramstring,","))
	variable GBound = 				str2num(stringfromlist( 9 ,s3d.paramstring,","))
	variable GLamSpd = 			str2num(stringfromlist( 10 ,s3d.paramstring,","))
	variable Xtalpct = 			str2num(stringfromlist( 11 ,s3d.paramstring,","))
	variable regiondens = 		str2num(stringfromlist( 12 ,s3d.paramstring,","))
	variable mat2inamorph = 		str2num(stringfromlist( 13 ,s3d.paramstring,","))/100
	variable edgethreshold = 	str2num(stringfromlist( 14 ,s3d.paramstring,","))/100
	variable Mat2Density = 		str2num(stringfromlist( 15 ,s3d.paramstring,","))/100
	variable swapfringe = 		str2num(stringfromlist( 16 ,s3d.paramstring,","))
	
	
	variable nxy = s3d.num
	
	// outputs are m1-4 and density1-3 
	//m1=lamella aligned, 
	//m2=lamella edge perp, 
	//m3=material 2 in amorphous phase, 
	//m4=vacuum (roughness)

	variable nnd = 3
	
	variable stopsum = 255-(2.55*Xtalpct) // average occupation level to stop at (out of 255) (~5 %)
	variable zratio = 100 // ratio of z probability to x and y (level of in plane alignment is here
	
	string foldersave = getdatafolder(1)
	killwindow /z progresspnl
	if(s3d.movie)
		killwindow /z alignmentmap
	endif
	make /d/o/n=(nz,nxy,nxy,4) dens1,mat,mat2 // the density and m waves required for calculations
	newdatafolder /o/s working
	// make system
	make /d/o/n=(nxy,nxy,nz) alignx=0, aligny=0, alignz=0, nucmap, alignmag=0, thickness=0, lcx=0, lcy=0, lcz=0, edge=0, gradx, grady, gradz
	make /o/U/B /n=(nxy,nxy,nz) empty=255
	
	make /n=((nnd*2+1)^3 -1) /o/d nndist, nnx,nny,nnz, nnoffx, nnoffy, nnoffz, annx, anny, annz, nnvalid, nnprop, nnlcx,nnlcy,nnlcz,nnmag, nnaligned
			// nndist,nnoffx, nnoffy, nnoffz  // fixed lookuptables for all nearest neighbors - do not change
			// nnx,nny,nnz //  calculated nearest neighbor locations
			// annx,anny,annz // alignment at the nearest neighbor location
			// nnaligned // is the nn aligned with the pixel? - if not we want to make sure the result is amorphous
			// nnvalid // nearest neighbor is valid
			// nnmag // nearest neighbor is crystalline
			// nnprop // the ideality of propagating from this pixel (has to be valid, have good alignment, be closest to direction of propagtion)
			// nnlcx,y,z // store the crystalline center from each nn pixel
			
	variable emptysum=255, edgesum=0, j, xloc, yloc, zloc, xoff, yoff, zoff, xeff, yeff, zeff, angle
	variable tempsum, temppnts, tempapnts, gx, gy, gz, xtalend, xtallen, xtallenpos
	variable n1x, n1y, n1z, n2x, n2y, n2z, n1off, n2off, poff, bx, by, bz, xtalavg, pindex, newthickness, temp, localalignment, amorphprop, noprogress=0,lastedgeloc=-1
	for(xoff=-nnd;xoff<=nnd;xoff+=1)
		for(yoff=-nnd;yoff<=nnd;yoff+=1)
			for(zoff=-nnd;zoff<=nnd;zoff+=1)
				if((xoff==0) && (yoff==0) && (zoff==0))
					continue
				endif
				nnoffx[j] = xoff
				nnoffy[j] = yoff
				nnoffz[j] = zoff
				nndist[j] = sqrt(xoff^2 + yoff^2 + zoff^2)
				j+=1
			endfor
		endfor
	endfor

	variable txloc, tyloc, tzloc, itteration=0, interfacialpixels=0, starttime, elapsedtime, lasttime, lastpercent=0
	variable /g percentdone
	starttime = ticks/60.15 // seconds when we started
	lasttime = ticks/60.15
	
	string dfolder = getdatafolder(1)
	string/g timeleftstr="?", timeelapsedstr=time2str(ticks/60.15-starttime)
	string normalvecs, centervec
	variable timerref=startmstimer
	make /n=15 /o timerdebug=0
	make /n=15 /t/o timerdebugname=""
	dowindow /k TimingDebug
	display /k=1 /W=(1273,45,1667,350) /n=TimingDebug timerdebug vs timerdebugname as "debug of timing"
	ModifyGraph /w=TimingDebug tkLblRot(bottom)=45,lsize=2,rgb=(0,0,0),useBarStrokeRGB=1
	// loop
	addtologbook(s3d,"Beginning Lamella Creation")
	do
	//nucleation stage
		// create a map of points at least 2 points away from aligned cells (for nucleation)
		duplicate /o empty, tempempty
		imageMorphology /E=501 /i=1 /W=1 /O Erosion tempempty
		// nucleate with some probability
		nucmap[][][] = (enoise(50) > 50-nucprob) && tempempty[p][q][r] ? 255 : 0
		alignx[][][] = nucmap[p][q][r] ? gnoise(1) : alignx[p][q][r] // most of the alignment is in the  and y directions, (not so much in the z)
		aligny[][][] = nucmap[p][q][r] ? gnoise(1) : aligny[p][q][r] 
		alignz[][][] = nucmap[p][q][r] ? gnoise(1/zratio) : alignz[p][q][r]
		alignmag[][][] = nucmap[p][q][r] ? 1 : alignmag[p][q][r] // the magnitute of alignment within crystals is 1
		// normalize alignment to be unity magnitute
		duplicate /o alignx, magnitude
		magnitude[][][] = sqrt(alignx[p][q][r]^2 + aligny[p][q][r]^2 + alignz[p][q][r]^2)
		alignx[][][] = magnitude[p][q][r] > 0.0001 ? alignx[p][q][r]/magnitude[p][q][r] : 0 // all pixels of align are normalized
		aligny[][][] = magnitude[p][q][r] > 0.0001 ? aligny[p][q][r]/magnitude[p][q][r] : 0
		alignz[][][] = magnitude[p][q][r] > 0.0001 ? alignz[p][q][r]/magnitude[p][q][r] : 0
		nucmap[][][] = magnitude[p][q][r] > 0.0001 ? nucmap[p][q][r] : 0
		lcx[][][] = nucmap[p][q][r] ? p : lcx[p][q][r] // if nucleated, store the nucleation site as the center of nucleation
		lcy[][][] = nucmap[p][q][r] ? q : lcy[p][q][r]
		lcz[][][] = nucmap[p][q][r] ? r : lcz[p][q][r]
		empty -= nucmap
	// edge finding
		duplicate /o empty, edge
		redimension /s edge
		Imagefilter /o /n=(nnd*2+1) Gauss3D edge 
		edge = edge < 252 + min(gnoise(.5),2.99) ? 0 : 255
		redimension /u/b edge
		edge = !edge && empty ? 1 : 0
		edge[0][][]   = edge[nxy-1][q][r]==1 && empty? 1 : edge[0][q][r] 
		edge[nxy-1][][] = edge[0][q][r]==1   && empty? 1 : edge[nxy-1][q][r] 
		edge[][0][]   = edge[p][nxy-1][r]==1 && empty? 1 : edge[p][0][r] 
		edge[][nxy-1][] = edge[p][0][r]==1   && empty? 1 : edge[p][nxy-1][r] //wrap edges around
		edgesum = sum(edge) // count the number of cells at an edge of aligned point
		// for all edge points, 
		if(edgesum==0)
			emptysum = mean(empty)
			
			continue
		endif
		tempempty = empty // tempempty is a floating point, whereas empty is a unsigned integer, so gradients doesn't make sense
		imagefilter /O /n=(1+nnd*2) gauss3D tempempty
		differentiate /dim=0  /EP=0 tempempty /D=gradx
		differentiate /dim=1  /EP=0 tempempty /D=grady
		differentiate /dim=2  /EP=0 tempempty /D=gradz
		magnitude[][][] = sqrt(gradx[p][q][r]^2 + grady[p][q][r]^2 + gradz[p][q][r]^2)
		gradx[][][] = magnitude[p][q][r] > 0.005 ? gradx/magnitude[p][q][r] : 0
		grady[][][] = magnitude[p][q][r] > 0.005 ? grady/magnitude[p][q][r] : 0
		gradz[][][] = magnitude[p][q][r] > 0.005 ? gradz/magnitude[p][q][r] : 0 // the gradients will allow each edge location to have an idea of where it is propagating from
		timerdebug[0]+=stopmstimer(timerref)
		timerdebugname[0] = "Nucleation and Setup"
		timerref = startmstimer
		lastedgeloc = -1
	// loop through edge points
		for(j=0;j<edgesum;j+=1)
			//find the next edge location
			findvalue /s=(lastedgeloc+1) /i=1 edge
			lastedgeloc = v_value
			timerdebug[13]+=stopmstimer(timerref)
			timerdebugname[13] = "finding edge"
			timerref = startmstimer
			xloc = mod(lastedgeloc,nxy)
			yloc = mod(lastedgeloc-xloc,nxy^2)/nxy
			zloc = (lastedgeloc-xloc-nxy*yloc)/nxy^2
			if(xloc==-1)
				addtologbook(s3d,"warning! no empty point found")
				break // something is wrong with edge finding - there are non-1s or 0s somewhow
			endif
			timerdebug[9]+=stopmstimer(timerref)
			timerdebugname[9] = "finding locs"
			timerref = startmstimer
			// remove from the edge list
			edge[xloc][yloc][zloc]=0
			// get the local gradinet (g)
			gx = gradx[xloc][yloc][zloc]
			gy = grady[xloc][yloc][zloc]
			gz = gradz[xloc][yloc][zloc] // the direction of edge propagation
			timerdebug[10]+=stopmstimer(timerref)
			timerdebugname[10] = "setting gz"
			timerref = startmstimer
			
			temppnts=0
			tempapnts=0
			xeff=0
			yeff=0
			zeff=0
			// first go through nearest neighbors to check for continuous crystalline propagation, finding the closest neighbor (must be continuous) which has the same orientation, then calculate the propagated thickness
			nnx = xloc + nnoffx  // arrays of nearest neightbor locations
			nny = yloc + nnoffy
			nnz = zloc + nnoffz
			timerdebug[11]+=stopmstimer(timerref)
			timerdebugname[11] = "Nearest Neighbor locs"
			timerref = startmstimer
			nnx += nnx[p]<0 ? nxy : 0 // wrap around x and y
			nny += nny[p]<0 ? nxy : 0
			timerdebug[14]+=stopmstimer(timerref)
			timerdebugname[14] = "Nearest Neighbor wrap 0"
			timerref = startmstimer
			nnx -= nnx[p]>nxy-1 ? nxy : 0
			nny -= nny[p]>nxy-1 ? nxy : 0
			timerdebug[12]+=stopmstimer(timerref)
			timerdebugname[12] = "Nearest Neighbor wrap nxy"
			timerref = startmstimer
			
			
			
			
			nnvalid = (nnz[p]>=0 && nnz[p]<=nz-1) ? 1 : 0 // the nn is baseline valid if it is within the simulation   ///nnx>=0 && nnx<=nxy-1 && nny>=0 && nny<=nxy-1 && // these should be taken account of
			// the nn also needs to be filled to be a valid propagator (processed as either amorphous or crystalline)
			nnvalid*= nnvalid[p] ? (255-empty[nnx[p]][nny[p]][nnz[p]])/255 : 0 //  (empty is 255 if empty, 0 if filled, valid is (255-255)/255 = 0 if empty, (255-0)/255 = 1 if filled)
			timerdebug[2]+=stopmstimer(timerref)
			timerdebugname[2] = "Nearest Neighbor valid"
			timerref = startmstimer
			annx = nnvalid? alignx[nnx][nny][nnz] : nan // alignment parameters of nearest neighbors (will error out if outside the simulation, which is why we need the validity check)
			anny = nnvalid? aligny[nnx][nny][nnz] : nan
			annz = nnvalid? alignz[nnx][nny][nnz] : nan
			timerdebug[3]+=stopmstimer(timerref)
			timerdebugname[3] = "Nearest Neighbor ann"
			timerref = startmstimer
			nnlcx = nnvalid? lcx[nnx][nny][nnz]-nnx : nan
			nnlcy = nnvalid? lcy[nnx][nny][nnz]-nny : nan
			nnlcz = nnvalid? lcz[nnx][nny][nnz]-nnz : nan // vector from pixel to center location of crystal
			timerdebug[4]+=stopmstimer(timerref)
			timerdebugname[4] = "Nearest Neighbor lc"
			timerref = startmstimer
			nnmag = nnvalid? alignmag[nnx][nny][nnz] : nan //(0 - amorphous, 1 - crystalline, nan - non valid)
			nnprop = nnvalid? ((-nnoffx * gx - nnoffy * gy - nnoffz * gz) + nnmag)/nndist : nan
			//  the maximum number is the closest valid pixel directly in line with the propagation, so use that one to propagate
			// if the point isn't valid, (will be nan)
			// can only be crystalline 
			
			timerdebug[5]+=stopmstimer(timerref)
			timerdebugname[5] = "Nearest Neighbor prop"
			timerref = startmstimer
			
			nnprop += gnoise(.2)
			wavestats /m=1/q/z nnprop
			pindex = v_maxRowLoc
			if(pindex<0)
				empty[xloc][yloc][zloc] = 255
				interfacialpixels+=1
				continue
			endif
			amorphprop = (nnmag[pindex]!=1) ? 1 : 0 // is the chosen propogator amorphous?
			if(!testprop(amorphprop,glamspd,nndist[pindex])) // if lamella propogation is favored (or not) this will fail to propogate depending on the propogator's distance (effectively slowing down propogation) and if it is amorphous
				empty[xloc][yloc][zloc]=255
				continue
			endif
			timerdebug[6]+=stopmstimer(timerref)
			timerdebugname[6] = "Propogation 1"
			timerref = startmstimer
			// find the new projected center and thickness of this pixel using this propagator pixel
			centervec = getlocalcenter(nnlcx[pindex],nnlcy[pindex],nnlcz[pindex],xloc,yloc,zloc,nnoffx[pindex],nnoffy[pindex],nnoffz[pindex],annx[pindex],anny[pindex],annz[pindex])
			if(cmpstr("fail",centervec)) // if true, then the calculation is ok //  if false, the best propagator didn't work, so we are in trouble
				lcx[xloc][yloc][zloc] = str2num(stringfromlist(0,centervec,","))
				lcy[xloc][yloc][zloc] = str2num(stringfromlist(1,centervec,","))
				lcz[xloc][yloc][zloc] = str2num(stringfromlist(2,centervec,","))
				thickness[xloc][yloc][zloc] = str2num(stringfromlist(3,centervec,","))
				bx=annx[pindex] // take the alignment of the propagator
				by=anny[pindex]
				bz=annz[pindex] 
				nnaligned =	 nnmag? annx*bx+anny*by+annz*bz : nan // alignment is the same as the propogator
				wavestats /m=1/q/z nnaligned
				localalignment = v_min
				empty[xloc][yloc][zloc] = 0
				if(nnmag[pindex]==1  && localalignment > .8 ) //crystaline and in aligned region (propogating from crystalline region)
					if((thickness[xloc][yloc][zloc] <= (lamthickness + gnoise(lamuncert))/2) ) // we are within an crystalline region and neighbors are all aligned, so just propagate (lamella -> lamella)
						alignmag[xloc][yloc][zloc] = 1
						bx=annx[pindex]+gnoise(XtalAngUncert) // take the alignment of the propagator
						by=anny[pindex]+gnoise(XtalAngUncert)
						bz=annz[pindex]+gnoise(XtalAngUncert/zratio)// add some noise but less than interface
						
					else
						// the crystalline region is too large or not aligned, so we are amorphous  (Lamella -> amorphous)
						alignmag[xloc][yloc][zloc] = 0
						bx=annx[pindex]+gnoise(angularuncert) // take the alignment of the propagator
						by=anny[pindex]+gnoise(angularuncert)
						bz=annz[pindex]+gnoise(angularuncert/zratio) // add some noise to the propagation because it is at an interface
						temp = sqrt(bx^2 + by^2 + bz^2)
						bx/=temp
						by/=temp
						bz/=temp
				
					endif
				elseif(nnmag[pindex]==1) // crystalline, but not aligned region -> this is an interface  (lamella -> fringe)
					alignmag[xloc][yloc][zloc] = 0
					bx=annx[pindex]+gnoise(angularuncert) // take the alignment of the propagator
					by=anny[pindex]+gnoise(angularuncert)
					bz=annz[pindex]+gnoise(angularuncert/zratio) // add some noise to the propagation because it is at an interface
					temp = sqrt(bx^2 + by^2 + bz^2)
					bx/=temp
					by/=temp
					bz/=temp
				else // amorphous //  if we allow this, amorphous regions keep growing
				
					if(enoise(1)>.5 && nndist[pindex]<1.5 && thickness[xloc][yloc][zloc] <= (lamthickness + gnoise(lamuncert))/2 + amorphthickness + gnoise(amorphuncert)) // we are within an amorphous region, so just propagate (amorphous -> amorphous)
						alignmag[xloc][yloc][zloc] = 0
						bx=annx[pindex]+gnoise(angularuncert) // take the alignment of the propagator
						by=anny[pindex]+gnoise(angularuncert)
						bz=annz[pindex]+gnoise(angularuncert/zratio) // add some noise to the propagation because it is amorphous
						temp = sqrt(bx^2 + by^2 + bz^2)
						bx/=temp
						by/=temp
						bz/=temp
					elseif(enoise(1)>.95 && nndist[pindex]<1.5 && thickness[xloc][yloc][zloc] > (lamthickness + gnoise(lamuncert))/2 + amorphthickness + gnoise(amorphuncert))
//						// the amorphous region is too large, so we propagate the alignment in a new crystal (amorphous -> new crystal)
						newthickness = (lamthickness + gnoise(lamuncert))/2
						lcx[xloc][yloc][zloc] = xloc + newthickness * (xloc-lcx[xloc][yloc][zloc])/thickness[xloc][yloc][zloc] // project a new crystal center
						lcy[xloc][yloc][zloc] = yloc + newthickness * (yloc-lcy[xloc][yloc][zloc])/thickness[xloc][yloc][zloc]
						lcz[xloc][yloc][zloc] = zloc + newthickness * (zloc-lcz[xloc][yloc][zloc])/thickness[xloc][yloc][zloc]
						thickness[xloc][yloc][zloc] = newthickness
						alignmag[xloc][yloc][zloc] = 1
						bx=annx[pindex]+gnoise(angularuncert) // take the alignment of the propagator
						by=anny[pindex]+gnoise(angularuncert)
						bz=annz[pindex]+gnoise(angularuncert/zratio) // add some noise to the propagation because it is at an interface
						temp = sqrt(bx^2 + by^2 + bz^2)
						bx/=temp
						by/=temp
						bz/=temp
					else
						empty[xloc][yloc][zloc]=255
						continue
					endif
				endif
				alignx[xloc][yloc][zloc] = bx
				aligny[xloc][yloc][zloc] = by
				alignz[xloc][yloc][zloc] = bz
				timerdebug[7]+=stopmstimer(timerref)
				timerdebugname[7] = "Propogation 2"
				timerref = startmstimer
				//timeelapsedstr=time2str(ticks/60.15-starttime)
				continue
			endif
			// there are no propagators - this is a statistical probability
			empty[xloc][yloc][zloc] = 255
			interfacialpixels+=1
		endfor
		s3d.progressbar =100- 100*(emptysum-stopsum)/(255-stopsum)
		percentdone = s3d.progressbar
		if(lastpercent == percentdone)
			noprogress+=1
		else
			noprogress=0
		endif
		if(noprogress>20) // three itterations with no progress - we are done!
			emptysum = mean(empty)
			emptysum=0
		endif
		
		if(percentdone-lastpercent > .25) //doupdate
			itteration +=1
			timeleftstr = time2str(min((ticks/60.15-starttime)*(100/percentdone - 1),(ticks/60.15-lasttime)*(100/(percentdone - lastpercent) - 1)))
			if(s3d.movie)
				if(itteration==1)
					tempalignmap(alignx,aligny,alignz,alignmag,0,0)
					ModifyGraph rgb=(3,52428,1)
					//appendimage thickness
					//ModifyImage thickness plane=3
				else
					tempalignmap(alignx,aligny,alignz,alignmag,0,1)
				endif
				doupdate
				savepict /p=_PictGallery_ /E=-5 /N=alignmentmap/w=(0,0,800,800) /o as "LamellaFrame"
				addmovieframe /pict=LamellaFrame
				addtologbook(s3d,num2str(round(10*percentdone)/10) + "% complete - " + timeleftstr + " Expected remaining time")
			endif

			
			updateprogress(s3d)
			
			lastpercent = round(10*percentdone)/10
			lasttime = ticks/60.15
			if(interfacialpixels > 50000)
				addtologbook(s3d,"More than 50000 weird pixels tried, something is wrong, stopping")
				return 0
			endif
		endif
		
		emptysum = mean(empty)
		timerdebug[8]+=stopmstimer(timerref)
		timerdebugname[8] = "EndLoop"
		timerref = startmstimer
			
	while(emptysum>stopsum)
	
	make /n=(dimsize(alignx,0)*1.8,dimsize(alignx,1)*1.8,dimsize(alignx,2))/o alignxext = alignx[mod(p+.6*dimsize(alignx,0),dimsize(alignx,0))][mod(q+.6*dimsize(alignx,1),dimsize(alignx,1))][r]
	make /n=(dimsize(aligny,0)*1.8,dimsize(aligny,1)*1.8,dimsize(alignx,2))/o alignyext = aligny[mod(p+.6*dimsize(aligny,0),dimsize(aligny,0))][mod(q+.6*dimsize(aligny,1),dimsize(aligny,1))][r]
	make /n=(dimsize(aligny,0)*1.8,dimsize(aligny,1)*1.8,dimsize(alignx,2))/o alignzext = alignz[mod(p+.6*dimsize(aligny,0),dimsize(aligny,0))][mod(q+.6*dimsize(aligny,1),dimsize(aligny,1))][r]
	make /n=(dimsize(aligny,0)*1.8,dimsize(aligny,1)*1.8,dimsize(alignx,2))/o alignmagext = alignmag[mod(p+.6*dimsize(aligny,0),dimsize(aligny,0))][mod(q+.6*dimsize(aligny,1),dimsize(aligny,1))][r]
	
	duplicate/o alignxext, edgesx
	duplicate/o alignyext, edgesy
	
	duplicate/o edgesx, dotpn
	dotpn = .5*abs(edgesx * edgesx[min(p+1,nxy*1.8)][q] + edgesy * edgesy[min(p+1,nxy*1.8)][q])
	dotpn += .5*abs(edgesx * edgesx[p][min(q+1,nxy*1.8)] + edgesy * edgesy[p][min(q+1,nxy*1.8)])
	if(GBound>2)
		imagefilter /n=(2*floor(GBound/2)+1) /o gauss3D dotpn
	elseif(Gbound<1)
		dotpn=1 // no edges
	endif
	dotpn = dotpn>min(1,max(0,edgethreshold)) ? 1 : 0
	duplicate /o alignmagext, alignmagsave
	alignmagext *= dotpn
	imagefilter /n=3 /o Gauss3D alignmagext
	
	// get the unaligned edges
	duplicate/o alignmagext, unalignmagext
	unalignmagext = 1-alignmagext
	for(j=0;j<dimsize(alignmagext,2);j+=1)
		imageedgedetection /P=(j) /M=1 frei alignmagext
		wave M_ImageEdges
		M_ImageEdges = M_ImageEdges==255 ? 0 :1
		
		if(swapfringe)
			duplicate /free M_imageedges, tempedges
			imagefilter/o /n=(max(13,amorphAlignWid+8)) gauss3D M_imageedges
			imagefilter/o /n=(max(5,amorphAlignWid)) gauss3D tempedges
			unalignmagext[][][j] *=min(1,max(0,abs(amorphalignment)*(M_imageedges[p][q])-20*tempedges[p][q])) // unalignmagext is the alignment within the unaligned portion
			
		
		else
			imagefilter/o /n=(max(3,amorphAlignWid)) gauss3D M_imageedges
			unalignmagext[][][j] *=min(1,abs(amorphalignment)*M_imageedges[p][q]) // unalignmagext is the alignment within the unaligned portion
		endif
	endfor
	duplicate/o alignx, palignx, paligny, palignz // the 90 degree rotated alignment vectors(for edge alignment)
	if(amorphalignment<0)
		make /free/n=(dimsize(alignx,0),dimsize(alignx,1),dimsize(alignx,2)) /t alignnormal
		alignnormal = normvecs(alignx,aligny,alignz,rand=1)
		palignx = str2num(stringfromlist(0,alignnormal[p][q][r],","))
		paligny = str2num(stringfromlist(1,alignnormal[p][q][r],","))
		palignz = str2num(stringfromlist(2,alignnormal[p][q][r],","))
		duplicate /free palignx, palignnorm
		palignnorm = palignx*paligny*palignz*0==0? 1 : 0
		palignx = palignnorm ? palignx : 0
		paligny = palignnorm ? paligny : 0
		palignz = palignnorm ? palignz : 0
	else
		palignx=alignx
		paligny=aligny
		palignz=alignz
	endif
		
	alignmag = alignmagext[p+.4*dimsize(alignmag,0)][q+.4*dimsize(alignmag,1)][r]
	
	duplicate/o alignmag, unalignmag
	unalignmag = unalignmagext[p+.4*dimsize(alignmag,0)][q+.4*dimsize(alignmag,1)][r]
	//if(swapfringe)
	//	unalignmag = 1-unalignmag - alignmag
	//endif
	
	//interlamella alignment (uses amorphalignment and amorphAlignWid)
	
	mat[][][][0] = alignz[q][r][p]*sqrt(alignmag[q][r][p]) // alignment along the beam direction (z in model, x for X-rays)
	mat[][][][1] = alignx[q][r][p]*sqrt(alignmag[q][r][p]) // alignment along the film plane (x in model, y for X-rays)
	mat[][][][2] = aligny[q][r][p]*sqrt(alignmag[q][r][p]) // alignment along the film plane (y in model, z in X-rays)
	
	mat2[][][][0] = palignz[q][r][p]*sqrt(abs(unalignmag[q][r][p]))
	mat2[][][][1] = palignx[q][r][p]*sqrt(abs(unalignmag[q][r][p]))
	mat2[][][][2] = paligny[q][r][p]*sqrt(abs(unalignmag[q][r][p]))
	mat2[][][][3] = 0
	
	mat[][][][3]=1-mat[p][q][r][0]^2 -mat[p][q][r][1]^2 -mat[p][q][r][2]^2 -mat2[p][q][r][0]^2 -mat2[p][q][r][1]^2 -mat2[p][q][r][2]^2
	 // mat is entirely either amorphous or crystalline at this point, we want to put a overall thicknessmap on top of all of this
	
	duplicate /o alignmag, tempdens
	
	duplicate /o dotpn, nonxtal
	
	imagefilter /o /n=11 gauss3d nonxtal
	nonxtal = nonxtal > .5 ? 1 : 0
	imagefilter /o /n=11 gauss3d nonxtal
	nonxtal = nonxtal > .5 ? 1 : 0
	imagefilter /o /n=11 gauss3d nonxtal
	nonxtal = nonxtal > .5 ? 1 : 0
	imagefilter /o /n=21 gauss3d nonxtal // nonxtal is 0 in amorphous regions, 1 in crystalline regions
	
	tempdens = 1-(1- regiondens/100)*(1-nonxtal[p+.4*dimsize(tempdens,0)][q+.4*dimsize(tempdens,1)][r]) // now this is an overall density filter for the film (no change in crystalline regions, some percent change in amorphous regions)
	
	
	mat[][][][0,2] *= sqrt(tempdens[q][r][p])
	mat2[][][][0,2] *= sqrt(tempdens[q][r][p])
	mat[][][][3] *=tempdens[q][r][p]
	 // add in the local difference for the crystalline / amorphous differences
	if(xtaldens>0)
		mat[][][][3] *= 1-(xtaldens/100)
		mat2[][][][0,2] *= sqrt(1-(xtaldens/100))
	else
		mat[][][][0,2] *= sqrt(1-(xtaldens/100))
	endif
	if(Mat2Density > 1)
		mat *= sqrt(1-(Mat2Density-1))
	else
		mat2[][][][0,2] *= sqrt(Mat2Density)
	endif
	
	setdatafolder foldersave
	
	
	make /n=(s3d.thickness,s3d.num,s3d.num,4) /o m1=mat,m2=mat2, m3=0, m4=0
	m2[][][][3] = 0 // material 2 is only the aligned second material
	m1[][][][3] = (1-mat2inamorph) * mat[p][q][r][3] // this is the unaligned portion of material 1
	m3[][][][3] = mat2inamorph * mat[p][q][r][3] // assigning half of the unaligned portion to the second material (setting materail 3 to the same as material 1 and 2 is default)
	make /n=(s3d.thickness,s3d.num,s3d.num) /o density1 = m1[p][q][r][3] + m1[p][q][r][0]^2 +m1[p][q][r][1]^2 +m1[p][q][r][2]^2
	make /n=(s3d.thickness,s3d.num,s3d.num) /o density2 = m2[p][q][r][3] + m2[p][q][r][0]^2 +m2[p][q][r][1]^2 +m2[p][q][r][2]^2 // can make different density for aligned vs unaligned
	make /n=(s3d.thickness,s3d.num,s3d.num) /o density3 = m3[p][q][r][3]
	m4[][][][3] = 1-density1[p][q][r]-density2[p][q][r]-density3[p][q][r] // the remaining space will fill with material 4 (usually this would be vacuum, so this is a film thickness mask)
	dowindow /k alignmentmap

	
	wave s3d.m1=m1, s3d.m2=m2, s3d.m3=m3, s3d.m4 = m4
	wave s3d.density1 = density1
	wave s3d.density2 = density2
	wave s3d.density3 = density3
	duplicate/o m1, m1save
	duplicate/o m2, m2save
	duplicate/o m3, m3save
	duplicate/o m4, m4save
	doalignmentmap("Current")
	if(s3d.movie)
		doupdate
		savepict /p=_PictGallery_ /E=-5 /N=Current_map/w=(0,0,800,800) /o as "LamellaFrame"
		addmovieframe /pict=LamellaFrame
		addmovieframe /pict=LamellaFrame
		addmovieframe /pict=LamellaFrame
		addmovieframe /pict=LamellaFrame
		addmovieframe /pict=LamellaFrame
		addmovieframe /pict=LamellaFrame
		addmovieframe /pict=LamellaFrame
		addmovieframe /pict=LamellaFrame
		addmovieframe /pict=LamellaFrame
		addmovieframe /pict=LamellaFrame
		addmovieframe /pict=LamellaFrame
		addmovieframe /pict=LamellaFrame
		addmovieframe /pict=LamellaFrame
		addmovieframe /pict=LamellaFrame
		addmovieframe /pict=LamellaFrame
		addmovieframe /pict=LamellaFrame
		addmovieframe /pict=LamellaFrame
		addtologbook(s3d,"100 % complete - lamella creation finished")
	endif
end

function testprop(amorph,glamspeed,dist)
	// find the probability of alignment spreading from one voxel to another, based on if it is amorphous
	variable amorph, glamspeed, dist
	if(glamspeed<1)
		amorph = amorph? 0 : 1
		glamspeed = 1/glamspeed
	endif
	if(glamspeed<0)
		glamspeed = 1
	endif
	if(!amorph) // if not amorphous, just propogate
		return 1
	endif
	if(statsnormalcdf(3/glamspeed,dist,.38+1/glamspeed) > enoise(.5)+.5) // 
	
		return 1
	endif
	return 0
end

function /s getlocalcenter(lcx,lcy,lcz,xloc,yloc,zloc,nnoffx,nnoffy,nnoffz,ax,ay,az) // finds the new local center given an old local center and old point
	variable xloc,yloc,zloc // location of new pixel
	variable lcx,lcy,lcz // vector from pixel in question to its local crystal center
	variable nnoffx,nnoffy,nnoffz // vector from new pixel to pixel
	variable ax,ay,az // alignment vector at the new pixel location (can be different than the old pixel alignment)
	// local crystalline center old crystalline center (loc + nn + lc) - projection of nn onto plane (nn.r1 r1 + nn.r2 r2)
	// where r1 and r2 are the unit vectors perpindicular to the NEW pixel's alignment
	// this is in case the new pixel has a different alignment, it's own thickness relative to it's center is correct (ie the crystal is bending, the center line must also bend)
	// for large bend angles, this will cause problems, but ptobably just a thickness that is too big (so it won't be crystalline) which is about right anyways
	
	// old crystalline center (yes, we already have this info, but it's easy enough to calculate)
	variable cx = xloc + nnoffx + lcx
	variable cy = yloc + nnoffy + lcy
	variable cz = zloc + nnoffz + lcz
	string r1r2 = normvecs(ax, ay, az)
	if(!cmpstr(r1r2,"fail"))
		return "fail"
	endif
	variable r1x = str2num(stringfromlist(0,stringfromlist(0,r1r2,";"),","))
	variable r1y = str2num(stringfromlist(1,stringfromlist(0,r1r2,";"),","))
	variable r1z = str2num(stringfromlist(2,stringfromlist(0,r1r2,";"),","))
	variable r2x = str2num(stringfromlist(0,stringfromlist(1,r1r2,";"),","))
	variable r2y = str2num(stringfromlist(1,stringfromlist(1,r1r2,";"),","))
	variable r2z = str2num(stringfromlist(2,stringfromlist(1,r1r2,";"),","))
	
	variable nlcx = cx - ( (nnoffx*r1x + nnoffy*r1y + nnoffz*r1z)*r1x + (nnoffx*r2x + nnoffy*r2y + nnoffz*r2z)*r2x)
	variable nlcy = cy - ( (nnoffx*r1x + nnoffy*r1y + nnoffz*r1z)*r1y + (nnoffx*r2x + nnoffy*r2y + nnoffz*r2z)*r2y)
	variable nlcz = cz - ( (nnoffx*r1x + nnoffy*r1y + nnoffz*r1z)*r1z + (nnoffx*r2x + nnoffy*r2y + nnoffz*r2z)*r2z)
	
	// thickness is the magnitute of new local crystalline center - new pixel
	variable thickness = sqrt((nlcx-xloc)^2 + (nlcy-yloc)^2 + (nlcz-zloc)^2)
	
	
	string outvec
	sprintf outvec, "%g,%g,%g,%g", nlcx,nlcy,nlcz,thickness
	return outvec
end

function /s normvecs(xin, yin, zin,[rand])
	variable xin, yin, zin,rand
	rand = paramisdefault(rand) ? 0 : rand
	variable nx1, nx2, ny1, ny2, nz1, nz2, mag
	string outvecs // the output of the function as a list of two vectors "x,y,z;a,b,c"
	if(rand)
		make /free win={xin, yin, zin}, wo1={enoise(1),enoise(1),enoise(1)}, wo2={-wo1[1],wo1[0],wo1[3]+3}
		wo1 /=norm(wo1)
		wo2 /=norm(wo2)
	else
		make /free win={xin, yin, zin}, wo1={1,0,0}, wo2={0,1,0}
	endif
	cross /z/free /Dest=test win, wo1
	if(waveexists(test))
		if(norm(test)>0)
			wo1 = test / norm(test)
			cross/z/free/Dest=test wo1,win
			if(norm(test)>0)
				wo2 = test / norm(test)
				sprintf outvecs, "%g,%g,%g;%g,%g,%g", wo1[0],wo1[1],wo1[2], wo2[0],wo2[1],wo2[2]
				return outvecs
			endif
		else
			cross /z/free /Dest=test win, wo2
			if(norm(test)>0)
				wo1 = test / norm(test)
				cross/z/free/Dest=test wo1,win
				if(norm(test)>0)
					wo2 = test / norm(test)
					sprintf outvecs, "%g,%g,%g;%g,%g,%g", wo1[0],wo1[1],wo1[2], wo2[0],wo2[1],wo2[2]
					return outvecs
				endif
			endif
		endif
	endif
	//failed
	return "fail"
end


function tempalignmap(xalign, yalign,zalign,alignmag,slicez,update)
	wave xalign, yalign,zalign, alignmag
	variable slicez, update
	variable scale=1
	string foldersave = getdatafolder(1)
	newdatafolder /o/s root:alignmap
	wave /z xloc, yloc
	if(update &&waveexists(xloc)&&waveexists(yloc))
		wave xcomp, ycomp, zcomp, arrowsyay
		xloc = scale*mod(p,dimsize(xalign,0))
		yloc = scale*floor(p/dimsize(xalign,0))
	else
		make/o /n=(dimsize(xalign,0)*dimsize(xalign,1)) xloc = scale*mod(p,dimsize(xalign,0)),yloc = scale*floor(p/dimsize(xalign,0))
		make/o /n=(dimsize(xalign,0)*dimsize(xalign,1)) xcomp,ycomp,zcomp
		make/o /n=(dimsize(xalign,0)*dimsize(xalign,1),2) arrowsyay
		update=0
	endif
	multithread xcomp = xalign[xloc/scale][yloc/scale][slicez] * alignmag[xloc/scale][yloc/scale][slicez]
	multithread ycomp = yalign[xloc/scale][yloc/scale][slicez] * alignmag[xloc/scale][yloc/scale][slicez]
	multithread zcomp = zalign[xloc/scale][yloc/scale][slicez] * alignmag[xloc/scale][yloc/scale][slicez]
	multithread arrowsyay[][1] = atan(ycomp[p]/xcomp[p])
	multithread arrowsyay[][0] = sqrt(xcomp[p]^2 + ycomp[p]^2) * 5000/dimsize(xalign,0)
	//arrowsyay[][0] = arrowsyay[p][0]< 3 &&  sqrt(zcomp[p]^2 + ycomp[p]^2+ xcomp[p]^2)>.5? 3 : arrowsyay[p][0]
	arrowsyay[][0] = arrowsyay[p][1]*0!=0 ? nan : arrowsyay[p][0]
	if(!update)
		dowindow /k alignmentmap
		Display /k=1/n=alignmentmap/W=(-1861,-486,-457,918)/K=1  yloc vs xloc as "Slice of Alignment"
		ModifyGraph/w=alignmentmap mode=3,gfSize=12
		ModifyGraph/w=alignmentmap rgb=(0,0,0)
		ModifyGraph/w=alignmentmap msize=0.5,marker=42,tlOffset=-5
		ModifyGraph/w=alignmentmap mrkThick=0.1,margin(left)=26,margin(bottom)=26
		ModifyGraph/w=alignmentmap arrowMarker(yloc)={arrowsyay,.75,0,0,1}
		ModifyGraph/w=alignmentmap mirror=0,margin(top)=1,margin(right)=1
		ModifyGraph/w=alignmentmap axOffset(left)=-4.46154,axOffset(bottom)=-0.962963,height={Plan,1,left,bottom}
		ModifyGraph /w=alignmentmap tick=2,mirror=1,axisOnTop=1,standoff=0,tkLblRot(left)=90
		Label /w=alignmentmap left "Size Scale [px]";DelayUpdate
		Label /w=alignmentmap bottom "Size Scale [px]"
	endif
	setdatafolder foldersave
end


function /s variables_spheresln()
	string variables
//	return "Interpenetration [pixels],SetVariable,2;Minimum Seperation,SetVariable,.1;Number of Particles (Max),SetVariable,500;PolyDispursity (sigma of radiuses),setvariable,5;Maximum Radius,SetVariable,8;Noise,SetVariable,0;"
	return "Interpenetration [pixels],SetVariable,2;Minimum Seperation,SetVariable,.1;Number of Particles (Max),SetVariable,500;PolyDispursity (sigma of radiuses),setvariable,5;Minimum Radius,SetVariable,1;Maximum Radius,SetVariable,8;Volume Fraction (<1),SetVariable,.65;Noise,SetVariable,0;"
end
function model3D_Spheresln(s3d)
	//Creates a spherical system, with two components, aligned 
	// utilizes extra parameters, number of particles, polydispursity (a sigma of radiuses) and populates the 3D space with 
		//aligned particles
	// material 1 is insize the spheres, material 2 is outside (although the concentrations are determined by m2volfrac
	//paramstring = 
		//MORPHOLOGY parameters
		//0 interpenetration of spheres (roughness smaller than wavelength)
		//1 minimum seperation (when placing a new sphere into the system, it needs to be this far away from any others)
			//	this can be a fraction of the radius (if less than 1) or a flat number of pixels (if greater than 1)
		//2 number of particles // the maximum number of particles . the program will stop when it cannot fit anymore
		//3 polydispursity (sigma of radiuses)
		//4 thickness (this can be different than size)
		//5 noise (percentage of vacuum - density variations)
	//num = resolution (in each of 3 dimensions)
	//size = average size of sphere (radius)
	
	struct ThreeDSystem &s3d
	if(itemsinlist(s3d.paramstring,",")<5)
		return -1
	endif
	newdatafolder /o/s CreatingSpheres
	variable interpenetration = 	str2num(stringfromlist( 0 ,s3d.paramstring,","))
	variable minsep = 			str2num(stringfromlist( 1 ,s3d.paramstring,","))
	variable particlesnum = 	str2num(stringfromlist( 2 ,s3d.paramstring,","))
	variable pd = 				str2num(stringfromlist( 3 ,s3d.paramstring,","))
	variable thickness = 		s3d.thickness
	variable minsize = 		str2num(stringfromlist( 4 ,s3d.paramstring,","))
	variable maxsize = 		str2num(stringfromlist( 5 ,s3d.paramstring,","))
	variable volfrac = 			str2num(stringfromlist( 6 ,s3d.paramstring,","))
	variable noise = 			str2num(stringfromlist( 7 ,s3d.paramstring,","))
	
	make /o /n=(thickness,s3d.num,s3d.num) mat=1,xwave, ywave, zwave
	if(minsep<1)
		make/B/U /o /n=(thickness,s3d.num,s3d.num,maxsize) exmat =  (p <= t/(1+minsep)) || (q <= t/(1+minsep)) || (r <= t/(1+minsep) ) || (p >= thickness-t/(1+minsep)) || (q >= s3d.num-t/(1+minsep)) || (r >= s3d.num-t/(1+minsep)) ? 0 : 1
	else
		make/B/U /o /n=(thickness,s3d.num,s3d.num,maxsize) exmat =  (p <= t-minsep) || (q <= t-minsep) || (r <= t-minsep) || (p >= thickness-t+minsep) || (q >= s3d.num-t+minsep) || (r >= s3d.num-t+minsep) ? 0 : 1
	endif
	make/B/U /o /n=(thickness,s3d.num,s3d.num) tempwave
	if(s3d.movie)
		Execute("Spheres3Ddisp(" +num2str(s3d.num)+", \""+getwavesdatafolder(mat,2)+"\")")
		Execute("exportgizmo wave=\"testimage\"   ;Spinoidal3DLayout();Spinoidal3DImage(\""+getdatafolder(1)+"testimage\")")
	endif
	setscale /i x, -thickness/2, thickness/2, mat, exmat,xwave, ywave, zwave
	setscale /i y, -s3d.num/2, s3d.num/2, mat, exmat,xwave, ywave, zwave
	setscale /i z, -s3d.num/2, s3d.num/2, mat, exmat,xwave, ywave, zwave
	xwave = x
	ywave = y
	zwave = z
	redimension /n=(thickness*s3d.num*s3d.num) xwave, ywave, zwave
	variable testcx,testcy,testcz,i,radius, orad, cx, cy, cz, failed, fnum =0, xmn,xmx,ymn,ymx,zmn,zmx, loc, qfnum=0
	fnum=0
	if(minsep<1)
		orad = s3d.size*(1+minsep/2)
	else
		orad = s3d.size + minsep/2
	endif
	wave radiuses = getradiusdistributionln(orad,s3d.num,thickness,pd,minsep,minsize,maxsize,volfrac)
	if(minsep<1)
		radiuses /= (1+minsep/2)
	else
		radiuses -= minsep/2
	endif
	sort /R radiuses, radiuses
	for(i=0;i<numpnts(radiuses);i+=1)
		radius = radiuses[i]
		radius =radius < 1 ? 1 : radius
		if(minsep<1)
			orad = radius*(1+minsep/2)
		else
			orad = radius + minsep/2
		endif
		qfnum=0
		do
			failed = 0
			testcx = enoise((thickness-orad)/2)
			testcy = enoise((s3d.num-orad)/2)
			testcz = enoise((s3d.num-orad)/2)
			if(exmat(testcx)(testcy)(testcz)[ceil(orad)] ==0)
				failed=1
				qfnum+=1
			else
				cx = testcx
				cy = testcy
				cz = testcz
			endif
		while(failed==1 && qfnum<300)
		
		if(failed==1)
			//redimension /n=(thickness,s3d.num,s3d.num) tempwave
			multithread tempwave[][][] = exmat[p][q][r][ceil(orad)]
	  		
	  		if(wavemax(tempwave)<1)
			//there are no possible locations for this radius, find another
				print "cannot fit " +num2str(radius) + "sized particle"
				fnum +=1
				if(fnum>10)
					print "warning : failed - only " +num2str(i-fnum) +" particles created"		
					break
				endif
				continue // get out of the loop
			endif
			// randomly pick a pixel that is good for the center
			//redimension /n=(thickness*s3d.num*s3d.num) tempwave
			integrate tempwave /D=intwave
			loc = binarysearch(intwave, enoise(wavemax(intwave)/2)+wavemax(intwave)/2)
			cx = xwave[loc]
			cy = ywave[loc]
			cz = zwave[loc]
		endif
		
		// subtract out this sphere from the matrix  // matrix starts at 1s, within this sphere, multiply this by 0, outside multiply by 1
		multithread mat*= (x-cx)^2 + (y-cy)^2 + (z-cz)^2 < radius^2 ? 0 : 1 
		multithread exmat*= (x-cx)^2 + (y-cy)^2 + (z-cz)^2 <= (orad+t)^2 ? 0 : 1 
		print num2str(i) + " Particles created out of "+ num2str(numpnts(radiuses)) + " failed quick spheres - " + num2str(qfnum)
		if(s3d.movie)
			execute("ModifyGizmo /n=Spheres3D update=2")
			doupdate
			Execute "exportgizmo wave=\"testimage\"   "
			//TextBox/w=Spinoidal3DLayout/C/N=text0/A=LT/X=0.00/Y=0.00 "\Z32" + time2str2(ttot)
			doupdate
			savepict /p=_PictGallery_ /E=-5 /N=Spinoidal3DLayout /o as "Frame3D"
			addmovieframe /pict=Frame3D
		endif
	endfor
	setdatafolder ::
	imagefilter /n=(interpenetration)/o gauss3d mat
	duplicate /o mat,s3d.density1 // this returns the density matrix of material 1 (the matrix) for alignment etc later on
end
function /wave getradiusdistributionln(size,length,thickness, pd,minsep,minsize,maxsize,volfrac)
	variable size, length, thickness, pd,minsep, maxsize, minsize, volfrac
	variable vol = length^2 * thickness, orad
	
	if(minsep<1)
		orad = size*(1+minsep/2)
		//maxsize /=(1+minsep/2)
	else
		orad = size + minsep/2
		//maxsize -=minsep/2
	endif
	
	make/o /n=(10000) /o prob, cumprob
	make/o/n=100000 radiuses=0
	setscale /i x, minsize, maxsize, prob, cumprob
	
	prob = (1/(x*ln(pd)))*exp(-(ln(x)-ln(size))^2 /(2*ln(pd)^2))
	prob[0]=0
	integrate prob /D=cumprob
	variable i=-1, totalvol=0, radius
	do
		i+=1
		radius = BinarySearchInterp(cumprob, enoise(wavemax(cumprob)/2)+wavemax(cumprob)/2)*(maxsize-minsize)/numpnts(cumprob) + minsize
		radiuses[i] = radius
		if(minsep<1)
			orad = radius*(1+minsep/2)
		else
			orad = radius + minsep/2
		endif
		totalvol += (4/3)*pi*orad^3
	while(totalvol/vol < volfrac)
	redimension /n=(i-1) radiuses
	return radiuses
end
function /s variables_Cylsln()
	string variables
	return "Interpenetration [pixels],SetVariable,2;Minimum Seperation,SetVariable,.1;Number of Particles (Max),SetVariable,500;PolyDispursity (sigma of radiuses),setvariable,5;Minimum Radius,SetVariable,1;Maximum Radius,SetVariable,8;Volume Fraction (<1),SetVariable,.65;Noise,SetVariable,0;"
end
function model3D_Cylsln(s3d)
	//Creates a spherical system, with two components, aligned 
	// utilizes extra parameters, number of particles, polydispursity (a sigma of radiuses) and populates the 3D space with 
		//aligned particles
	// material 1 is insize the spheres, material 2 is outside (although the concentrations are determined by m2volfrac
	//paramstring = 
		//MORPHOLOGY parameters
		//0 interpenetration of spheres (roughness smaller than wavelength)
		//1 minimum seperation (when placing a new sphere into the system, it needs to be this far away from any others)
			//	this can be a fraction of the radius (if less than 1) or a flat number of pixels (if greater than 1)
		//2 number of particles // the maximum number of particles . the program will stop when it cannot fit anymore
		//3 polydispursity (sigma of radiuses)
		//4 thickness (this can be different than size)
		//5 noise (percentage of vacuum - density variations)
	//num = resolution (in each of 3 dimensions)
	//size = average size of sphere (radius)
	
	struct ThreeDSystem &s3d
	if(itemsinlist(s3d.paramstring,",")<5)
		return -1
	endif
	newdatafolder /o/s CreatingSpheres
	variable interpenetration = 	str2num(stringfromlist( 0 ,s3d.paramstring,","))
	variable minsep = 			str2num(stringfromlist( 1 ,s3d.paramstring,","))
	variable particlesnum = 	str2num(stringfromlist( 2 ,s3d.paramstring,","))
	variable pd = 				str2num(stringfromlist( 3 ,s3d.paramstring,","))
	variable thickness = 		s3d.thickness
	variable minsize = 		str2num(stringfromlist( 4 ,s3d.paramstring,","))
	variable maxsize = 		str2num(stringfromlist( 5 ,s3d.paramstring,","))
	variable volfrac = 			str2num(stringfromlist( 6 ,s3d.paramstring,","))
	variable noise = 			str2num(stringfromlist( 7 ,s3d.paramstring,","))
	

	
	make /o /n=(s3d.num,s3d.num) mat=1, ywave, zwave
	if(minsep<1)
		make/B/U /o /n=(s3d.num,s3d.num,30) exmat =  (p <= r/(1+minsep)) || (q <= r/(1+minsep) ) || (p >= s3d.num-r/(1+minsep)) || (q >= s3d.num-r/(1+minsep)) ? 0 : 1
	else
		make/B/U /o /n=(s3d.num,s3d.num,30) exmat =   (p <= r-minsep) || (q <= r-minsep) || (p >= s3d.num-r+minsep) || (q >= s3d.num-r+minsep) ? 0 : 1
	endif
	make/B/U /o /n=(s3d.num,s3d.num) tempwave
	if(s3d.movie)
		//Execute("Spheres3Ddisp(" +num2str(s3d.num)+", \""+getwavesdatafolder(mat,2)+"\")")
		//Execute("exportgizmo wave=\"testimage\"   ;Spinoidal3DLayout();Spinoidal3DImage(\""+getdatafolder(1)+"testimage\")")
		newimage /k=1 /n=cylimage mat
	endif

	setscale /i y, -s3d.num/2, s3d.num/2, mat, exmat, ywave, zwave
	setscale /i z, -s3d.num/2, s3d.num/2, mat, exmat, ywave, zwave
	
	ywave = x// 2D y-> x and z->y
	zwave = y
	redimension /n=(s3d.num*s3d.num) ywave, zwave
	variable i,radius, orad, cx, cy, cz, failed, fnum =0, xmn,xmx,ymn,ymx,zmn,zmx, loc
	fnum=0
	if(minsep<1)
		orad = s3d.size*(1+minsep/2)
	else
		orad = s3d.size + minsep/2
	endif
	wave radiuses = get2dradiusdistributionln(orad,s3d.num,pd,minsep,minsize,maxsize,volfrac)
	if(minsep<1)
		radiuses /= (1+minsep/2)
	else
		radiuses -= minsep/2
	endif
	sort /R radiuses, radiuses
	for(i=0;i<numpnts(radiuses);i+=1)
		radius = radiuses[i]
		radius =radius < 1 ? 1 : radius
		if(minsep<1)
			orad = radius*(1+minsep/2)
		else
			orad = radius + minsep/2
		endif
		redimension /n=(s3d.num,s3d.num) tempwave
		multithread tempwave[][] = exmat[p][q][ceil(orad)]
  		
  		if(wavemax(tempwave)<1)
		//there are no possible locations for this radius, find another
			print "cannot fit " +num2str(radius) + "sized particle"
			fnum +=1
			if(fnum>10)
				print "warning : failed - only " +num2str(i-fnum) +" particles created"		
				break
			endif
			continue // get out of the loop
		endif
		// randomly pick a pixel that is good for the center
		redimension /n=(s3d.num*s3d.num) tempwave
		integrate tempwave /D=intwave
		loc = binarysearch(intwave, enoise(wavemax(intwave)/2)+wavemax(intwave)/2)
		cx = ywave[loc]
		cy = zwave[loc]
		
		// subtract out this sphere from the matrix  // matrix starts at 1s, within this sphere, multiply this by 0, outside multiply by 1
		multithread mat*= (x-cx)^2 + (y-cy)^2 < radius^2 ? 0 : 1 
		multithread exmat*= (x-cx)^2 + (y-cy)^2 <= (orad+r)^2 ? 0 : 1 
		if(s3d.movie)
			//execute("ModifyGizmo /n=Spheres3D update=2")
			doupdate
			//Execute "exportgizmo wave=\"testimage\"   "
			//TextBox/w=Spinoidal3DLayout/C/N=text0/A=LT/X=0.00/Y=0.00 "\Z32" + time2str2(ttot)
			//doupdate
			//savepict /p=_PictGallery_ /E=-5 /N=Spinoidal3DLayout /o as "Frame3D"
			dowindow /F cylimage
			addmovieframe // /pict=Frame3D
		endif
	endfor
	setdatafolder ::
	make /o /n=(thickness,s3d.num,s3d.num) s3d.density1=mat[q][r]
	setscale /i x, -thickness/2, thickness/2, s3d.density1
	setscale /i y, -s3d.num/2, s3d.num/2, s3d.density1
	setscale /i z, -s3d.num/2, s3d.num/2, s3d.density1 // this returns the density matrix of material 1 (the matrix) for alignment etc later on
	s3d.density1 +=gnoise(noise)
	s3d.density1 = abs(s3d.density1)
	imagefilter /n=(interpenetration)/o gauss3d s3d.density1
end
function /wave get2dradiusdistributionln(size,length,pd,minsep,minsize,maxsize,volfrac)
	variable size, length, pd,minsep, maxsize,volfrac,minsize
	variable vol = length^2, orad
	
	if(minsep<1)
		orad = size*(1+minsep/2)
		//maxsize /=(1+minsep/2)
	else
		orad = size + minsep/2
		//maxsize -=minsep/2
	endif
	
	make/o /n=(10000) /o prob, cumprob
	make/o/n=100000 radiuses=0
	setscale /i x, minsize, maxsize, prob, cumprob
	
	prob = (1/(x*ln(pd)))*exp(-(ln(x)-ln(size))^2 /(2*ln(pd)^2))
	prob[0]=0
	integrate prob /D=cumprob
	variable i=-1, totalvol=0, radius
	do
		i+=1
		radius = BinarySearchInterp(cumprob, enoise(wavemax(cumprob)/2)+wavemax(cumprob)/2)*(maxsize-minsize)/numpnts(cumprob) + minsize 
		radiuses[i] = radius
		if(minsep<1)
			orad = radius*(1+minsep/2)
		else
			orad = radius + minsep/2
		endif
		totalvol += pi*orad^2
	while(totalvol/vol < volfrac)
	redimension /n=(i-1) radiuses
	return radiuses
end
function /wave getradiusdistribution2ln(size,size2,length,thickness, pd,pd2)
	variable size, length, thickness, pd, size2, pd2
	variable vol = length^2 * thickness
	
	make/o /n=(10000) /o prob, cumprob
	make/o/n=100000 radiuses=0
	setscale /i x, 0, thickness/2, prob, cumprob
	
	prob = (1/(x*ln(pd)))*exp(-(ln(x)-ln(size))^2 /(2*ln(pd)^2)) + (1/(x*ln(pd2)))*exp(-(ln(x)-ln(size2))^2 /(2*ln(pd2)^2))
	prob[0]=0
	integrate prob /D=cumprob
	variable i=-1, totalvol=0, radius
	do
		i+=1
		radius = BinarySearchInterp(cumprob, enoise(wavemax(cumprob)/2)+wavemax(cumprob)/2)*thickness/(2*10000) 
		radiuses[i] = radius
		totalvol += (4/3)*pi*radius^3
	while(totalvol/vol < .4)
	redimension /n=(i-1) radiuses
	return radiuses
end
function /s variables_spheres2ln()
	string variables
	variables += "Minimum Seperation,SetVariable,.1;"
	variables += "PolyDispursity (sigma of radius),setvariable,1.2;"
	variables += "Noise,SetVariable,0;"
	variables += "Radius 2,SetVariable,12;"
	variables += "PolyDispursity (sigma of radius 2),setvariable,1.2;"
	variables  =  "Interpenetration [pixels],SetVariable,2;"
	return variables
end
function model3D_Spheres2ln(s3d)
	//Creates a spherical system, with two components, aligned 
	// utilizes extra parameters, number of particles, polydispursity (a sigma of radiuses) and populates the 3D space with 
		//aligned particles
	// material 1 is insize the spheres, material 2 is outside (although the concentrations are determined by m2volfrac
	//paramstring = 
		//MORPHOLOGY parameters
		//0 interpenetration of spheres (roughness smaller than wavelength)
		//1 minimum seperation (when placing a new sphere into the system, it needs to be this far away from any others)
			//	this can be a fraction of the radius (if less than 1) or a flat number of pixels (if greater than 1)
		//2 number of particles // the maximum number of particles . the program will stop when it cannot fit anymore
		//3 polydispursity (sigma of radiuses)
		//4 thickness (this can be different than size)
		//5 noise (percentage of vacuum - density variations)
	//num = resolution (in each of 3 dimensions)
	//size = average size of sphere (radius)
	
	struct ThreeDSystem &s3d
	if(itemsinlist(s3d.paramstring,",")<5)
		return -1
	endif
	newdatafolder /o/s CreatingSpheres
	variable minsep = 		str2num(stringfromlist( 0 ,s3d.paramstring,","))
	variable pd = 				str2num(stringfromlist( 1 ,s3d.paramstring,","))
	variable thickness = 		s3d.thickness
	variable noise = 			str2num(stringfromlist( 2 ,s3d.paramstring,","))
	variable rad2 = 			str2num(stringfromlist( 3 ,s3d.paramstring,","))
	variable pd2 = 			str2num(stringfromlist( 4 ,s3d.paramstring,","))
	variable interpenetration = 	str2num(stringfromlist( 5 ,s3d.paramstring,","))
	

	
	make /o /n=(thickness,s3d.num,s3d.num) mat=1,xwave, ywave, zwave
	
	make/B/U /o /n=(thickness,s3d.num,s3d.num,30) exmat= (p <= t) || (q <= t) || (r <= t) || (p >= thickness-t) || (q >= s3d.num-t) || (r >= s3d.num-t) ? 0 : 1
	make/B/U /o /n=(thickness,s3d.num,s3d.num) tempwave
	if(s3d.movie)
		Execute("Spheres3Ddisp(" +num2str(s3d.num)+", \""+getwavesdatafolder(mat,2)+"\")")
		Execute("exportgizmo wave=\"testimage\"   ;Spinoidal3DLayout();Spinoidal3DImage(\""+getdatafolder(1)+"testimage\")")
	endif
	setscale /i x, -thickness/2, thickness/2, mat, exmat,xwave, ywave, zwave
	setscale /i y, -s3d.num/2, s3d.num/2, mat, exmat,xwave, ywave, zwave
	setscale /i z, -s3d.num/2, s3d.num/2, mat, exmat,xwave, ywave, zwave
	xwave = x
	ywave = y
	zwave = z
	redimension /n=(thickness*s3d.num*s3d.num) xwave, ywave, zwave
	variable i,radius, orad, cx, cy, cz, failed, fnum =0, xmn,xmx,ymn,ymx,zmn,zmx, loc
	fnum=0

	wave radiuses = getradiusdistribution2ln(s3d.size, rad2,s3d.num,thickness,pd, pd2)
	
	sort /R radiuses, radiuses
	for(i=0;i<numpnts(radiuses);i+=1)
		radius = radiuses[i]
		radius =radius < 1 ? 1 : radius
		if(minsep<1)
			orad = radius*(1+minsep/2)
		else
			orad = radius + minsep/2
		endif
		redimension /n=(thickness,s3d.num,s3d.num) tempwave
		multithread tempwave[][][] = exmat[p][q][r][ceil(orad)]
  		
  		if(wavemax(tempwave)<1)
		//there are no possible locations for this radius, find another
			fnum +=1
			if(fnum>10)
				print "warning : failed to find unoccupied location - only " +num2str(i-fnum) +" particles created"		
				break
			endif
			continue // get out of the loop
		endif
		// randomly pick a pixel that is good for the center
		redimension /n=(thickness*s3d.num*s3d.num) tempwave
		integrate tempwave /D=intwave
		loc = binarysearch(intwave, enoise(wavemax(intwave)/2)+wavemax(intwave)/2)
		cx = xwave[loc]
		cy = ywave[loc]
		cz = zwave[loc]
		
		// subtract out this sphere from the matrix  // matrix starts at 1s, within this sphere, multiply this by 0, outside multiply by 1
		multithread mat*= (x-cx)^2 + (y-cy)^2 + (z-cz)^2 < radius^2 ? 0 : 1 
		multithread exmat*= (x-cx)^2 + (y-cy)^2 + (z-cz)^2 <= (orad+t)^2 ? 0 : 1 
		if(s3d.movie)
			execute("ModifyGizmo /n=Spheres3D update=2")
			doupdate
			Execute "exportgizmo wave=\"testimage\"   "
			//TextBox/w=Spinoidal3DLayout/C/N=text0/A=LT/X=0.00/Y=0.00 "\Z32" + time2str2(ttot)
			doupdate
			savepict /p=_PictGallery_ /E=-5 /N=Spinoidal3DLayout /o as "Frame3D"
			addmovieframe /pict=Frame3D
		endif
	endfor
	setdatafolder ::
	imagefilter /n=(interpenetration)/o gauss3d mat
	duplicate /o mat,s3d.density1 // this returns the density matrix of material 1 (the matrix) for alignment etc later on
end

function calcchordanalysis()
	setdatafolder root:
	newdatafolder /s/o spheres
	duplicate/o root:packages:ScatterSim3D:testln2:CreatingSpheres:mat root:spheres:a
	wave a
	fft /mags /dest=fa a
	make /o/n=(256,256) ap = fa[0][p][q]
	newimage ap
	ModifyImage ap log=1
	setscale/i x, -1/2,1/2 , ap
	setscale/i y, -1/2, 1/2, ap
	wave radialintensity = radialintegratew(ap,0,360,"radialintensity")
	
	radialintensity *=x^2
	fft /PAD=256 /real /dest=corr radialintensity
	display corr
	smooth 2, corr
	ModifyGraph log(left)=0
	ModifyGraph log=0
	differentiate corr /D=d1corr
	display d1corr
	ModifyGraph log(left)=0
	differentiate d1corr /D=d2corr
	display d2corr
end

function make3darrows(size,thick,ox,oy,oz,amat,dmat)
	wave amat,dmat
	variable size,thick, ox, oy, oz
	
	make /n=(size^2*thick,3)/o arrows3dloc, arrows3dsize
	make /n=(size^2*thick,4)/o arrows3drot
	arrows3dloc[][0]=mod(p,thick) +ox
	arrows3dloc[][1]=mod(floor(p/thick),size) + oy
	arrows3dloc[][2]=floor(p/(size*thick)) + oz
	
	arrows3dsize[][] = sqrt(amat[arrows3dloc[p][0]][arrows3dloc[p][1]][arrows3dloc[p][2]][0]^2 + amat[arrows3dloc[p][0]][arrows3dloc[p][1]][arrows3dloc[p][2]][1]^2 + amat[arrows3dloc[p][0]][arrows3dloc[p][1]][arrows3dloc[p][2]][2]^2) /5
	arrows3drot[][] = cang(0,0,1,amat[arrows3dloc[p][0]][arrows3dloc[p][1]][arrows3dloc[p][2]][0],amat[arrows3dloc[p][0]][arrows3dloc[p][1]][arrows3dloc[p][2]][1],amat[arrows3dloc[p][0]][arrows3dloc[p][1]][arrows3dloc[p][2]][2],q)
	
	make/o/n=(thick,size,size) arrowsurface
	arrowsurface = dmat[ox+p][oy+q][oz+z]
	setscale /p x, ox, 1, arrowsurface
	setscale /p y, oy, 1, arrowsurface
	setscale /p z, oz, 1, arrowsurface
	string datafolder = getdatafolder(1)
	execute("Alignment3D(\""+datafolder+"\")")

end

function cang(v1,v2,v3,v4,v5,v6,elem)
	variable v1,v2,v3,v4,v5,v6,elem
	string returns
	make /n=3 /o vec1 = {v1,v2,v3}, vec2 = {v4,v5,v6}
	
	variable mag = sqrt(vec1[0]^2 + vec1[1]^2 + vec1[2]^2)
	vec1/=mag
	mag = sqrt(vec2[0]^2 + vec2[1]^2 + vec2[2]^2)
	vec2/=mag
	make /n=4 /o w_cang
	cross /T/Z vec2, vec1
	wave w_cross
	mag = sqrt(w_cross[0]^2 + w_cross[1]^2 + w_cross[2]^2)
	if(mag>0)
		w_cross /=mag
		mag = asin(mag)*180/pi
		if(v6>0)
			w_cang[0]= 180-mag
		else
			w_cang[0]= mag
		endif
		w_cang[1,]=w_cross[p-1]
		//print w_cang
		return w_cang[elem]
	else
		w_cang = {0,1,0,0}
		//print "invalid result"
		return w_cang[elem]
	endif
end
Window Alignment3D(loc) : GizmoPlot
	string loc
	PauseUpdate; Silent 1	// Building Gizmo 6 window...

	// Do nothing if the Gizmo XOP is not available.
	if(exists("NewGizmo")!=4)
		DoAlert 0, "Gizmo XOP must be installed"
		return
	endif

	dowindow /F Alignment3D
	if(v_flag)
		return
	endif


	NewGizmo/K=1/N=Alignment3D/T="Alignment 3D View"
	ModifyGizmo startRecMacro
	MoveWindow 60,46.25,934.5,686.75 
	ModifyGizmo opName=scale0, operation=scale,data={0.2,1,1}
	AppendToGizmo Scatter=$(loc+"arrows3dloc"),name=scatter0
	ModifyGizmo ModifyObject=scatter0 property={ scatterColorType,0}
	ModifyGizmo ModifyObject=scatter0 property={ markerType,0}
	ModifyGizmo ModifyObject=scatter0 property={ sizeType,1}
	ModifyGizmo ModifyObject=scatter0 property={ rotationType,1}
	ModifyGizmo ModifyObject=scatter0 property={ Shape,10}
	ModifyGizmo ModifyObject=scatter0 property={ size,1}
	ModifyGizmo ModifyObject=scatter0 property={ color,0,0,0.996109,0.5}
	ModifyGizmo ModifyObject=scatter0 property={ rotationWave,$(loc+"arrows3drot")}
	ModifyGizmo ModifyObject=scatter0 property={ sizeWave,$(loc+"arrows3dsize")}
	AppendToGizmo isoSurface=$(loc+"arrowsurface"),name=isoSurface0
	ModifyGizmo ModifyObject=isoSurface0 property={ surfaceColorType,1}
	ModifyGizmo ModifyObject=isoSurface0 property={ lineColorType,1}
	ModifyGizmo ModifyObject=isoSurface0 property={ lineWidthType,0}
	ModifyGizmo ModifyObject=isoSurface0 property={ fillMode,2}
	ModifyGizmo ModifyObject=isoSurface0 property={ lineWidth,1}
	ModifyGizmo ModifyObject=isoSurface0 property={ isoValue,0.75}
	ModifyGizmo ModifyObject=isoSurface0 property={ frontColor,1,0,0,0}
	ModifyGizmo ModifyObject=isoSurface0 property={ backColor,0,0.597665,0,0}
	ModifyGizmo ModifyObject=isoSurface0 property={ lineColor,0,0,0,1}
	ModifyGizmo modifyObject=isoSurface0 property={calcNormals,1}
	AppendToGizmo light=Directional,name=light0
	ModifyGizmo light=light0 property={ liveUpdate,1}
	ModifyGizmo light=light0 property={ position,-0.607562,-0.336777,-0.719340,0.000000}
	ModifyGizmo light=light0 property={ direction,-0.607562,-0.336777,-0.719340}
	ModifyGizmo light=light0 property={ specular,1.000000,1.000000,1.000000,1.000000}
	AppendToGizmo light=Directional,name=light1
	ModifyGizmo light=light1 property={ liveUpdate,1}
	ModifyGizmo light=light1 property={ position,-0.332028,-0.651641,-0.681998,0.000000}
	ModifyGizmo light=light1 property={ direction,-0.332028,-0.651641,-0.681998}
	ModifyGizmo light=light1 property={ specular,1.000000,1.000000,1.000000,1.000000}
	AppendToGizmo light=Directional,name=light2
	ModifyGizmo light=light2 property={ liveUpdate,1}
	ModifyGizmo light=light2 property={ position,0.584553,0.297844,-0.754709,0.000000}
	ModifyGizmo light=light2 property={ direction,0.584553,0.297844,-0.754709}
	ModifyGizmo light=light2 property={ specular,1.000000,1.000000,1.000000,1.000000}
	ModifyGizmo setDisplayList=0, opName=ortho0, operation=ortho, data={-1.48235,1.48235,-1.48235,1.48235,-2,2}
	ModifyGizmo setDisplayList=1, object=light2
	ModifyGizmo setDisplayList=2, object=scatter0
	ModifyGizmo setDisplayList=3, object=light1
	ModifyGizmo setDisplayList=4, object=light0
	ModifyGizmo setDisplayList=5, object=isoSurface0
	ModifyGizmo SETQUATERNION={-0.102261,-0.307889,0.638372,0.698024}
	ModifyGizmo autoscaling=1
	ModifyGizmo currentGroupObject=""
	ModifyGizmo compile

	ModifyGizmo showInfo
	ModifyGizmo infoWindow={1129,172,1629,426}
	ModifyGizmo bringToFront
	ModifyGizmo userString={wmgizmo_df,"DF_Gizmo00"}
	ModifyGizmo endRecMacro
End


function azimuthplot(en, minq, maxq)
	variable minq, maxq, en
	nvar enmin = ::startxrayenergy
	nvar enmax = ::endxrayenergy
	wave scatter3dsave
	setscale /i z, enmin, enmax, scatter3dsave
	duplicate /o /r=()()(en) scatter3dsave, tempwave
	duplicate /o tempwave, distwave, anglewave
	distwave = sqrt(x^2 + y^2)
	anglewave = y>0 ? atan(x/y) : pi+ atan(x/y)
	anglewave *= 180/pi
	anglewave +=90
	redimension /n=(numpnts(distwave)) distwave, anglewave, tempwave
	distwave = distwave <minq || distwave> maxq ? 0 : 1
	tempwave *=distwave
	anglewave *= distwave
	make/o /n=360 azhist, thist
	setscale /p x,1,1,azhist, thist
	histogram /b=2 anglewave, thist
	histogram /b=2 /w=tempwave anglewave, azhist
	azhist /=thist
end
function azploter(en)
	variable en
	variable qi
	make/o /n=(360,.5/.003) azplot
	setscale /p x,1,1,azplot
	setscale /p y, .003,.003, azplot
	variable i
	for( i = 1 ; i <= dimsize(azplot,1) ; i+=1)
		azimuthplot(en,i*.003,i*.003+.003)
		wave azhist
		azplot[][i-1] = azhist[p]
	endfor
	azplot = azplot*0 != 0 ? 0 : azplot
	smooth /dim=0 10, azplot
	smooth /dim=1 10, azplot
	dowindow /k azplotwindow
	dowindow /k azgraphwindow
	azgraphwindowf()
//	azplotwindowf()
//	smooth /dim=0 10, azplot
end
function azgraphwindowf()
	wave azplot
	Display/k=1 /n=azgraphwindow /W=(637.2,74,811.2,224) azplot[*][20],azplot[*][30],azplot[*][40],azplot[*][50]
	AppendToGraph azplot[*][60],azplot[*][70],azplot[*][80],azplot[*][90],azplot[*][100]
	ModifyGraph margin(right)=10, margin(left)=35, margin(top)=10,gfSize=10
	ModifyGraph mode=7
	ModifyGraph rgb(azplot)=(0,0,0),rgb(azplot#1)=(16049,0,0),rgb(azplot#2)=(33436,0,0)
	ModifyGraph rgb(azplot#3)=(49485,0,0),rgb(azplot#5)=(65535,11762,0),rgb(azplot#6)=(65535,31927,0)
	ModifyGraph rgb(azplot#7)=(65535,53772,0),rgb(azplot#8)=(65535,65535,0)
	ModifyGraph hbFill=5
	ModifyGraph log(left)=1
	ModifyGraph tick=2
	ModifyGraph mirror=1
	ModifyGraph standoff=0
	ModifyGraph axOffset(left)=-2.22222,axOffset(bottom)=-0.666667
	ModifyGraph axisOnTop=1
	Label left "Scattering Intensity"
	Label bottom "Azimuthal Angle [deg]"
//	ColorScale/C/N=text0/A=RC/X=0.00/Y=0.00/E  ctab={0.6,3,YellowHot,0}, height=150
//	ColorScale/C/N=text0 width=8, tickLen=2, lblMargin=0, axisRange={3,0.6,0}
//	AppendText "Momentum Transfer [nm\\S-1\\M]"
End

function azplotwindowf()
	Display/k=1 /n=azplotwindow /W=(24,46.4,510.6,474.2)
	wave azplot
	AppendImage/T azplot
	ModifyImage azplot ctab= {*,*,SpectrumBlack,0}
	ModifyGraph margin(left)=28,margin(bottom)=14,margin(top)=26,margin(right)=14
	ModifyGraph mirror(left)=2,mirror(top)=0
	ModifyGraph nticks(left)=4,nticks(top)=8
	ModifyGraph minor=1
	ModifyGraph fSize=8
	ModifyGraph standoff=0
	ModifyGraph tkLblRot(left)=90
	ModifyGraph btLen=3
	ModifyGraph tlOffset=-2
	Label left "Momentum Transfer [nm\\S-1\\M]"
	Label top "Azimuthal Angle [deg]"
	SetAxis/A/R left
	ColorScale/C/N=text0/A=RC/X=3.21/Y=-11.22/E image=azplot
End
function nalignmap(m, density,slicen,rev)
	wave m, density
	variable slicen,rev
	string foldersave = getdatafolder(1)
	newdatafolder /o/s root:alignmap
	make/o /n=(dimsize(m,1)*dimsize(m,2)) xloc = 5*mod(p,dimsize(m,1)),yloc = 5*floor(p/dimsize(m,1))
	make/o /n=(dimsize(m,1)*dimsize(m,2)) xcomp,ycomp,zcomp
	make/o /n=(dimsize(m,1),dimsize(m,2)) slice = density[slicen][p][q][0]
	make/o /n=(dimsize(m,1)*dimsize(m,2),2) arrowsyay
	multithread xcomp = m[slicen][xloc/5][yloc/5][0]
	multithread ycomp = m[slicen][xloc/5][yloc/5][1]
	multithread zcomp = m[slicen][xloc/5][yloc/5][2]
	multithread arrowsyay[][1] = atan(zcomp[p]/ycomp[p])
	multithread arrowsyay[][0] = sqrt(zcomp[p]^2 + ycomp[p]^2) * 5
	//arrowsyay[][0] = arrowsyay[p][0]< 3 &&  sqrt(zcomp[p]^2 + ycomp[p]^2+ xcomp[p]^2)>.5? 3 : arrowsyay[p][0]
	yloc = arrowsyay[p][0] <4 ? nan : yloc[p]
	dowindow /k alignmentmap
	Display /k=1/n=alignmentmap/W=(270.6,59.6,518.4,308.6)/K=1  yloc vs xloc as "Slice of Alignment"
	//AppendImage/w=alignmentmap slice
	//if(rev)
	//	ModifyImage/w=alignmentmap slice ctab= {1,0,yellowhot,0}
	//else
	//	ModifyImage/w=alignmentmap slice ctab= {0,1,yellowhot,0}
	//endif
	ModifyGraph/w=alignmentmap mode=3,gfSize=12
	ModifyGraph/w=alignmentmap rgb=(0,0,0)
	ModifyGraph/w=alignmentmap msize=0.5,marker=42,tlOffset=-5
	ModifyGraph/w=alignmentmap mrkThick=0.1,margin(left)=26,margin(bottom)=26
	ModifyGraph/w=alignmentmap arrowMarker(yloc)={arrowsyay,.75,0,0,1}
	ModifyGraph/w=alignmentmap mirror=0,margin(top)=1,margin(right)=1
	ModifyGraph/w=alignmentmap axOffset(left)=-4.46154,axOffset(bottom)=-0.962963
	setscale /p x, 0,5,slice
	setscale /p y, 0,5,slice
	setaxis /w=alignmentmap left 0,min(1000,5*dimsize(slice,0))
	setaxis /w=alignmentmap bottom 0,min(1000,5*dimsize(slice,1))
	ModifyGraph /w=alignmentmap tick=2,mirror=1,axisOnTop=1,standoff=0,tkLblRot(left)=90
	Label /w=alignmentmap left "Size Scale [nm]";DelayUpdate
	Label /w=alignmentmap bottom "Size Scale [nm]"
	if(rev)
		make /n=3 /o tics = {1,.5,0}
	else
		make /n=3 /o tics = {0,.5,1}
	endif
	//make /n=3 /t /o ticlabels = {"P3HT", "Mixture", "PCBM"}
	//ColorScale /w=alignmentmap/C/N=text0/A=RB image=slice;DelayUpdate
	//ColorScale /w=alignmentmap/C/N=text0 userTicks={tics,ticlabels}
	//ColorScale /w=alignmentmap/C/N=text0 vert=0,side=2,tickLblRot=0
	//ColorScale /w=alignmentmap/C/N=text0/X=3.00/Y=3.00 width=100,height=15
	//ColorScale /w=alignmentmap/C/N=text0 lblMargin=5,tickLen=3.00
	setdatafolder foldersave
end


Window alignmentmap() : Graph
	PauseUpdate; Silent 1		// building window...
	String fldrSav0= GetDataFolder(1)
	SetDataFolder root:alignmap:
	Display /W=(-1861,-486,-457,918)/K=1  yloc vs xloc as "Slice of Alignment"
	AppendImage ::packages:ScatterSim3D:lamglong2:working:thickness
	ModifyImage thickness ctab= {*,*,Grays,0}
	ModifyImage thickness plane= 3
	SetDataFolder fldrSav0
	ModifyGraph margin(left)=26,margin(bottom)=26,margin(top)=1,margin(right)=1,gfSize=12
	ModifyGraph height={Plan,1,left,bottom}
	ModifyGraph mode=3
	ModifyGraph marker=42
	ModifyGraph rgb=(3,52428,1)
	ModifyGraph msize=0.5
	ModifyGraph mrkThick=0.1
	ModifyGraph arrowMarker(yloc)={:alignmap:arrowsyay,0.75,0,0,1}
	ModifyGraph tick=2
	ModifyGraph mirror=1
	ModifyGraph standoff=0
	ModifyGraph axOffset(left)=-4.46154,axOffset(bottom)=-0.962963
	ModifyGraph tkLblRot(left)=90
	ModifyGraph axisOnTop=1
	ModifyGraph tlOffset=-5
	Label left "Size Scale [px]"
	Label bottom "Size Scale [px]"
EndMacro
function doalignmentmap(folder,[addtolayout])
	string folder
	variable addtolayout
	
	addtolayout = paramisdefault(addtolayout) ? 0 : addtolayout
	if(!stringmatch("Current",folder))
		string foldersave = getdatafolder(1)
		setdatafolder root:Packages:ScatterSim3D:$folder
	endif
	wave m=m3save
	wave /z m2 = m1save
	wave density1
	wave /z density2
	wave /z density3
	wave /z density4
	make/o /n=(dimsize(m,1)*dimsize(m,2)) xloc = 5*mod(p,dimsize(m,1)),yloc = 5*floor(p/dimsize(m,1))
	duplicate /o yloc, yloc2
	make/o /n=(dimsize(m,1)*dimsize(m,2)) xcomp,ycomp,zcomp
	make/o /n=(dimsize(m,1)*dimsize(m,2),2) arrowsyay, arrowsyay2
	MAKE/O /N=(dimsize(m,1),dimsize(m,2)) slice = density1[5][p][q]
	if(waveexists(density2))
		slice += density2[5][p][q]
	endif
	if(waveexists(density3))
		slice += density1[5][p][q]
	endif
	setscale /i x, 0,dimsize(m,1)*5,slice
	setscale /i y, 0,dimsize(m,1)*5,slice
	
	multithread xcomp = m[5][xloc/5][yloc/5][0]
	multithread ycomp = m[5][xloc/5][yloc/5][1]
	multithread zcomp = m[5][xloc/5][yloc/5][2]
	multithread arrowsyay[][1] = atan(zcomp[p]/ycomp[p])
	multithread arrowsyay[][0] = sqrt(zcomp[p]^2 + ycomp[p]^2) * 20
	if(waveexists(m2))
		multithread xcomp = m2[5][xloc/5][yloc/5][0]
		multithread ycomp = m2[5][xloc/5][yloc/5][1]
		multithread zcomp = m2[5][xloc/5][yloc/5][2]
		multithread arrowsyay2[][1] = atan(zcomp[p]/ycomp[p])
		multithread arrowsyay2[][0] = sqrt(zcomp[p]^2 + ycomp[p]^2) * 20
		yloc2 = arrowsyay2[p][0] <4 ? nan : yloc2[p]
	else
		arrowsyay2=nan
	endif
	//arrowsyay[][0] = arrowsyay[p][0]< 3 &&  sqrt(zcomp[p]^2 + ycomp[p]^2+ xcomp[p]^2)>.5? 3 : arrowsyay[p][0]
	yloc = arrowsyay[p][0] <4 ? nan : yloc[p]
	
	string alignmapname = cleanupname((folder + "_map"),0)
	dowindow /k $alignmapname
	Display /k=1/n=$alignmapname/W=(-1861,-486,-457,918) /l /b  yloc vs xloc as "Slice of Alignment"
	appendtograph /w=$alignmapname /l /b  yloc2 vs xloc
	AppendImage/w=$alignmapname /l/b slice
	//if(rev)
		ModifyImage/w=$alignmapname slice ctab= {1,0,grays,0}
	//else
	//	ModifyImage/w=alignmentmap slice ctab= {0,1,yellowhot,0}
	//endif
	ModifyGraph/w=$alignmapname mode=3,gfSize=12
	ModifyGraph/w=$alignmapname rgb(yloc2)=(65535,32768,32768),rgb(yloc)=(16385,65535,65535)
	ModifyGraph/w=$alignmapname msize=0.5,marker=42,tlOffset=-5
	ModifyGraph/w=$alignmapname mrkThick=0.1,margin(left)=26,margin(bottom)=26
	ModifyGraph/w=$alignmapname arrowMarker(yloc)={arrowsyay,.75,0,0,1}
	ModifyGraph/w=$alignmapname arrowMarker(yloc2)={arrowsyay2,.75,0,0,1}
	ModifyGraph/w=$alignmapname mirror=0,margin(top)=1,margin(right)=1
	ModifyGraph/w=$alignmapname axOffset(left)=-4.46154,axOffset(bottom)=-0.962963
	setaxis /w=$alignmapname left 0,min(1000,5*dimsize(m,1))
	setaxis /w=$alignmapname bottom 0,min(1000,5*dimsize(m,2))
	ModifyGraph /w=$alignmapname mirror(bottom)=0,nticks=0,noLabel=2,axisOnTop=0
	ModifyGraph /w=$alignmapname margin=1,height={Plan,1,left,bottom}
	
	if(addtolayout)
		appendlayoutobject /F=0 graph $alignmapname
		dowindow /HIDE=1 $alignmapname
	endif
	if(!stringmatch("Current",folder))
		setdatafolder foldersave
	endif
end
function scatterimage(folder,en,[addtolayout,qpwr,append2graph])
	string folder
	variable en
	variable addtolayout
	variable qpwr, append2graph
	qpwr = paramisdefault(qpwr)? 2:qpwr
	append2graph = paramisdefault(append2graph)? 0:append2graph
	addtolayout = paramisdefault(addtolayout) ? 0 : addtolayout
	string foldersave = getdatafolder(1)
	setdatafolder root:Packages:ScatterSim3D:$folder
	wave scatter3dsave
	wave Int3DvsEn
	wave perp3Dvsen
	wave para3dvsen
	wave enwave
	if(!waveexists(scatter3dsave))
		print "no scattering in this directory"
		return 0
	endif
	if(!waveexists(enwave))
		wave enwavedisp
		if(!waveexists(enwavedisp))
			print "no energy reference in directory"
			return 0
		endif
		duplicate enwavedisp, enwave
		deletepoints 0,1,enwave
	endif
	string scatname = cleanupname("Scatter"+num2str(en),0)
	string scat1Dname = cleanupname("Scatter1D"+num2str(en),0)
	string Para1Dname = cleanupname("Para1D"+num2str(en),0)
	string Perp1Dname = cleanupname("Perp1D"+num2str(en),0)
	
	make /n=(1000,1000) /o $scatname
	//make /n=1000 /o $scat1Dname, $Para1Dname, $Perp1Dname
	wave scatterdisp = $scatname
	//setscale /i x,.015,1, scat1D, para1D, perp1d
	setscale /i x,-1,1, scatterdisp
	setscale /i y,-1,1, scatterdisp
	//setscale /i z,260,320,scatter3dsave
	//setscale /i x,-.5,.5,scatter3dsave
	//setscale /i y,-.5,.5,scatter3dsave
	variable enstep = binarysearch(enwave,en)
	scatterdisp = sqrt(x^2 +y^2) < .005 ? 0 : scatter3dsave(x)(y)[enstep]
	//scat1D = Int3DvsEn(x)
	//perp1d = perp3DvsEn(x)
	//para1D = para3Dvsen(x)
	matrixfilter /n=7 gauss scatterdisp
	
	radialintegratew(scatterdisp,0,90,scat1dname)
	radialintegratew(scatterdisp,0,20,para1Dname)
	radialintegratew(scatterdisp,70,90,perp1Dname)
	
	
	wave scat1D = $scat1dname
	wave para1D = $para1Dname
	wave perp1d = $perp1Dname
	scat1D *=x^qpwr
	para1D *=x^qpwr
	perp1D *=x^qpwr
	string svname = cleanupname(folder+num2str(en),0)
	string sgname = cleanupname(folder+"1D"+num2str(en),0)
	if(append2graph)
		appendtograph para1D /tn=$(folder + para1Dname)
		appendtograph perp1D /tn=$(folder + perp1Dname)
		appendtograph scat1d /tn=$(folder + scat1dname)
		
		ModifyGraph mode($(folder + para1Dname))=7,useNegPat($(folder + para1Dname))=1,toMode($(folder + para1Dname))=1,hbFill($(folder + para1Dname))=2,hBarNegFill($(folder + para1Dname))=2,rgb($(folder + para1Dname))=(0,0,65535)
		ModifyGraph rgb($(folder + scat1dname))=(0,0,0),lsize($(folder + para1Dname))=0,lsize($(folder + perp1Dname))=0,lsize($(folder + scat1dname))=2,log=1
		ModifyGraph grid=1,tick=2,gfSize=20,axThick=2,standoff=0,gridStyle=3,gridRGB=(43690,43690,43690),useNegRGB($(folder + para1Dname))=1,usePlusRGB($(folder + para1Dname))=1,plusRGB=(65535,0,0),negRGB($(folder + para1Dname))=(3,1,52428)
	
	else
		dowindow /k $sgname
		display /W=(887,57,1623,553)/n=$sgname /k=1 para1D, perp1D, scat1d as (folder + " " + num2str(en) + " eV 1D Scattering")
		SetAxis /w=$sgname /A=2 left //0.00001,50
		SetAxis /w=$sgname bottom 0.02,1
		ModifyGraph /w=$sgname mode($para1Dname)=7,useNegPat($para1Dname)=1,toMode($para1Dname)=1,hbFill($para1Dname)=2,hBarNegFill($para1Dname)=2,rgb($para1Dname)=(0,0,65535)
		ModifyGraph /w=$sgname rgb($scat1dname)=(0,0,0),lsize($para1Dname)=0,lsize($perp1Dname)=0,lsize($scat1dname)=2,log=1
	
	//Label /w=$sgname left "Scattering Intensity";DelayUpdate
	//Label /w=$sgname bottom "Momentum Transfer Q [nm\\S-1\\M]"
		ModifyGraph /w=$sgname grid=1,tick=2,gfSize=20,axThick=2,standoff=0,gridStyle=3,gridRGB=(43690,43690,43690),useNegRGB($para1Dname)=1,usePlusRGB($para1Dname)=1,plusRGB=(65535,0,0),negRGB($para1Dname)=(3,1,52428)
	//Legend /w=$sgname/C/N=text0/J/F=0/B=1 "\\s("+scat1dname+") Azimuthal Integration\r\\s("+para1Dname+") Positive Anisotropy"
	endif
	dowindow /k $svname
	display /k=1 /n=$svname /W=(40,45,486,479) as (folder + " " + num2str(en) + " eV 2D Scattering")
	appendimage /w=$svname scatterdisp
	SetAxis /w=$svname left -0.4,0.4
	SetAxis /w=$svname bottom -0.4,0.4
	ModifyImage /w=$svname $scatname ctab= {.0001,5,Terrain,0},log=1
	ModifyGraph /w=$svname height={Plan,1,left,bottom}
	//ColorScale/w=$svname /C/N=text0/A=RC/X=0.00/Y=0.00 log=1,image=$scatname
	//ColorScale/w=$svname /C/N=text0 "Scattering Intensity"
	TextBox/w=$svname/C/N=text1/F=0/A=LT/X=10.00/Y=10.00/B=1/G=(0,0,0) "\\F'Arial Black'\\Z20"+num2str(en)+" eV"
	setdatafolder foldersave
	
	if(addtolayout && !append2graph)
		appendlayoutobject /F=0 graph $svname
		appendlayoutobject /F=0 graph $sgname
		dowindow /HIDE=1 $svname
		dowindow /HIDE=1 $sgname
	endif
	//azploter(en)
end


function doratiograph(folder,[addtolayout])
	string folder
	variable addtolayout
	addtolayout = paramisdefault(addtolayout) ? 0 : addtolayout
	string foldersave = getdatafolder(1)
	setdatafolder root:Packages:ScatterSim3D:$folder
	duplicate/o ratio3dvsen, ratio3dtemp
	//make/o /n=(dimsize(ratio3dtemp,0)) /o normalizewave
	//normalizewave = ratio3dtemp[p][0]
	//ratio3dtemp -= normalizewave[p]
	wave enwavedisp
	dowindow /k $("ratiograph_"+folder)
	Display /W=(103,49,476,547) /n=$("ratiograph_"+folder)
	AppendImage/T ratio3dtemp vs {*,enwavedisp}
	ModifyImage ratio3dtemp ctab= {-1,1,RedWhiteBlue,0}
	ModifyGraph margin(left)=43,margin(bottom)=9,margin(top)=28,margin(right)=9,gfSize=14
	ModifyGraph nticks(top)=8
	ModifyGraph standoff=0
	ModifyGraph tlOffset=-2
	Label left "X-ray Energy [eV]"
	SetAxis/R left 275,295
	SetAxis top 0.01,1
	smooth 5, ratio3dtemp
	ModifyGraph log(top)=1
	ModifyGraph mirror=2
	ModifyGraph minor=1
	ModifyGraph fSize=8
	ModifyGraph standoff=0
	ModifyGraph btLen=3
	ModifyGraph tlOffset=-2
	ModifyImage ratio3dtemp ctab= {-1,1,RedWhiteGreen,0}
	ModifyGraph fSize=0
	ModifyGraph gfSize=14
	ModifyGraph tkLblRot=0
	ModifyGraph nticks(left)=6
	ModifyGraph tick=2,btLen=5
	ModifyGraph ftLen(left)=5
	ModifyGraph sep(left)=2
	dowindow /k $("Ingraph_"+folder)
	Display /W=(103,49,476,547) /n=$("Ingraph_"+folder)
	appendimage int3DvsEn vs {*,enwavedisp}
	ModifyImage int3DvsEn ctab= {1e-08,0.1,YellowHot,0},log=1
	ModifyGraph log(bottom)=1
	SetAxis bottom 0.015,1
	SetAxis/R left 275,295
	setdatafolder foldersave
	if(addtolayout)
		appendlayoutobject /F=0 graph $("ratiograph_"+folder)
		appendlayoutobject /F=0 graph $("Ingraph_"+folder)
		dowindow /HIDE=1 $("Ingraph_"+folder)
		dowindow /HIDE=1 $("ratiograph_"+folder)
	endif
	
end

function sysrotate(mext,mout,angle, offsetx, offsety)
	wave mext
	wave mout  // the output wave, which should be a duplicate of the original m wave
	variable angle
	variable offsetx
	variable offsety
	if(angle==0)
		mout = mext[p][q+offsetx][r+offsety][s]
		return 0
	endif
	imagetransform /TM4D=2418 transpose4D mext
	wave M_4DTranspose
	make /o /n=(dimsize(M_4DTranspose,0),dimsize(M_4DTranspose,1),dimsize(M_4DTranspose,2)) w0=M_4DTranspose[p][q][r][0],w1=M_4DTranspose[p][q][r][1],w2=M_4DTranspose[p][q][r][2],w3=M_4DTranspose[p][q][r][3] , wtot
//	imageinterpolate /func=spline /TRNS={scaleshift,0,2,0,2} resample w0
//	duplicate /o M_InterpolatedImage, w0
	//wavestats/q w0
	multithread wtot = w0[p][q][r]^2 + w1[p][q][r]^2 + w2[p][q][r]^2 + w3[p][q][r]
	imagerotate /o/E=(0)/A=(angle) wtot
	imagerotate /o/E=(0)/A=(angle) w0
	//w0 = w0[p][q][r]<v_min? v_min : w0[p][q][r]
	//w0 = w0[p][q][r]>v_max? v_max : w0[p][q][r]
//	imageinterpolate /PXSZ={2,2} pixelate w0
//	duplicate /o M_PixelatedImage, w0
//	imageinterpolate /func=spline /TRNS={scaleshift,0,2,0,2} resample w1
//	duplicate /o M_InterpolatedImage, w1
	//wavestats/q w1
	imagerotate /o/E=(nan)/A=(angle) w1
	//w1 = w1[p][q][r]<v_min? v_min : w1[p][q][r]
	//w1 = w1[p][q][r]>v_max? v_max : w1[p][q][r]
//	imageinterpolate /PXSZ={2,2} pixelate w1
//	duplicate /o M_PixelatedImage, w1
//	imageinterpolate /func=spline /TRNS={scaleshift,0,2,0,2} resample w2
//	duplicate /o M_InterpolatedImage, w2
	//wavestats/q w2
	imagerotate /o/E=(nan)/A=(angle) w2
	//w2 = w2[p][q][r]<v_min? v_min : w2[p][q][r]
	//w2 = w2[p][q][r]>v_max? v_max : w2[p][q][r]
//	imageinterpolate /PXSZ={2,2} pixelate w2
//	duplicate /o M_PixelatedImage, w2
//	imageinterpolate /func=spline /TRNS={scaleshift,0,2,0,2} resample w3
//	duplicate /o M_InterpolatedImage, w3
	//wavestats/q w3
	imagerotate /o/E=(nan)/A=(angle) w3
	//w3 = w3[p][q][r]<v_min? v_min : w3[p][q][r]
	//w3 = w3[p][q][r]>v_max? v_max : w3[p][q][r]
//	imageinterpolate /PXSZ={2,2} pixelate w3
//	duplicate /o M_PixelatedImage, w3
//	imagefilter /n=5/o hybridmedian w0
//	imagefilter /n=5/o hybridmedian w1
//	imagefilter /n=5/o hybridmedian w2
//	imagefilter /n=5/o hybridmedian w3
	duplicate/o w0, wtotnew
	multithread wtotnew = w0[p][q][r]^2 + w1[p][q][r]^2 + w2[p][q][r]^2 + w3[p][q][r]
	multithread wtotnew = (wtotnew[p][q][r]<wtot[p][q][r] || wtotnew[p][q][r] >wtotnew[p][q][r]) && wtot[p][q][r] >0 ? wtotnew[p][q][r]/wtot[p][q][r] : 1
	multithread w0 /= sqrt(abs(wtotnew[p][q][r]))
	multithread w1 /= sqrt(abs(wtotnew[p][q][r]))
	multithread w2 /= sqrt(abs(wtotnew[p][q][r]))
	multithread w3 /= abs(wtotnew[p][q][r])
	
	offsetx += floor((dimsize(w0,0)-dimsize(mout,1))/2)
	offsety += floor((dimsize(w0,1)-dimsize(mout,2))/2)
	
	//multithread mout[][][][0] = w0[q+floor((dimsize(w0,0)-dimsize(mout,1))/2)][r+floor((dimsize(w0,1)-dimsize(mout,2))/2)][p]
	//make /n=(dimsize(mout,0),dimsize(mout,1),dimsize(mout,2)) /o tempy
	//multithread tempy = w1[q+floor((dimsize(w0,0)-dimsize(mout,1))/2)][r+floor((dimsize(w0,1)-dimsize(mout,2))/2)][p]
	//multithread mout[][][][2] = w2[q+floor((dimsize(w0,0)-dimsize(mout,1))/2)][r+floor((dimsize(w0,1)-dimsize(mout,2))/2)][p]
	//multithread mout[][][][3] = w3[q+floor((dimsize(w0,0)-dimsize(mout,1))/2)][r+floor((dimsize(w0,1)-dimsize(mout,2))/2)][p]
	
	multithread mout[][][][0] = w0[q+offsetx][r+offsety][p]
	make /n=(dimsize(mout,0),dimsize(mout,1),dimsize(mout,2)) /o tempy
	multithread tempy = w1[q+offsetx][r+offsety][p]
	multithread mout[][][][2] = w2[q+offsetx][r+offsety][p]
	multithread mout[][][][3] = w3[q+offsetx][r+offsety][p]
	
	
	
	
	
	multithread mout[][][][1] = cos(angle*pi/180)*tempy[p][q][r] - sin(angle*pi/180) *mout[p][q][r][2]// y' = cos(th) y - sin(th) z
	multithread mout[][][][2] = cos(angle*pi/180)*mout[p][q][r][2] + sin(angle*pi/180) *tempy[p][q][r]// z' = cos(th) z + sin(th) y

end


function /s variables_coreshell()
	string variables = "Total Radius is Length scale above,String,^;"
	variables +=  "Core Interpenetration [pixels],SetVariable,0;"
	variables += "Shell Interpenetration [pixels],SetVariable,1;"
	variables += "Minimum Seperation,SetVariable,.1;"
	variables += "Number of Particles (Max),SetVariable,500;"
	variables += "sigma of radiuses,setvariable,.1;"
	variables += "Noise,SetVariable,0;"
	variables += "Shell Width [px],SetVariable,1;"
	variables += "Maximum radius [px],SetVariable,15;"
	variables += "Minimum radius [px],SetVariable,2;"
	return variables
end
function model3D_coreshell(s3d)
	//Creates a spherical system, with two components, aligned 
	// utilizes extra parameters, number of particles, polydispursity (a sigma of radiuses) and populates the 3D space with 
		//aligned particles
	// material 1 is insize the spheres, material 2 is outside (although the concentrations are determined by m2volfrac
	//paramstring = 
		//MORPHOLOGY parameters
		//0 interpenetration of spheres (roughness smaller than wavelength)
		//1 minimum seperation (when placing a new sphere into the system, it needs to be this far away from any others)
			//	this can be a fraction of the radius (if less than 1) or a flat number of pixels (if greater than 1)
		//2 number of particles // the maximum number of particles . the program will stop when it cannot fit anymore
		//3 polydispursity (sigma of radiuses)
		//4 thickness (this can be different than size)
		//5 noise (percentage of vacuum - density variations)
	//num = resolution (in each of 3 dimensions)
	//size = average size of sphere (radius)
	
	struct ThreeDSystem &s3d
	if(itemsinlist(s3d.paramstring,",")<5)
		return -1
	endif
	newdatafolder /o/s CreatingSpheres
	variable coreinterpenetration = 	str2num(stringfromlist( 1 ,s3d.paramstring,","))
	variable shellinterpenetration = 	str2num(stringfromlist( 2 ,s3d.paramstring,","))
	variable minsep = 			str2num(stringfromlist( 3 ,s3d.paramstring,","))
	variable particlesnum = 	str2num(stringfromlist( 4 ,s3d.paramstring,","))
	variable pd = 				str2num(stringfromlist( 5 ,s3d.paramstring,","))
	variable thickness = 		s3d.thickness
	variable noise = 			str2num(stringfromlist( 6 ,s3d.paramstring,","))
	variable shellwidth = 		str2num(stringfromlist( 7 ,s3d.paramstring,","))
	variable maxradius = 		min(str2num(stringfromlist( 8 ,s3d.paramstring,",")),thickness/2)
	variable minradius = 		str2num(stringfromlist( 9 ,s3d.paramstring,","))
	

	
	make /o /n=(thickness,s3d.num,s3d.num) mat=1, core=0,xwave, ywave, zwave
	
	make/B/U /o /n=(thickness,s3d.num,s3d.num,maxradius) exmat= (p <= t) || (q <= t) || (r <= t) || (p >= thickness-t) || (q >= s3d.num-t) || (r >= s3d.num-t) ? 0 : 1
	make/B/U /o /n=(thickness,s3d.num,s3d.num) tempwave
	if(s3d.movie)
		Execute("Spheres3Ddisp(" +num2str(s3d.num)+", \""+getwavesdatafolder(mat,2)+"\")")
		Execute("exportgizmo wave=\"testimage\"  ;Spinoidal3DLayout();Spinoidal3DImage(\""+getdatafolder(1)+"testimage\")")
	endif
	setscale /i x, -thickness/2, thickness/2, mat,core, exmat,xwave, ywave, zwave
	setscale /i y, -s3d.num/2, s3d.num/2, mat,core, exmat,xwave, ywave, zwave
	setscale /i z, -s3d.num/2, s3d.num/2, mat,core, exmat,xwave, ywave, zwave
	xwave = x
	ywave = y
	zwave = z
	redimension /n=(thickness*s3d.num*s3d.num) xwave, ywave, zwave
	variable i,radius, orad, cx, cy, cz, failed, fnum =0, xmn,xmx,ymn,ymx,zmn,zmx, loc
	
	for(i=0;i<particlesnum;i+=1)
		fnum=0
		do
			failed = 0
			radius = abs(gnoise(pd)+s3d.size)
			radius =radius < shellwidth ? shellwidth : radius
			radius = min(maxradius,max(max(shellwidth,minradius),radius))
			if(minsep<1)
				orad = radius*(1+minsep/2)
			else
				orad = radius + minsep/2
			endif
			//duplicate/o /r=()()()(ceil(2*orad)) exmat,tempwave
			redimension /n=(thickness,s3d.num,s3d.num) tempwave
			multithread tempwave[][][] = exmat[p][q][r][ceil(orad)]
//			imagefilter /n=(ceil(2*orad)) /o min3d tempwave
//			multithread tempwave = tempwave<1? 0 : 1 
//			multithread tempwave = (p <= orad) || (q <= orad) || (r <= orad) || (p >= thickness-orad) || (q >= s3d.num-orad) || (r >= s3d.num-orad) ? 0 : tempwave[p][q][r]
  			if(wavemax(tempwave)<1)
				//there are no possible locations for this radius, find another
				failed=1
			else
				// randomly pick a pixel that is good for the center
				redimension /n=(thickness*s3d.num*s3d.num) tempwave
				integrate tempwave /D=intwave
				loc = binarysearch(intwave, enoise(wavemax(intwave)/2)+wavemax(intwave)/2)
				cx = xwave[loc]
				cy = ywave[loc]
				cz = zwave[loc]
			endif
			
			if(failed && (fnum>10))
				print "warning : failed to find unoccupied location - only " +num2str(i) +" particles created"		
				break // get out of the loop
				//we have failed 30 times to create a radius or find a location
			endif
			fnum +=1
		while(failed==1)
		if(failed)
			break // do not add this sphere to the system
					// we are done adding spheres
		endif
		// subtract out this sphere from the matrix  // matrix starts at 1s, within this sphere, multiply this by 0, outside multiply by 1
		multithread mat*= (x-cx)^2 + (y-cy)^2 + (z-cz)^2 < radius^2 ? 0 : 1 
		//multithread core*= (x-cx)^2 + (y-cy)^2 + (z-cz)^2 < radius^2 && (x-cx)^2 + (y-cy)^2 + (z-cz)^2 > (radius-shellwidth)^2 ? 1 : 0 
		multithread core+= (x-cx)^2 + (y-cy)^2 + (z-cz)^2 < (radius-shellwidth)^2 ? 1 : 0
		multithread exmat*= (x-cx)^2 + (y-cy)^2 + (z-cz)^2 <= (orad+t)^2 ? 0 : 1 
		
		if(s3d.movie)
			execute("ModifyGizmo /n=Spheres3D update=2")
			doupdate
			Execute "exportgizmo wave=\"testimage\"   "
			//TextBox/w=Spinoidal3DLayout/C/N=text0/A=LT/X=0.00/Y=0.00 "\Z32" + time2str2(ttot)
			doupdate
			savepict /p=_PictGallery_ /E=-5 /N=Spinoidal3DLayout /o as "Frame3D"
			addmovieframe /pict=Frame3D
		endif
	endfor
	setdatafolder ::
	imagefilter /n=(shellinterpenetration)/o gauss3d mat
	imagefilter /n=(coreinterpenetration)/o gauss3d core
	duplicate /o mat,s3d.density1 // this returns the density matrix of material 1 (the matrix) for alignment etc later on
	duplicate /o core,s3d.density2 // this returns the density matrix of material 1 (the matrix) for alignment etc later on
end