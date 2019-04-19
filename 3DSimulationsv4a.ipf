#pragma rtGlobals=1		// Use modern global access method.
#include "opticalconstantsDB"

structure ThreeDSystem
	wave m1 //4 dimensional density of material 1 in (x direction, y direction, z direction, unaligned)  sum = relative concentration sum(between -1 - 1 for x, y and z)^2 + unaligned = total material density
	wave m2 //4 dimensional density of material 2 in (x direction, y direction, z direction, unaligned)  sum = relative concentration sum(between -1 - 1 for x, y and z)^2 + unaligned = total material density
	wave m3 //4 dimensional density of material 3 in (x direction, y direction, z direction, unaligned)  sum = relative concentration sum(between -1 - 1 for x, y and z)^2 + unaligned = total material density
	wave m4 //4 dimensional density of material 4 in (x direction, y direction, z direction, unaligned)  sum = relative concentration sum(between -1 - 1 for x, y and z)^2 + unaligned = total material density
	wave m5 //4 dimensional density of material 5 in (x direction, y direction, z direction, unaligned)  sum = relative concentration sum(between -1 - 1 for x, y and z)^2 + unaligned = total material density
	
	wave /d/c n //3D complex index of refraction map (depends on energy) // this is now outdated
				// in the current version of the program, this will be a single complex number for each material which will be set at each energy step, and used in calculations
	
	wave /d /c px // induced polarization in the x direction for real space, x,y,z
	wave /d /c pxFFT // fourier transform of above, in qx, qy, qz space
	wave /d /c py // same in y direction
	wave /d /c pyFFT// fourier transform of above, in qx, qy, qz space
	wave /d /c pz // same in z direction
	wave /d /c pzFFT// fourier transform of above, in qx, qy, qz space
	//wave /d qtensor // the coeffifients of the PFFT which are used to calculate EScat components below // this is calculated one for the simulation
	// the x,y,z components of qtensor correspond to it's native units of qx, qy, qz
		//qtensor[][][][0] = - k * (2 * x + y + z) - x * (x + y + z)
		//qtensor[][][][1] =   k^2 - k * y - y * (x + y + z)
		//qtensor[][][][2] =   k^2 - k * z - z * (x + y + z)
		// This is matrix multiplied by ((Pxfft,Pyfft,Pzfft (which are each 3D scalar fields with coordinates in qx, qy, and qz)) producing a 3 component vector field in qx, qy, qz components.
	wave /d /c Escatx, Escaty, Escatz // the result of the above line
		// Escatx = pxFFT[x][y][z] * qtensor[x][y][z][0]
		// Escatx = pxFFT[x][y][z] * qtensor[x][y][z][0]
		// Escatx = pxFFT[x][y][z] * qtensor[x][y][z][0]
		// the elastic scattering condition is applied to these components (interpolate and magsqr 2D qy, qz components, where qx = sqrt(k^2-qz^2-qy^2)-k)
	wave /d EscatSqr // the sum of the magsquared values as described above, before interpolation
//	wave/d nr  // real component of n // depricated in this version, not useful
//	wave/d ni  // imaginary component of n // depricated in this version, not useful
//	wave /d/c np // 2D complex index of refraction projection (or slice) (along x axis) // depricated in this version, not useful
	
	wave density1 // the density of material 1 output from the morphology module, this will be changed into
	wave density2 // optional density of material 2 (implies 3 materials)
	wave density3 // implies 4 materials
	wave density4 // implies 5 materials
	
	wave/d scatter3D // 2d scattering simulation from 3D data
	//wave/d scatter2D // 2d scattering simulation from projected data // no longer calculating this in the new version
	
	wave/d int3D // radially summed intensity
	// wave/d int2D // radially summed intensity  // no longer calculating this in the new version
	wave/d int3Dp0 // radially averaged intensity
	// wave/d int2Dp0 // radially averaged intensity  // no longer calculating this in the new version
	wave/d int3Dpara // radially summed intensity in direction of electric field
	// wave/d int2Dpara // radially summed intensity in direction of electric field  // no longer calculating this in the new version
	wave/d int3Dp0para // radially averaged intensity in direction of electric field
	// wave/d int2Dp0para // radially averaged intensity in direction of electric field  // no longer calculating this in the new version
	wave/d int3Dperp // radially summed intensity in direction normal to the electric field
	// wave/d int2Dperp // radially summed intensity in direction normal to the electric field  // no longer calculating this in the new version
	wave/d int3Dp0perp // radially averaged intensity in direction normal to the electric field
	// wave/d int2Dp0perp // radially averaged intensity in direction normal to the electric field  // no longer calculating this in the new version
	
	
	wave/d Scatter3DSave // Full 2D Scattering Scattering Pattern Saved vs Energy
	// wave /d int2DvsEn // two dimensional average intensity vs momentum transfer and energy
	wave /d int3DvsEn // same but for the 3D calculated data
	// wave /d ratio2DvsEn // the ratio of parallel to perpindicular scattering intensity vs Q and Energy
	wave /d ratio3DvsEn // same but for the 3D calculated data
	wave /d Para3DvsEn // same but for the 3D calculated data
	wave /d Perp3DvsEn // same but for the 3D calculated data
	
	variable step // current energy step (integer for storing data)
	variable en // current energy = 1239.8 / wavelength =  ( 1239.8/(2*pi)) * wavevector
	variable wavelength // current wavelength = 1239.8 / energy = 2*pi / wavevector (in nanometers)
	variable k // current wavevector magnitude (2*pi / wavelength) = (2 * pi / 1239.8 ) * Energy (in inverse nanometers)
	
	variable num // resolution of simulation (in each of three dimensions)
	variable thickness
	variable voxelsize // the size of each pixel in the simulation (in nanometers)
	variable size // the dominant size of features (ie the radius of the sphere, if we're simulating spheres)
	variable materialnum
	string paramstring //specific parameters for creating density matricies according to the model
	string modelname
	string materials // list of materials
	string efield // two components which add to 1 which indicate the unit direction of the efield in the (ydirection, z direction)
				// x-ray propogation is always in the x direction
	string materialalignmentstring // list of list 
	string path // the path to put important things for this calculation
	variable movie
	variable timer
	variable timesofar
endstructure

function model3D(modelname,voxelsize,sizescale,resolution,thickness,paramstring,materiallist,materialalignmentstring,efield,energymin,energymax,energynum[movie, save3d])
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
	string modelname,paramstring,materiallist,materialalignmentstring,efield
	variable sizescale, resolution, energymin,energymax,energynum, movie, save3d, thickness,voxelsize
	struct ThreeDSystem s
	variable result, timer1
	for(timer1=0;timer1<11;timer1+=1)
		s.timesofar = stopmstimer(timer1)
	endfor
	s.voxelsize = voxelsize
	s.timesofar = 0
	s.timer = startmstimer
	s.thickness = thickness
	s.num = resolution
	s.size = sizescale
	s.materials = materiallist
	s.efield = efield
	s.modelname = modelname
	s.paramstring = paramstring
	s.movie = movie
	s.materialalignmentstring = materialalignmentstring
	s.materialnum = itemsinlist(materiallist,",")
	dowindow /k simulation_layout
	execute "simulation_layout()"
	if(s.movie)
		Variable DebugEnab
		DebuggerOptions
		DebugEnab = V_debugOnError //check for debug on error
		if (DebugEnab) //if it is on,
			DebuggerOptions debugOnError=0 //turn it off
		endif
		try
			closemovie;AbortOnRTE
		catch
			Print "Creating new movie"
		endtry
		if(DebugEnab)
			DebuggerOptions debugOnError=1
		endif
		SavePict/O/WIN=Simulation_Layout /E=-5/P=_PictGallery_/w=(0,0,100*11,100*8) as "SimPict"
		Newmovie /PICT=SimPict
	endif
	s.path = getdatafolder(1)
//	if(stringmatch("Spheres2",modelname))
//		print "Creating system of aligned spheres"
//		model3D_Spheres2(s)
//	elseif(stringmatch("existing",modelname))
//		print "Using Existing Density Waves"
//		model3D_Existing(s)
//	elseif(stringmatch("spinoidal",modelname))
//		print "Creating a spinoidal depomposition system"
//		model3D_spinoidal(s)
//	else
//		print "no recognized model could be created"
//		return -1
//	endif 
	if(exists("model3d_"+modelname)==6)
		funcref model3d_existing creatfunc=$("model3d_"+modelname)
		creatfunc(s)
	else
		print "no recognized model could be created"
		return -1
	endif
	variable alignmentincluded
	if(exists("special_"+modelname)==6)
		funcref special_existing specfunc=$("special_"+modelname)
		alignmentincluded = stringmatch(specfunc(),"*IncludesAlignment*")
	endif
	
	s.timesofar =stopmstimer(s.timer)/10^6
	s.timer = startmstimer
	//calculate alignment of system
	if( stringmatch(materialalignmentstring,"none") )
		Print "Time : "+time2str(s.timesofar) +"  -   Loading Existing material Alignment"
		loadexistingmaterialalignment(s)
	elseif(  !alignmentincluded   )
		Print "Time : "+time2str(s.timesofar) +"  -   Calculating alignment in system"
		Align3Dsystem(s)
	endif
	s.timesofar +=stopmstimer(s.timer)/10^6
	s.timer = startmstimer
	Print "Time : "+time2str(s.timesofar) +"  -   Alignment Loaded - Checking System Density for consistancy"
	if(sum3dsystem(s)<0)
		s.timesofar +=stopmstimer(s.timer)/10^6
		s.timer = startmstimer
		Print "Time : "+time2str(s.timesofar) +"  -   warning : System Creation failed to conserve density throughout system"
	elseif(sum3dsystem(s)>0)
		s.timesofar +=stopmstimer(s.timer)/10^6
		s.timer = startmstimer
		print  "Time : "+time2str(s.timesofar) +"  -   Error : System Creation failed no density matrices are found"
		return -1
	else
		s.timesofar +=stopmstimer(s.timer)/10^6
		s.timer = startmstimer
		Print "Time : "+time2str(s.timesofar) +"  -  Check of system density checks out"
	endif
	
	// make waves to store integrations and ratios
	make/d/o/n=(floor(s.num/sqrt(2)),energynum) int3DvsEn=0,ratio3DvsEn=0, para3dvsen, perp3dvsen
	make/d/o/n=(s.num,s.num,energynum) scatter3DSave
	wave s.scatter3DSave=scatter3DSave, s.int3DvsEn=int3DvsEn, s.ratio3DvsEn=ratio3DvsEn, s.para3dvsen=para3dvsen, s.perp3dvsen=perp3dvsen
	setscale /p x,pi/s.num,pi/s.num, s.int3DvsEn,s.ratio3DvsEn, s.para3dvsen, s.perp3dvsen // in physics units of q
	setscale /i x, -pi/s.voxelsize,pi/s.voxelsize, s.scatter3DSave
	setscale /i y, -pi/s.voxelsize,pi/s.voxelsize, s.scatter3DSave
	setscale /i y, energymin, energymax, s.int3DvsEn,s.ratio3DvsEn, s.para3dvsen, s.perp3dvsen
	setscale /i z, energymin, energymax, s.scatter3DSave
	
	
	variable en
	variable enstep = (energymax-energymin)/(energynum-1)
	s.step=0
	for( en = energymin ; en <=energymax ; en += enstep)
		s.timesofar +=stopmstimer(s.timer)/10^6
		s.timer = startmstimer
		s.en = en
		s.wavelength = 1239.8/en
		s.k = 2*pi/s.wavelength
		Print "Time : "+time2str(s.timesofar) +"  -  Creating Indexwave for energy : " + num2str(s.en)
		result = CalculatePolarizationWave(s)
		if(result<0)
			print "quitting because of critical errors in calculating polarization"
			return -1
		elseif(result>0 && s.en == energymin)
			print "warning - no aligned optical constants for "+stringfromlist(result-1, s.materials,",")+" were found, using normal Optical Constants"
		endif

		s.timesofar +=stopmstimer(s.timer)/10^6
		s.timer = startmstimer
		Print "Time : "+time2str(s.timesofar) +"  -  Creating 3D Scattering Pattern"
		Calculate3DScatter(s)
		RadialIntegratesystem(s)
		StoreIntegrations(s)
		// visualizescatter(s) // updates displays of scattering profiles
		//(re)create layout of scattering for an energy
		CreateSimLayout()
		// update the layout (projections, 2D scattering pattern, label, and cuts)
		TextBox/w=Simulation_Layout/C/N=text5/F=0/A=LB/X=33.80/Y=91.12 "\\F'Arial Black'\\Z24"+num2str(s.en)+"eV"
//		SetDrawLayer /w=RatioIntensity /k userfront
//		SetDrawEnv /w=RatioIntensity ycoord= left,linefgc= (65280,0,0),dash= 1,linethick= 3.00;DelayUpdate
//		DrawLine/w=RatioIntensity 0,s.en,1,s.en
//		SetDrawLayer /w=ScatteringIntensity /k userfront
//		SetDrawEnv /w=ScatteringIntensity ycoord= left,linefgc= (65280,0,0),dash= 1,linethick= 3.00;DelayUpdate
//		DrawLine/w=ScatteringIntensity 0,s.en,1,s.en
		doupdate
		if(s.movie)
			SavePict/O/WIN=Simulation_Layout /E=-5/P=_PictGallery_/w=(0,0,100*11,100*8) as "SimPict"
			AddMovieFrame/PICT=SimPict
		endif
		// add ratio and integration to 2D Q-energy plots
		s.step +=1
	endfor
	if(s.movie)
		s.timesofar +=stopmstimer(s.timer)/10^6
		s.timer = startmstimer
		Print "Time : "+time2str(s.timesofar) +"  -  Finishing Movie, and closing file"
		closemovie
	endif
	s.timesofar +=stopmstimer(s.timer)/10^6
	Print "Time : "+time2str(s.timesofar) +"  -  Complete."
	//display 2D rato and integrations vs q and energy
end
function createsimlayout()
	dowindow /k Simulation_Layout
	ScatteringIntensitydisp()
	RatioIntensitydisp()
	ScatteringPatterndisp()
	ParaPerpIntdisp()
	//DeltaProjectiondisp()
	//BetaProjectiondisp()
	execute("Simulation_Layout()")
	doupdate
end

function align3dsystem(s)
	struct ThreeDSystem &s
	string matstr = s.materialalignmentstring
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
			make /o/n=(dimsize(s.density1,0),dimsize(s.density1,1),dimsize(s.density1,2),4) s.m1=0
			s.m1[][][][3] = s.density1[p][q][r]
		else
			splitstring/e="^[^,]*,([^,]*),([^,]*)" mat, alignmentwidth, alignmentstrength
			strength = str2num(alignmentstrength)
			width = str2num(alignmentwidth)
			s.timesofar += stopmstimer(s.timer)/1e6
			s.timer = startmstimer
			wave alignmentw = createalignmentdensity(s.density1,width,strength, align, s.movie, s.timesofar)
			//createbinarysystem(s.density1,width,strength, 1, 0, 1- volfrac, 0, 1,s=s)
			duplicate/o alignmentw, s.m1
			//s.m1 *=s.density1[p][q][r]
			s.m1 *=sqrt(s.density1[p][q][r])
			s.m1[][][][3] *=sqrt(s.density1[p][q][r])
		endif
		mat = stringfromlist(1,matstr)
		if(stringmatch( mat,"*Face-on*")  || stringmatch(mat, "1*") )
			align = 1
		elseif(stringmatch( mat,"*Edge-on*") || stringmatch(mat, "0*") )
			align=0
		else
			duplicate/o s.m1, s.m2
			s.m2[][][][0,2]=0
			s.m2[][][][3] = 1 - s.m1[p][q][r][0]^2 - s.m1[p][q][r][1]^2 - s.m1[p][q][r][2]^2 -s.m1[p][q][r][3] // all density that isn't accounted for in m1 is unaligned m2
			align=-1
		endif
		if(align==-1)
			// we've already created both m1 and m2, so we're done!
		else
			splitstring/e="^[^,]*,([^,]*),([^,]*)" mat, alignmentwidth, alignmentstrength
			strength = str2num(alignmentstrength)
			width = str2num(alignmentwidth)
			volfrac = mean(s.density1)
			s.timesofar += stopmstimer(s.timer)/1e6
			s.timer = startmstimer
			wave alignmentw = createalignmentdensity(s.density1,width,strength, align, s.movie, s.timesofar)
		//	wave alignmentw = root:packages:ScatterSim3D:test:sn
			//createbinarysystem(s.density1,width,strength, 1, 0, 1- volfrac, 0, 1,s=s)
			duplicate/o alignmentw, s.m2
			s.m2 *=sqrt(1-s.density1[p][q][r])
			s.m2[][][][3] *=sqrt(1-s.density1[p][q][r])
		endif
	else
		for(i=0;i<nummaterials;i+=1)
			mat = stringfromlist(i,matstr)
			if(stringmatch("*Face-on*", mat))
				align = 1
			elseif(stringmatch("*Edge-on*", mat))
				align=0
			else
				continue
			endif
			splitstring/e="^[^,]*,([^,]*),([^,]*)" mat, alignmentwidth, alignmentstrength
			strength = str2num(alignmentstrength)
			width = str2num(alignmentwidth)
			if(strength * align * width * 0 !=0)
				continue // there is a nan, which will not work
			endif
			if(i==0)
				make /n=(dimsize(s.density1,0),dimsize(s.density1,1),dimsize(s.density1,2),4) s.m1
				wave materialwave = s.m1
			elseif(i==1)
				make /n=(dimsize(s.density2,0),dimsize(s.density2,1),dimsize(s.density2,2),4) s.m2
				wave materialwave = s.m2
			elseif(i==2)
				make /n=(dimsize(s.density3,0),dimsize(s.density3,1),dimsize(s.density3,2),4) s.m3
				wave materialwave = s.m3
			elseif(i==3)
				make /n=(dimsize(s.density4,0),dimsize(s.density4,1),dimsize(s.density4,2),4) s.m4
				wave materialwave = s.m4
			else
				continue
			endif
				
		//	createalignment(
		
		endfor
	endif
end
function storeintegrations(s)
	struct ThreeDSystem &s
	s.int3DvsEn[][s.step] =  s.int3d[p]
	s.ratio3DvsEn[][s.step] = (s.int3dpara(x) - s.int3dperp(x)) / (s.int3dpara(x) + s.int3dperp(x))
	s.scatter3dSave[][][s.step] = s.scatter3d(x)(y)
	s.perp3Dvsen[][s.step] =  s.int3dperp(x)
	s.para3Dvsen[][s.step] = s.int3dpara(x)
	
end
function radialintegratesystem(s)
	//radially integrate the 2D and 3D transforms and populate the corresponding waves in s (see structure declaration for details)
	struct ThreeDSystem &s
	variable ey =0// str2num(stringfromlist(0,s.efield,","))
	variable ez =1// str2num(stringfromlist(1,s.efield,","))
	variable pe = atan(ey/ez)*180/pi
	variable pa = atan(-ez/ey)*180/pi
	wave s.int3dp0 = radialintegratew(s.scatter3d,0,90,"int3dp0")
	wave s.int3dp0para = radialintegratew(s.scatter3d,pa-10,pa+10,"int3dp0para")
	wave s.int3dp0perp = radialintegratew(s.scatter3d,pe-10,pe+10,"int3dp0perp")
	duplicate/o s.int3dp0, s.int3d
	duplicate/o s.int3dp0para, s.int3dpara
	duplicate/o s.int3dp0perp, s.int3dperp
	s.int3d *= x^2
	s.int3dpara *= x^2
	s.int3dperp *= x^2
end



function Calculate3DScatter(s)
	struct ThreeDSystem &s
	variable thicknum = dimsize(s.px,0) , num = s.num
	if(num/2 != round(num/2) )
		num +=1
	endif
	if(thicknum/2 != round(thicknum/2) )
		thicknum +=1
	endif
	wave px = s.px, py=s.py, pz=s.pz
	FFT  /pad={1*thicknum,1*num,1*num} /Dest=pxfft px
	FFT  /pad={1*thicknum,1*num,1*num} /Dest=pyfft py
	FFT  /pad={1*thicknum,1*num,1*num} /Dest=pzfft pz
	wave/c s.pxfft = pxfft, s.pyfft=pyfft, s.pzfft=pzfft
	//make /n=(dimsize(s.pxfft,0),dimsize(s.pxfft,1),dimsize(s.pxfft,2),9) /d/o s.qtensor
	make /n=(dimsize(s.pxfft,0),dimsize(s.pxfft,1),dimsize(s.pxfft,2)) /d/o s.EscatSqr =0
	wave EscatSqr=s.EscatSqr
	setscale /i x, -pi/s.voxelsize, pi/s.voxelsize, pxfft, pyfft, pzfft, EscatSqr // changing from math units (1/voxelsize) to physics units (2pi/voxelsize)
	setscale /i y, -pi/s.voxelsize, pi/s.voxelsize, pxfft, pyfft, pzfft, EscatSqr
	setscale /i z, -pi/s.voxelsize, pi/s.voxelsize, pxfft, pyfft, pzfft, EscatSqr // all three dimensions have the same range.  The x dimension may have a lower number of pixels, but they cover the same range
//	s.qtensor[][][][0] = x * (2 * s.k + x)
//	s.qtensor[][][][1] = y * (s.k + x)
//	s.qtensor[][][][2] = z * (s.k + x)
//	s.qtensor[][][][3] = y * (s.k + x)
//	s.qtensor[][][][4] = y^2 - s.k^2
//	s.qtensor[][][][5] = y * z
//	s.qtensor[][][][6] = z * (s.k + x)
//	s.qtensor[][][][7] = y * z
//	s.qtensor[][][][8] = z^2 - s.k^2
	variable k = s.k
	multithread EscatSqr = magsqr(x * (2 * k - x) * pxfft[p][q][r] + y * (k - x) * pyfft[p][q][r] + z * (k - x) * pzfft[p][q][r])
	multithread EscatSqr += magsqr(y * (k - x)* pxfft[p][q][r]- (y^2 - k^2) * pyfft[p][q][r] - y * z * pzfft[p][q][r])
	multithread EscatSqr += magsqr(z * (k - x)* pxfft[p][q][r]- y * z* pyfft[p][q][r]- (z^2 - k^2 ) * pzfft[p][q][r])
	make /d/n=(s.num,s.num)/o s.scatter3D // this will hold the final scattering result
	wave scatter3d = s.scatter3d
	setscale /i x, -pi/s.voxelsize, pi/s.voxelsize, s.scatter3D // changing from math units (1/voxelsize) to physics units (2pi/voxelsize)
	setscale /i y, -pi/s.voxelsize, pi/s.voxelsize, s.scatter3D
	multithread scatter3D = k^2 > x^2 + y^2 ? interp3d(EscatSqr,-k + sqrt(k^2 - x^2 -y^2),x,y) : nan
	
	duplicate /o s.EscatSqr, pxreal
	//pxreal = magsqr(s.pxfft)
	imagetransform/meth=2 xprojection, pxreal
	duplicate/o m_xprojection, pxfftproj
	
	return 0
end

function CalculatePolarizationWave(s)
	struct ThreeDSystem &s
	variable result
	//create polarization waves
	make /o/c/d/n=(dimsize(s.m1,0),s.num,s.num) s.px = 0, s.py = 0, s.pz = 0
	setscale /p x, dimoffset(s.m1,0), dimdelta(s.m1,0), s.px, s.py, s.pz
	setscale /p y, dimoffset(s.m1,1), dimdelta(s.m1,1), s.px, s.py, s.pz
	setscale /p z, dimoffset(s.m1,2), dimdelta(s.m1,2), s.px, s.py, s.pz
	wave/c px= s.px, py=s.py, pz=s.pz
	//for each materials in list
	string material
	variable i, matnum = itemsinlist(s.materials,",")
	variable delpara, betapara, delperp, betaperp,ey,ez
	variable /c nperp, npara
	ey = 0 //str2num(stringfromlist(0,s.efield,","))
	ez = 1 //str2num(stringfromlist(1,s.efield,","))  // the complication of adding depolarization, makes the equations only work with incident polarization in the z direction
	if(ez+ey != 1)
		print "error - Electric field components are invalid : "+ s.efield
		return -1
	endif
	for(i=0;i<matnum;i+=1)
		material = stringfromlist(i,s.materials,",")
		//if there exist both para and perp optical constants load the corresponding values for that energy, otherwise throw warning and load same OC for both
		if(ocsexist(material)==2)
			delpara = getdelta(material,1239.8/s.en,aligned=1)
			betapara = getbeta(material,1239.8/s.en,aligned=1)
			delperp = getdelta(material,1239.8/s.en,aligned=0)
			betaperp = getbeta(material,1239.8/s.en,aligned=0)
			result = 0
		elseif(ocsexist(material)==1)
			delpara = getdelta(material,1239.8/s.en)
			betapara = getbeta(material,1239.8/s.en)
			delperp = getdelta(material,1239.8/s.en)
			betaperp = getbeta(material,1239.8/s.en)
			result = i+1
		else
			print "error - optical constants for " + material+ " were not found"
			return -1
		endif
		nperp = 1-delperp + sqrt(-1) * betaperp
		npara = 1-delpara + sqrt(-1) * betapara
		// select the right material polarized density wave
		if(i==0)
			wave m = s.m1
		elseif(i==1)
			wave m = s.m2
		elseif(i==2)
			wave m = s.m3
		elseif(i==3)
			wave m = s.m4
		elseif(i==4)
			wave m = s.m5
		else
			print " Error - invalid number of materials specified"
			return -1
		endif
		multithread px[][][] +=   ( (nperp^2 - npara^2)/(4*pi) ) * m[p][q][r][2]* m[p][q][r][0]
		multithread py[][][] +=   ( (nperp^2 - npara^2)/(4*pi) ) * m[p][q][r][2]* m[p][q][r][1]
		multithread pz[][][] +=   (m[p][q][r][0]^2 + m[p][q][r][1]^2 + m[p][q][r][2]^2 + m[p][q][r][3])/(4*pi) -  (m[p][q][r][3]/(12*pi)) * (npara^2 + 2*nperp^2) - (npara^2 * m[p][q][r][2]^2 + nperp^2 * (m[p][q][r][0]^2 +m[p][q][r][1]^2) )/(4*pi)
	endfor
	duplicate /o s.density1, pxreal
	pxreal = magsqr(s.px)
	imagetransform/meth=2 xprojection, pxreal
	duplicate/o m_xprojection, pxrealproj
	pxreal = magsqr(s.py)
	imagetransform/meth=2 xprojection, pxreal
	duplicate/o m_xprojection, pyrealproj
	pxreal = magsqr(s.pz)
	imagetransform/meth=2 xprojection, pxreal
	duplicate/o m_xprojection, pzrealproj
	return result
end

function sum3dsystem(s)
	struct ThreeDsystem &s
	make /o/n=(dimsize(s.m1,0),s.num,s.num) sumwave
	if(waveexists(s.m5))
		multithread sumwave = s.m1[p][q][r][0]^2 + s.m1[p][q][r][1]^2 + s.m1[p][q][r][2]^2 
		multithread sumwave += s.m2[p][q][r][0]^2 + s.m2[p][q][r][1]^2 + s.m2[p][q][r][2]^2 
		multithread sumwave += s.m3[p][q][r][0]^2 + s.m3[p][q][r][1]^2 +s.m3[p][q][r][2]^2  
		multithread sumwave += s.m4[p][q][r][0]^2 + s.m4[p][q][r][1]^2 +s.m4[p][q][r][2]^2  
		multithread sumwave += s.m5[p][q][r][0]^2 + s.m5[p][q][r][1]^2 +s.m5[p][q][r][2]^2 
		multithread sumwave += s.m1[p][q][r][3]+s.m2[p][q][r][3]+s.m3[p][q][r][3]+s.m4[p][q][r][3]+s.m5[p][q][r][3]
	elseif(waveexists(s.m4))
		multithread sumwave = s.m1[p][q][r][0]^2 + s.m1[p][q][r][1]^2 + s.m1[p][q][r][2]^2 
		multithread sumwave += s.m2[p][q][r][0]^2 + s.m2[p][q][r][1]^2 + s.m2[p][q][r][2]^2 
		multithread sumwave += s.m3[p][q][r][0]^2 + s.m3[p][q][r][1]^2 +s.m3[p][q][r][2]^2  
		multithread sumwave += s.m4[p][q][r][0]^2 + s.m4[p][q][r][1]^2 +s.m4[p][q][r][2]^2  
		multithread sumwave +=s.m1[p][q][r][3]+s.m2[p][q][r][3]+s.m3[p][q][r][3]+s.m4[p][q][r][3]
	elseif(waveexists(s.m3))
		multithread sumwave = s.m1[p][q][r][0]^2 + s.m1[p][q][r][1]^2 + s.m1[p][q][r][2]^2 + s.m2[p][q][r][0]^2 + s.m2[p][q][r][1]^2 + s.m2[p][q][r][2]^2 + s.m3[p][q][r][0]^2 + s.m3[p][q][r][1]^2 +s.m3[p][q][r][2]^2 +s.m1[p][q][r][3]+s.m2[p][q][r][3]+s.m3[p][q][r][3]
	elseif(waveexists(s.m2))
		multithread sumwave = s.m1[p][q][r][0]^2 + s.m1[p][q][r][1]^2 + s.m1[p][q][r][2]^2 + s.m2[p][q][r][0]^2 + s.m2[p][q][r][1]^2 + s.m2[p][q][r][2]^2 +s.m1[p][q][r][3]+s.m2[p][q][r][3]
	elseif(waveexists(s.m1))
		multithread sumwave = s.m1[p][q][r][0]^2 + s.m1[p][q][r][1]^2 + s.m1[p][q][r][2]^2 +s.m1[p][q][r][3]
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
function model3D_Spheres2(s)
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
	
	struct ThreeDSystem &s
	if(itemsinlist(s.paramstring,",")<5)
		return -1
	endif
	newdatafolder /o/s CreatingSpheres
	variable interpenetration = 	str2num(stringfromlist( 0 ,s.paramstring,","))
	variable minsep = 			str2num(stringfromlist( 1 ,s.paramstring,","))
	variable particlesnum = 	str2num(stringfromlist( 2 ,s.paramstring,","))
	variable pd = 				str2num(stringfromlist( 3 ,s.paramstring,","))
	variable thickness = 		s.thickness
	variable noise = 			str2num(stringfromlist( 4 ,s.paramstring,","))
	

	
	make /o /n=(thickness,s.num,s.num) mat=1,xwave, ywave, zwave
	
	make/B/U /o /n=(thickness,s.num,s.num,30) exmat= (p <= t) || (q <= t) || (r <= t) || (p >= thickness-t) || (q >= s.num-t) || (r >= s.num-t) ? 0 : 1
	make/B/U /o /n=(thickness,s.num,s.num) tempwave
	if(s.movie)
		Execute("Spheres3Ddisp(" +num2str(s.num)+", \""+getwavesdatafolder(mat,2)+"\")")
		Execute("exportgizmo wave as \"testimage\";Spinoidal3DLayout();Spinoidal3DImage(\""+getdatafolder(1)+"testimage\")")
	endif
	setscale /i x, -thickness/2, thickness/2, mat, exmat,xwave, ywave, zwave
	setscale /i y, -s.num/2, s.num/2, mat, exmat,xwave, ywave, zwave
	setscale /i z, -s.num/2, s.num/2, mat, exmat,xwave, ywave, zwave
	xwave = x
	ywave = y
	zwave = z
	redimension /n=(thickness*s.num*s.num) xwave, ywave, zwave
	variable i,radius, orad, cx, cy, cz, failed, fnum =0, xmn,xmx,ymn,ymx,zmn,zmx, loc
	
	for(i=0;i<particlesnum;i+=1)
		fnum=0
		do
			failed = 0
			radius = abs(gnoise(pd)+s.size)
			radius =radius < 1 ? 1 : radius
			if(minsep<1)
				orad = radius*(1+minsep/2)
			else
				orad = radius + minsep/2
			endif
			//duplicate/o /r=()()()(ceil(2*orad)) exmat,tempwave
			redimension /n=(thickness,s.num,s.num) tempwave
			multithread tempwave[][][] = exmat[p][q][r][ceil(orad)]
//			imagefilter /n=(ceil(2*orad)) /o min3d tempwave
//			multithread tempwave = tempwave<1? 0 : 1 
//			multithread tempwave = (p <= orad) || (q <= orad) || (r <= orad) || (p >= thickness-orad) || (q >= s.num-orad) || (r >= s.num-orad) ? 0 : tempwave[p][q][r]
  			if(wavemax(tempwave)<1)
				//there are no possible locations for this radius, find another
				failed=1
			else
				// randomly pick a pixel that is good for the center
				redimension /n=(thickness*s.num*s.num) tempwave
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
		if(s.movie)
			execute("ModifyGizmo /n=Spheres3D update=2")
			doupdate
			Execute "exportgizmo wave as \"testimage\""
			//TextBox/w=Spinoidal3DLayout/C/N=text0/A=LT/X=0.00/Y=0.00 "\Z32" + time2str2(ttot)
			doupdate
			savepict /p=_PictGallery_ /E=-5 /N=Spinoidal3DLayout /o as "Frame3D"
			addmovieframe /pict=Frame3D
		endif
	endfor
	setdatafolder ::
	imagefilter /n=(interpenetration)/o gauss3d mat
	duplicate /o mat,s.density1 // this returns the density matrix of material 1 (the matrix) for alignment etc later on
end

function loadexistingmaterialalignment(s)
	// loads existing model into memory for the rest of the calculation
	// no parameter wave is required
	struct ThreeDSystem &s
	wave/z s.m1 = m1
	wave/z s.m2 = m2
	wave/z s.m3 = m3
	wave/z s.m4 = m4
	wave/z s.m5 = m5
end

function/s variables_existing()
	return ""
end
function/s special_existing()
	return ""
end
function model3D_existing(s)
	// loads existing model into memory for the rest of the calculation
	// no parameter wave is required
	struct ThreeDSystem &s
	wave/z s.density1 = density1
	if(!waveexists(s.density1))
		print "Density1 wave was not found"
		return -1
	endif
	wave/z s.density2 = density2
	if(!waveexists(s.density1))
		print "Warning : Density2 wave was not found - If the number of materials is less than 3, this is fine :)"
		return 1
	endif
	wave/z s.density3 = density3
	if(!waveexists(s.density1))
		print "Warning : Density3 wave was not found - If the number of materials is less than 4, this is fine :)"
		return 1
	endif
	wave/z s.density4 = density4
	if(!waveexists(s.density1))
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
	PauseUpdate; Silent 1		// building window...
	Layout/P=Landscape/C=1/K=1/W=(134.25,68.75,602.25,440.75) ScatteringIntensity(378,54,576,594)/O=1/F=0
	PrintSettings /i /w=Simulation_Layout margins={.25,.25,.25,.25}
	Append RatioIntensity(576,54,774,594)/O=1/F=0//,BetaProjection(24.75,194.25,192.75,345.75)/O=1/F=0
	//Append DeltaProjection(198.75,194.25,368.25,350.25)/O=1/F=0
	Append ScatteringPattern(39.75,29.25,252.75,189.75)/O=1/F=0
	Append ParaPerpInt(11.25,339.75,377.25,584.25)/O=1/F=0
	TextBox/C/N=text0/F=0/A=LB/X=28.63/Y=66.06 "Real part of Index of\rRefraction (projection)"
	TextBox/C/N=text1/F=0/A=LB/X=4.77/Y=66.32 "Imaginary part of Index of\rRefraction (projection)"
	TextBox/C/N=text2/F=0/A=LB/X=76.74/Y=93.21 "\\Z16Anisotropic Ratio \rvs Q and Energy"
	TextBox/C/N=text3/F=0/A=LB/X=51.69/Y=93.21 "\\Z16Total Scatter (*q^2) \rvs Q and Energy"
	TextBox/C/N=text4/F=0/A=LB/X=1.39/Y=97.13 "\\Z16Simulated Scattering Pattern"
	TextBox/C/N=text5/F=0/A=LB/X=33.80/Y=91.12 "\\F'Arial Black'\\Z24286eV"
	SetDrawLayer UserFront
	DrawLine 0,291,1,297
EndMacro

function ScatteringIntensitydisp()
	dowindow /k ScatteringIntensity
	Display/n=ScatteringIntensity /W=(819,176,991.5,665.75)/K=1 
	AppendImage/w=ScatteringIntensity /T int3DvsEn
	ModifyImage/w=ScatteringIntensity  int3DvsEn ctab= {1e-05,100,BlueHot,0}
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
	SetAxis/w=ScatteringIntensity /R left 320,260
	SetAxis/w=ScatteringIntensity /A  top 
	ColorScale/w=ScatteringIntensity /C/N=text0/A=MC/X=54.43/Y=-0.16 image=int3DvsEn, side=2, width=10
	ColorScale/w=ScatteringIntensity /C/N=text0 log=1
	SetDrawLayer/w=ScatteringIntensity  UserFront
	SetDrawEnv/w=ScatteringIntensity  xcoord= prel,ycoord= left,linethick= 3,linefgc= (65280,0,0),dash= 1
	DrawLine/w=ScatteringIntensity  0,319.999999999984,1,319.999999999984
End

function RatioIntensitydisp()
	dowindow /k RatioIntensity
	Display/n=RatioIntensity /W=(1002.75,176.75,1175.25,664.25)/K=1 
	AppendImage/w=RatioIntensity /T ratio3DvsEn
	ModifyImage/w=RatioIntensity  ratio3DvsEn ctab= {-1,1,RedWhiteGreen,0}
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
	SetAxis/w=RatioIntensity /R left 320,260
	SetAxis/w=RatioIntensity /A top
	ColorScale/w=RatioIntensity /C/N=text0/A=MC/X=56.33/Y=0.66 image=ratio3DvsEn, side=2, width=10
	SetDrawLayer/w=RatioIntensity  UserFront
	SetDrawEnv/w=RatioIntensity  xcoord= prel,ycoord= left,linethick= 3,linefgc= (65280,0,0),dash= 1
	DrawLine/w=RatioIntensity  0,319.999999999984,1,319.999999999984
End

function ScatteringPatterndisp()
	dowindow /k ScatteringPattern
	Display/n=ScatteringPattern /W=(1151.25,48.5,1361.25,209)/K=1 
	AppendImage/w=ScatteringPattern /T scatter3D
	ModifyImage/w=ScatteringPattern  scatter3D ctab= {0.0001,30000,Terrain,0}
	ModifyImage/w=ScatteringPattern  scatter3D ctabAutoscale=1
	ModifyImage/w=ScatteringPattern  scatter3D log= 1
	ModifyGraph/w=ScatteringPattern  margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=58
	ModifyGraph/w=ScatteringPattern  mirror=2
	ModifyGraph/w=ScatteringPattern  nticks=4
	ModifyGraph/w=ScatteringPattern  minor=1
	ModifyGraph/w=ScatteringPattern  fSize=8
	ModifyGraph/w=ScatteringPattern  standoff=0
	ModifyGraph/w=ScatteringPattern  tkLblRot(left)=90
	ModifyGraph/w=ScatteringPattern  btLen=3
	ModifyGraph/w=ScatteringPattern  tlOffset=-2
	SetAxis/w=ScatteringPattern /R left 0.5,-0.5
	SetAxis/w=ScatteringPattern  top -0.5,0.5
	ColorScale/w=ScatteringPattern /C/N=text0/F=0/A=MC/X=69.19/Y=2.81 image=scatter3D, log=1
End

function BetaProjectiondisp()
	dowindow/k BetaProjection
	Display/n=BetaProjection /W=(1130.25,433.25,1301.25,594.5)/K=1 
	AppendImage/w=BetaProjection /T npi
	ModifyImage/w=BetaProjection  npi ctab= {0,0.003,BlueHot,0}
	ModifyGraph/w=BetaProjection  margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=22
	ModifyGraph/w=BetaProjection  mirror=2
	ModifyGraph/w=BetaProjection  nticks=3
	ModifyGraph/w=BetaProjection  minor=1
	ModifyGraph/w=BetaProjection  fSize=8
	ModifyGraph/w=BetaProjection  standoff=0
	ModifyGraph/w=BetaProjection  tkLblRot(left)=90
	ModifyGraph/w=BetaProjection  btLen=3
	ModifyGraph/w=BetaProjection  tlOffset=-2
	SetAxis/w=BetaProjection /A/R left
	setaxis /w=BetaProjection left 0,100
	setaxis /w=BetaProjection top 0,100
	ColorScale/w=BetaProjection /C/N=text0/F=0/A=MC/X=44.20/Y=6.15 image=npi, nticks=2, lowTrip=1e-05
End

function ParaPerpIntdisp()
	dowindow/k ParaPerpInt
	Display/n=ParaPerpInt /W=(798,43.25,1173,213.5)/K=1  int3Dpara,int3Dperp
	ModifyGraph/w=ParaPerpInt lSize=2
	ModifyGraph/w=ParaPerpInt  rgb(int3Dpara)=(65280,16384,16384),rgb(int3Dperp)=(5808,0,5808)
	ModifyGraph/w=ParaPerpInt  log=1
	SetAxis/w=ParaPerpInt  left 1e-06,100
	Legend/w=ParaPerpInt /C/N=text0/J/F=0/A=MC/X=22.36/Y=-38.42 "\\Z14\\s(int3Dpara) Parallel with E-Field\r\\s(int3Dperp) Perpendicular to E-Field"
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
	AppendImage/w=deltaprojection/T npr
	ModifyImage/w=deltaprojection npr ctab= {0.997,1.003,BlueHot,0}
	ModifyGraph/w=deltaprojection margin(left)=14,margin(bottom)=14,margin(top)=14,margin(right)=22
	ModifyGraph/w=deltaprojection mirror=2
	ModifyGraph/w=deltaprojection nticks=3
	ModifyGraph/w=deltaprojection minor=1
	ModifyGraph/w=deltaprojection fSize=8
	ModifyGraph/w=deltaprojection standoff=0
	ModifyGraph/w=deltaprojection tkLblRot(left)=90
	ModifyGraph/w=deltaprojection btLen=3
	ModifyGraph/w=deltaprojection tlOffset=-2
	SetAxis/A/R/w=deltaprojection left
	setaxis /w=deltaprojection left 0,100
	setaxis /w=deltaprojection top 0,100
	ColorScale/w=deltaprojection/C/N=text0/F=0/A=MC/X=41.90/Y=5.59 image=npr, nticks=2, minor=1
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
			setvariable setefield,disable= (tab!=2)
			CheckBox SaveMovieChk,disable= (tab!=2)
			CheckBox Save2DDataEn,disable= (tab!=2)
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
		nvar sizescale, thickness, resolution, startxrayenergy, endxrayenergy, numensteps, movie, voxelsize, usealignment
		nvar alignmentsize, alignmentstrength, useprecalcalignment,SaveScattering1D, SaveParaPerp1D, SaveAnisotropic1D, SaveScattering2D
	else
		string/g funcnames="", morphology="existing", material="pcbm", alignment="None", SimName = "test", efield = " ( 1.0000 , 0.0000 )"
		string /g controllist="", extralist=""
		variable /g sizescale=5, thickness=16, resolution=128, startxrayenergy=260, endxrayenergy=320, numensteps=600, movie=0, voxelsize=5
		variable /g alignmentsize=10, alignmentstrength=1, useprecalcalignment=0, usemorphalignment=0
		variable /g  SaveScattering1D=1, SaveParaPerp1D=0, SaveAnisotropic1D=1, SaveScattering2D=0
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
	SetVariable SetEnSteps,limits={10,30000,0.05},value= root:Packages:ScatterSim3D:numensteps,live= 1
	GroupBox MorphOptions,pos={20,87},size={431,314},title="Variables for Morphology"
	Button RemoveMaterial,pos={40,216},size={134,39},disable=1,proc=RemoveMaterialButton,title="Remove Material"
	TabControl tab0,pos={8,7},size={458,443},proc=Model3dTabProc,tabLabel(0)="Morphology"
	TabControl tab0,tabLabel(1)="Materials",tabLabel(2)="Scattering",value= 0
	GroupBox MaterialsListText,pos={206,37},size={249,206},disable=1,title="List of Materials"
	SetVariable SetRes,pos={20,31},size={258,16},title="Resolution (Film Width and Length) [Pixels]"
	SetVariable SetRes,limits={10,5000,1},value= root:Packages:ScatterSim3D:resolution,live= 1
	SetVariable SetSize,pos={180,60},size={137,16},title="Length Scale [Voxels]"
	SetVariable SetSize,limits={1,5000,1},value= root:packages:ScatterSim3D:sizescale,live= 1
	Button CalcMorphBut,pos={262,409},size={192,31},title="Calculate Morphology Now"
	Button CalcAlignBut,pos={206,250},size={243,38},disable=1,title="Calculate Material Alignment Now"
	CheckBox SaveMovieChk,pos={234,54},size={75,14},disable=1,title="Save Movie"
	CheckBox SaveMovieChk,variable= root:Packages:ScatterSim3D:movie
	CheckBox SaveQEn,pos={234,72},size={161,14},disable=1,title="Save 1D Scattering vs Energy"
	CheckBox SaveQEn,variable= root:Packages:ScatterSim3D:SaveScattering1D
	CheckBox Save2DDataEn,pos={234,121},size={187,14},disable=1,title="Save 2D Scattering Data vs Energy"
	CheckBox Save2DDataEn,variable= root:Packages:ScatterSim3D:SaveScattering2D
	CheckBox SaveParaPerp,pos={234,88},size={206,14},disable=1,title="Save Parallel and Perpindicular 1D data"
	CheckBox SaveParaPerp,variable= root:Packages:ScatterSim3D:SaveParaPerp1D
	CheckBox SaveAnisotropy,pos={234,104},size={177,14},disable=1,title="Save Anisotropic Ratio Vs Energy"
	CheckBox SaveAnisotropy,variable= root:Packages:ScatterSim3D:SaveAnisotropic1D
	SetVariable SetThickness,pos={291,33},size={146,16},title="Thickness (Pixels)"
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
function model3d_spinoidal(s)
	struct ThreeDSystem &s
	variable ndim = s.num
	variable thickness = s.thickness
	variable num = str2num(stringfromlist(0,s.paramstring,","))
	variable eps = str2num(stringfromlist(1,s.paramstring,","))
	variable texp = str2num(stringfromlist(2,s.paramstring,","))
	variable t0 = str2num(stringfromlist(3,s.paramstring,","))
//	setdatafolder root:packages:ScatterSim3D
	make /o/n=(s.thickness,s.num,s.num) /o s.density1
	newdatafolder /o/s modelCreation
	spinoidalLB3Dn(ndim, num,eps, thickness,texp,t0,movie = s.movie)
	wave realu
	s.density1 = realu[q][r][p] // density's smaller dimension is p, whereas realu's smallest dimension is layers
	killwaves /A/Z
	setdatafolder ::
	s.density1 = s.density1 > 0 ? 1 : 0
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
		Execute("exportgizmo wave as \"testimage\";Spinoidal3DLayout_2();Spinoidal3DImage(\""+getdatafolder(1)+"testimage\")")
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
			Execute "exportgizmo wave as \"testimage\""
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
window fibrils3Ddisp(size, densitywaveloc) : GizmoPlot
	variable size
	string densitywaveloc
	PauseUpdate; Silent 1	// Building Gizmo 6 window...
	dowindow /k fibrils3D
	// Do nothing if the Gizmo XOP is not available.
	if(exists("NewGizmo")!=4)
		DoAlert 0, "Gizmo XOP must be installed"
		return 0
	endif

	NewGizmo/N=fibrils3D/T="3D fibrils" /W=(30,30,1310,1054) /k=1
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
//			nvar useprecalcalignment
			svar efield
			string ycomp, zcomp
			splitstring /e="^[( ]*([1234567890.]*)[^,]{0,5},[^,0.1]{0,5}([1234567890.]*)[) ]*" efield, ycomp, zcomp
			newdatafolder /o/s $cleanupname(simname,0)
			model3D(morphology,voxelsize,sizescale,resolution,thickness,paramstring,materialstring,alignmentstring,ycomp+","+zcomp,startxrayenergy,endxrayenergy,numensteps,movie = movie, save3d = savescattering2D)
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
function /s variables_fibrils()
	string variables
	variables = "Interpenetration [voxels],SetVariable,2;"
	variables = variables + "Maximum umber of fibrils,SetVariable,500;"
	variables = variables + "Poly dispursity (radii),SetVariable,.1;"
	variables = variables + "Absolute minimum radius,SetVariable,2.5;"
	variables = variables + "Absolute maximum radius,SetVariable,3.5;"
	variables = variables + "Sigma (around x-plane in degrees),SetVariable,5;"
	variables = variables + "Material volume fraction (<1),SetVariable,.5;"
	variables = variables + "Noise (not used),SetVariable,0;"
	variables = variables + "Minimum seperation,SetVariable,0.1;"
	variables = variables + "Minimum fibril length,SetVariable,10;"
	variables = variables + "Maximum fibril length,SetVariable,100;"
	variables = variables + "Goal fibril volume fraction (<1),SetVariable,.1;"
	variables = variables + "Fibril relative density,SetVariable,1;"
	variables = variables + "Degree of alignment (0-1) ,SetVariable,0;"
	return variables
end
function /s special_fibrils()
	return "IncludesAlignment"
end
function model3D_fibrils(s)

	struct ThreeDSystem &s
	if(itemsinlist(s.paramstring,",")<11)
		return -1
	endif
	newdatafolder /o/s Creatingfibrils
	variable interpenetration = 	str2num(stringfromlist( 0 ,s.paramstring,","))
	variable fibrilnum = 		str2num(stringfromlist( 1 ,s.paramstring,","))
	variable pd = 				str2num(stringfromlist( 2 ,s.paramstring,","))
	variable thickness = 		s.thickness
	variable minsize = 		str2num(stringfromlist( 3 ,s.paramstring,","))
	variable maxsize = 		str2num(stringfromlist( 4 ,s.paramstring,","))
	variable angsigma = 		str2num(stringfromlist( 5 ,s.paramstring,","))
	variable volfrac = 			str2num(stringfromlist( 6 ,s.paramstring,","))
	variable noise = 			str2num(stringfromlist( 7 ,s.paramstring,","))
	variable minsep = 		str2num(stringfromlist( 8 ,s.paramstring,","))
	variable minlength = 		str2num(stringfromlist( 9 ,s.paramstring,","))
	variable maxlength = 		str2num(stringfromlist( 10 ,s.paramstring,","))
	variable Volfracgoal = 		str2num(stringfromlist( 11 ,s.paramstring,","))
	variable FibRelDensity = 	str2num(stringfromlist( 12 ,s.paramstring,","))
	variable DegAlign = 		str2num(stringfromlist( 13 ,s.paramstring,","))
	
	make /n=(thickness*s.num^2,3)/o rmat
	make /o /n=(thickness,s.num,s.num) mat=0,xwave, ywave, zwave
	make /n=( thickness,s.num,s.num,3)/o vecmat=0
	if(minsep<1)
		make/B/U /o /n=(thickness,s.num,s.num,ceil(maxsize+minsep/2)) exmat =  (p <= t/(1+minsep)) || (q <= t/(1+minsep)) || (r <= t/(1+minsep) ) || (p >= thickness-t/(1+minsep)) || (q >= s.num-t/(1+minsep)) || (r >= s.num-t/(1+minsep)) ? 0 : 1
	else
		make/B/U /o /n=(thickness,s.num,s.num,ceil(maxsize+minsep/2)) exmat =  (p <= t-minsep) || (q <= t-minsep) || (r <= t-minsep) || (p >= thickness-t+minsep) || (q >= s.num-t+minsep) || (r >= s.num-t+minsep) ? 0 : 1
	endif
	make/B/U /o /n=(thickness,s.num,s.num) tempwave, tempx
	if(s.movie)
		Execute("fibrils3Ddisp(" +num2str(s.num)+", \""+getwavesdatafolder(mat,2)+"\")")
		Execute("exportgizmo wave as \"testimage\";Spinoidal3DLayout();Spinoidal3DImage(\""+getdatafolder(1)+"testimage\")")
	endif

	xwave = x
	ywave = y
	zwave = z
	redimension /n=(thickness*s.num*s.num) xwave, ywave, zwave
	rmat[][0] = xwave[p]
	rmat[][1] = ywave[p]
	rmat[][2] = zwave[p]
	variable testcx,testcy,testcz,i,radius, orad, cx, cy, cz, failed, fnum =0, f2num=0, xmn,xmx,ymn,ymx,zmn,zmx, loc, qfnum=0, theta, phi
	variable tx,ty,tz,hit, mx=thickness-1, my=s.num-1, mz=s.num-1, length, totlength
	fnum=0
	make /o /n=3 vec
	variable fibrilvol
	for(i=0;i<fibrilnum;i+=1)
		radius=min(maxsize,max(minsize,gnoise(pd)+s.size))
		tempx = exmat[p][q][r][ceil(radius +minsep/2)]
		duplicate /o tempx, tempwave
		redimension /n=(numpnts(tempwave)) tempwave
		integrate tempwave /D=intwave
		if(wavemax(intwave) < 5)
			fnum+=1
			if(fnum<5)
				continue
			else
				print "warning:  can't fit in anymore fibrils, only " + num2str(i-fnum) + " fibrils were created"
				break
			endif
		endif
		loc = binarysearch(intwave, enoise(wavemax(intwave)/2)+wavemax(intwave)/2)
		
		cx = rmat[loc][0]
		cy = rmat[loc][1]
		cz = rmat[loc][2]
		
		theta = (90 + gnoise(angsigma)) * pi / 180
		
		if(DegAlign<=0)
			phi = enoise(pi)
		else
			phi = gnoise(.1/DegAlign)
		endif
		vec[2] = sin(theta)*cos(phi)
		vec[1] = sin(theta)*sin(phi)
		vec[0] = cos(theta)
		
		//forward direction
		hit=0
		length=0
		tx=cx
		ty=cy
		tz=cz
		do
			vecmat[max(tx-radius,0),min(tx+radius,mx)][max(ty-radius,0),min(ty+radius,my)][max(tz-radius,0),min(tz+radius,mz)][]= (p-tx)^2 + (q-ty)^2 + (r-tz)^2 < radius^2 ? vec[t] : vecmat
			exmat[max(tx-radius-maxsize,0),min(tx+radius+maxsize,mx)][max(ty-radius-maxsize,0),min(ty+radius+maxsize,my)][max(tz-radius-maxsize,0),min(tz+radius+maxsize,mz)][]= (p-tx)^2 + (q-ty)^2 + (r-tz)^2 < (radius+t+minsep/2)^2 ? 0 : exmat // added in minsep/2
			
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
		totlength = length
		length-=1
		do
			vecmat[max(tx-radius,0),min(tx+radius,mx)][max(ty-radius,0),min(ty+radius,my)][max(tz-radius,0),min(tz+radius,mz)][]= (p-tx)^2 + (q-ty)^2 + (r-tz)^2 < radius^2 ? vec[t] : vecmat
			exmat[max(tx-radius-maxsize,0),min(tx+radius+maxsize,mx)][max(ty-radius-maxsize,0),min(ty+radius+maxsize,my)][max(tz-radius-maxsize,0),min(tz+radius+maxsize,mz)][]= (p-tx)^2 + (q-ty)^2 + (r-tz)^2 < (radius+t+minsep/2)^2 ? 0 : exmat
			
			tx -= vec[0]
			ty -= vec[1]
			tz -= vec[2]
			hit = 1-tempx[Min(max(tx,0),mx)][Min(max(ty,0),my)][Min(max(tz,0),mz)] // if we hit another fibril or the edge, tempx will be 0 (hit will be 1)  -  until then, hit will equal 0
			length+=1
		while( length<maxlength && hit==0 && tx <= mx && ty <= my && tz <= mz && tx>=0 && ty>=0 && tz>=0 )
		totlength +=length
		mat = sqrt( vecmat[p][q][r][0]^2 + vecmat[p][q][r][1]^2 + vecmat[p][q][r][2]^2 )
		fibrilvol = mean(mat)
		imagefilter /n=3 /o gauss3d, mat
		if(s.movie)
			execute("ModifyGizmo /n=fibrils3D update=2")
			doupdate
			Execute "exportgizmo wave as \"testimage\""
			//TextBox/w=Spinoidal3DLayout/C/N=text0/A=LT/X=0.00/Y=0.00 "\Z32" + time2str2(ttot)
			doupdate
			savepict /p=_PictGallery_ /E=-5 /N=Spinoidal3DLayout /o as "Frame3D"
			addmovieframe /pict=Frame3D
		endif
		
		if(totlength < minlength)
			print num2str(i) + " fibrils created, volume fraction = "+num2str(fibrilvol) + ", but this is fibril number "+num2str(f2num+1)+" that is too short at "+num2str(totlength)+ " voxels in length. When num/10 + 5 fibrils are too short, creation will cease."
			f2num+=1
		else
			print num2str(i) + " fibrils created, volume fraction = "+num2str(fibrilvol)
		endif
		
		if(f2num>(5 + i/10)) // if ever more than 10% + 5 are short, then stop
			print "too many short fibrils have been created, ending fibril creation"
			break
		endif
		if(fibrilvol > Volfracgoal)
			Print "Volume Fraction goal met"
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
	mat = sqrt( vecmat[p][q][r][0]^2 + vecmat[p][q][r][1]^2 + vecmat[p][q][r][2]^2 )
	
	setdatafolder ::
	fibrilvol = mean(mat)
	print "Volume fraction of "+num2str(fibrilvol)+" fibrils created"
	variable rhomatrix = (volfrac - fibrilvol)/(1-fibrilvol)
	if(rhomatrix < 0 )
		print "Volume fraction is too low, make less fibrils"
		return -1
	endif
	make /n=(thickness,s.num,s.num,4) /o m1=0, m2=0
	wave s.m1=m1, s.m2=m2
	s.m1[][][][0,2] = vecmat[p][q][r][t] * sqrt(FibRelDensity) // sets the aligned parts (fibril parts)  to this density fraction
	
	s.m1[][][][3] = rhomatrix * (1-mat[p][q][r]^2)
	s.m2[][][][3] = (1-rhomatrix) * (1-mat[p][q][r]^2)
	
	
	duplicate /o mat,s.density1 // this returns the density matrix of material 1 (the matrix) for alignment etc later on
end
function /s variables_spheresln()
	string variables
//	return "Interpenetration [pixels],SetVariable,2;Minimum Seperation,SetVariable,.1;Number of Particles (Max),SetVariable,500;PolyDispursity (sigma of radiuses),setvariable,5;Maximum Radius,SetVariable,8;Noise,SetVariable,0;"
	return "Interpenetration [pixels],SetVariable,2;Minimum Seperation,SetVariable,.1;Number of Particles (Max),SetVariable,500;PolyDispursity (sigma of radiuses),setvariable,5;Minimum Radius,SetVariable,1;Maximum Radius,SetVariable,8;Volume Fraction (<1),SetVariable,.65;Noise,SetVariable,0;"
end
function model3D_Spheresln(s)
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
	
	struct ThreeDSystem &s
	if(itemsinlist(s.paramstring,",")<5)
		return -1
	endif
	newdatafolder /o/s CreatingSpheres
	variable interpenetration = 	str2num(stringfromlist( 0 ,s.paramstring,","))
	variable minsep = 			str2num(stringfromlist( 1 ,s.paramstring,","))
	variable particlesnum = 	str2num(stringfromlist( 2 ,s.paramstring,","))
	variable pd = 				str2num(stringfromlist( 3 ,s.paramstring,","))
	variable thickness = 		s.thickness
	variable minsize = 		str2num(stringfromlist( 4 ,s.paramstring,","))
	variable maxsize = 		str2num(stringfromlist( 5 ,s.paramstring,","))
	variable volfrac = 			str2num(stringfromlist( 6 ,s.paramstring,","))
	variable noise = 			str2num(stringfromlist( 7 ,s.paramstring,","))
	
	make /o /n=(thickness,s.num,s.num) mat=1,xwave, ywave, zwave
	if(minsep<1)
		make/B/U /o /n=(thickness,s.num,s.num,maxsize) exmat =  (p <= t/(1+minsep)) || (q <= t/(1+minsep)) || (r <= t/(1+minsep) ) || (p >= thickness-t/(1+minsep)) || (q >= s.num-t/(1+minsep)) || (r >= s.num-t/(1+minsep)) ? 0 : 1
	else
		make/B/U /o /n=(thickness,s.num,s.num,maxsize) exmat =  (p <= t-minsep) || (q <= t-minsep) || (r <= t-minsep) || (p >= thickness-t+minsep) || (q >= s.num-t+minsep) || (r >= s.num-t+minsep) ? 0 : 1
	endif
	make/B/U /o /n=(thickness,s.num,s.num) tempwave
	if(s.movie)
		Execute("Spheres3Ddisp(" +num2str(s.num)+", \""+getwavesdatafolder(mat,2)+"\")")
		Execute("exportgizmo wave as \"testimage\";Spinoidal3DLayout();Spinoidal3DImage(\""+getdatafolder(1)+"testimage\")")
	endif
	setscale /i x, -thickness/2, thickness/2, mat, exmat,xwave, ywave, zwave
	setscale /i y, -s.num/2, s.num/2, mat, exmat,xwave, ywave, zwave
	setscale /i z, -s.num/2, s.num/2, mat, exmat,xwave, ywave, zwave
	xwave = x
	ywave = y
	zwave = z
	redimension /n=(thickness*s.num*s.num) xwave, ywave, zwave
	variable testcx,testcy,testcz,i,radius, orad, cx, cy, cz, failed, fnum =0, xmn,xmx,ymn,ymx,zmn,zmx, loc, qfnum=0
	fnum=0
	if(minsep<1)
		orad = s.size*(1+minsep/2)
	else
		orad = s.size + minsep/2
	endif
	wave radiuses = getradiusdistributionln(orad,s.num,thickness,pd,minsep,minsize,maxsize,volfrac)
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
			testcy = enoise((s.num-orad)/2)
			testcz = enoise((s.num-orad)/2)
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
			//redimension /n=(thickness,s.num,s.num) tempwave
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
			//redimension /n=(thickness*s.num*s.num) tempwave
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
		if(s.movie)
			execute("ModifyGizmo /n=Spheres3D update=2")
			doupdate
			Execute "exportgizmo wave as \"testimage\""
			//TextBox/w=Spinoidal3DLayout/C/N=text0/A=LT/X=0.00/Y=0.00 "\Z32" + time2str2(ttot)
			doupdate
			savepict /p=_PictGallery_ /E=-5 /N=Spinoidal3DLayout /o as "Frame3D"
			addmovieframe /pict=Frame3D
		endif
	endfor
	setdatafolder ::
	imagefilter /n=(interpenetration)/o gauss3d mat
	duplicate /o mat,s.density1 // this returns the density matrix of material 1 (the matrix) for alignment etc later on
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
function model3D_Cylsln(s)
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
	
	struct ThreeDSystem &s
	if(itemsinlist(s.paramstring,",")<5)
		return -1
	endif
	newdatafolder /o/s CreatingSpheres
	variable interpenetration = 	str2num(stringfromlist( 0 ,s.paramstring,","))
	variable minsep = 			str2num(stringfromlist( 1 ,s.paramstring,","))
	variable particlesnum = 	str2num(stringfromlist( 2 ,s.paramstring,","))
	variable pd = 				str2num(stringfromlist( 3 ,s.paramstring,","))
	variable thickness = 		s.thickness
	variable minsize = 		str2num(stringfromlist( 4 ,s.paramstring,","))
	variable maxsize = 		str2num(stringfromlist( 5 ,s.paramstring,","))
	variable volfrac = 			str2num(stringfromlist( 6 ,s.paramstring,","))
	variable noise = 			str2num(stringfromlist( 7 ,s.paramstring,","))
	

	
	make /o /n=(s.num,s.num) mat=1, ywave, zwave
	if(minsep<1)
		make/B/U /o /n=(s.num,s.num,30) exmat =  (p <= r/(1+minsep)) || (q <= r/(1+minsep) ) || (p >= s.num-r/(1+minsep)) || (q >= s.num-r/(1+minsep)) ? 0 : 1
	else
		make/B/U /o /n=(s.num,s.num,30) exmat =   (p <= r-minsep) || (q <= r-minsep) || (p >= s.num-r+minsep) || (q >= s.num-r+minsep) ? 0 : 1
	endif
	make/B/U /o /n=(s.num,s.num) tempwave
	if(s.movie)
		//Execute("Spheres3Ddisp(" +num2str(s.num)+", \""+getwavesdatafolder(mat,2)+"\")")
		//Execute("exportgizmo wave as \"testimage\";Spinoidal3DLayout();Spinoidal3DImage(\""+getdatafolder(1)+"testimage\")")
		newimage /k=1 /n=cylimage mat
	endif

	setscale /i y, -s.num/2, s.num/2, mat, exmat, ywave, zwave
	setscale /i z, -s.num/2, s.num/2, mat, exmat, ywave, zwave
	
	ywave = x// 2D y-> x and z->y
	zwave = y
	redimension /n=(s.num*s.num) ywave, zwave
	variable i,radius, orad, cx, cy, cz, failed, fnum =0, xmn,xmx,ymn,ymx,zmn,zmx, loc
	fnum=0
	if(minsep<1)
		orad = s.size*(1+minsep/2)
	else
		orad = s.size + minsep/2
	endif
	wave radiuses = get2dradiusdistributionln(orad,s.num,pd,minsep,minsize,maxsize,volfrac)
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
		redimension /n=(s.num,s.num) tempwave
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
		redimension /n=(s.num*s.num) tempwave
		integrate tempwave /D=intwave
		loc = binarysearch(intwave, enoise(wavemax(intwave)/2)+wavemax(intwave)/2)
		cx = ywave[loc]
		cy = zwave[loc]
		
		// subtract out this sphere from the matrix  // matrix starts at 1s, within this sphere, multiply this by 0, outside multiply by 1
		multithread mat*= (x-cx)^2 + (y-cy)^2 < radius^2 ? 0 : 1 
		multithread exmat*= (x-cx)^2 + (y-cy)^2 <= (orad+r)^2 ? 0 : 1 
		if(s.movie)
			//execute("ModifyGizmo /n=Spheres3D update=2")
			doupdate
			//Execute "exportgizmo wave as \"testimage\""
			//TextBox/w=Spinoidal3DLayout/C/N=text0/A=LT/X=0.00/Y=0.00 "\Z32" + time2str2(ttot)
			//doupdate
			//savepict /p=_PictGallery_ /E=-5 /N=Spinoidal3DLayout /o as "Frame3D"
			dowindow /F cylimage
			addmovieframe // /pict=Frame3D
		endif
	endfor
	setdatafolder ::
	make /o /n=(thickness,s.num,s.num) s.density1=mat[q][r]
	setscale /i x, -thickness/2, thickness/2, s.density1
	setscale /i y, -s.num/2, s.num/2, s.density1
	setscale /i z, -s.num/2, s.num/2, s.density1 // this returns the density matrix of material 1 (the matrix) for alignment etc later on
	s.density1 +=gnoise(noise)
	s.density1 = abs(s.density1)
	imagefilter /n=(interpenetration)/o gauss3d s.density1
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
function model3D_Spheres2ln(s)
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
	
	struct ThreeDSystem &s
	if(itemsinlist(s.paramstring,",")<5)
		return -1
	endif
	newdatafolder /o/s CreatingSpheres
	variable minsep = 		str2num(stringfromlist( 0 ,s.paramstring,","))
	variable pd = 				str2num(stringfromlist( 1 ,s.paramstring,","))
	variable thickness = 		s.thickness
	variable noise = 			str2num(stringfromlist( 2 ,s.paramstring,","))
	variable rad2 = 			str2num(stringfromlist( 3 ,s.paramstring,","))
	variable pd2 = 			str2num(stringfromlist( 4 ,s.paramstring,","))
	variable interpenetration = 	str2num(stringfromlist( 5 ,s.paramstring,","))
	

	
	make /o /n=(thickness,s.num,s.num) mat=1,xwave, ywave, zwave
	
	make/B/U /o /n=(thickness,s.num,s.num,30) exmat= (p <= t) || (q <= t) || (r <= t) || (p >= thickness-t) || (q >= s.num-t) || (r >= s.num-t) ? 0 : 1
	make/B/U /o /n=(thickness,s.num,s.num) tempwave
	if(s.movie)
		Execute("Spheres3Ddisp(" +num2str(s.num)+", \""+getwavesdatafolder(mat,2)+"\")")
		Execute("exportgizmo wave as \"testimage\";Spinoidal3DLayout();Spinoidal3DImage(\""+getdatafolder(1)+"testimage\")")
	endif
	setscale /i x, -thickness/2, thickness/2, mat, exmat,xwave, ywave, zwave
	setscale /i y, -s.num/2, s.num/2, mat, exmat,xwave, ywave, zwave
	setscale /i z, -s.num/2, s.num/2, mat, exmat,xwave, ywave, zwave
	xwave = x
	ywave = y
	zwave = z
	redimension /n=(thickness*s.num*s.num) xwave, ywave, zwave
	variable i,radius, orad, cx, cy, cz, failed, fnum =0, xmn,xmx,ymn,ymx,zmn,zmx, loc
	fnum=0

	wave radiuses = getradiusdistribution2ln(s.size, rad2,s.num,thickness,pd, pd2)
	
	sort /R radiuses, radiuses
	for(i=0;i<numpnts(radiuses);i+=1)
		radius = radiuses[i]
		radius =radius < 1 ? 1 : radius
		if(minsep<1)
			orad = radius*(1+minsep/2)
		else
			orad = radius + minsep/2
		endif
		redimension /n=(thickness,s.num,s.num) tempwave
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
		redimension /n=(thickness*s.num*s.num) tempwave
		integrate tempwave /D=intwave
		loc = binarysearch(intwave, enoise(wavemax(intwave)/2)+wavemax(intwave)/2)
		cx = xwave[loc]
		cy = ywave[loc]
		cz = zwave[loc]
		
		// subtract out this sphere from the matrix  // matrix starts at 1s, within this sphere, multiply this by 0, outside multiply by 1
		multithread mat*= (x-cx)^2 + (y-cy)^2 + (z-cz)^2 < radius^2 ? 0 : 1 
		multithread exmat*= (x-cx)^2 + (y-cy)^2 + (z-cz)^2 <= (orad+t)^2 ? 0 : 1 
		if(s.movie)
			execute("ModifyGizmo /n=Spheres3D update=2")
			doupdate
			Execute "exportgizmo wave as \"testimage\""
			//TextBox/w=Spinoidal3DLayout/C/N=text0/A=LT/X=0.00/Y=0.00 "\Z32" + time2str2(ttot)
			doupdate
			savepict /p=_PictGallery_ /E=-5 /N=Spinoidal3DLayout /o as "Frame3D"
			addmovieframe /pict=Frame3D
		endif
	endfor
	setdatafolder ::
	imagefilter /n=(interpenetration)/o gauss3d mat
	duplicate /o mat,s.density1 // this returns the density matrix of material 1 (the matrix) for alignment etc later on
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

function make3darrows(size,ox,oy,oz,amat,dmat)
	wave amat,dmat
	variable size, ox, oy, oz
	
	make /n=(size^3,3)/o arrows3dloc, arrows3dsize
	make /n=(size^3,4)/o arrows3drot
	arrows3dloc[][0]=mod(p,size) +ox
	arrows3dloc[][1]=mod(floor(p/size),size) + oy
	arrows3dloc[][2]=floor(p/(size^2)) + oz
	
	arrows3dsize[][] = sqrt(amat[arrows3dloc[p][0]][arrows3dloc[p][1]][arrows3dloc[p][2]][0]^2 + amat[arrows3dloc[p][0]][arrows3dloc[p][1]][arrows3dloc[p][2]][1]^2 + amat[arrows3dloc[p][0]][arrows3dloc[p][1]][arrows3dloc[p][2]][2]^2) /5
	arrows3drot[][] = cang(0,0,1,amat[arrows3dloc[p][0]][arrows3dloc[p][1]][arrows3dloc[p][2]][0],amat[arrows3dloc[p][0]][arrows3dloc[p][1]][arrows3dloc[p][2]][1],amat[arrows3dloc[p][0]][arrows3dloc[p][1]][arrows3dloc[p][2]][2],q)
	
	make/o/n=(size,size,size) arrowsurface
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
	make/o /n=(360,.5/.006) azplot
	setscale /p x,1,1,azplot
	setscale /p y, .006,.006, azplot
	variable i
	for( i = 1 ; i <= dimsize(azplot,1) ; i+=1)
		azimuthplot(en,i*.006,i*.006+3*.006)
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
	Display/k=1 /n=azgraphwindow /W=(500,100,800,400) azplot[*][5],azplot[*][10],azplot[*][15],azplot[*][20]
	AppendToGraph azplot[*][25],azplot[*][30],azplot[*][35],azplot[*][40],azplot[*][45]
	ModifyGraph margin(right)=10, margin(left)=35, margin(top)=10,gfSize=10
	ModifyGraph mode=7
	ModifyGraph rgb(azplot)=(65535,65535,0),rgb(azplot#1)=(65535,53772,0),rgb(azplot#2)=(65535,31927,0)
	ModifyGraph rgb(azplot#3)=(65535,11762,0),rgb(azplot#5)=(49485,0,0),rgb(azplot#6)=(33436,0,0)
	ModifyGraph rgb(azplot#7)=(16049,0,0),rgb(azplot#8)=(0,0,0)
	ModifyGraph hbFill=5
	ModifyGraph lSize=2
	ModifyGraph log(left)=1
	ModifyGraph tick=2
	ModifyGraph mirror=1
	ModifyGraph standoff=0
	ModifyGraph axOffset(left)=-2.22222,axOffset(bottom)=-0.666667
	ModifyGraph axisOnTop=1
	Label left "Scattering Intensity"
	Label bottom "Azimuthal Angle [deg]"
	ModifyGraph gfSize=16
	ModifyGraph margin(bottom)=50,margin(top)=9,margin(right)=9, margin(left) = 57
	ColorScale/C/N=text0/A=RC/X=0.00/Y=0.00/E  ctab={0.03,0.27,YellowHot,1},axisRange={0.27,0.03}, height=150
	ColorScale/C/N=text0 width=8, tickLen=2, lblMargin=0
	AppendText "Momentum Transfer [nm\\S-1\\M]"
	SetAxis left 0.001,100
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
	multithread arrowsyay[][0] = sqrt(zcomp[p]^2 + ycomp[p]^2) * 10
	//arrowsyay[][0] = arrowsyay[p][0]< 3 &&  sqrt(zcomp[p]^2 + ycomp[p]^2+ xcomp[p]^2)>.5? 3 : arrowsyay[p][0]
	yloc = arrowsyay[p][0] <4 ? nan : yloc[p]
	dowindow /k alignmentmap
	Display /k=1/n=alignmentmap/W=(270.6,59.6,518.4,308.6)/K=1  yloc vs xloc as "Slice of Alignment"
	AppendImage/w=alignmentmap slice
	if(rev)
		ModifyImage/w=alignmentmap slice ctab= {1,0,yellowhot,0}
	else
		ModifyImage/w=alignmentmap slice ctab= {0,1,yellowhot,0}
	endif
	ModifyGraph/w=alignmentmap mode=3,gfSize=12
	ModifyGraph/w=alignmentmap rgb=(65535,65535,65535)
	ModifyGraph/w=alignmentmap msize=0.5,marker=42,tlOffset=-5
	ModifyGraph/w=alignmentmap mrkThick=0.1,margin(left)=26,margin(bottom)=26
	ModifyGraph/w=alignmentmap arrowMarker(yloc)={arrowsyay,.75,0,0,1}
	ModifyGraph/w=alignmentmap mirror=0,margin(top)=1,margin(right)=1
	ModifyGraph/w=alignmentmap axOffset(left)=-4.46154,axOffset(bottom)=-0.962963
	setscale /p x, 0,5,slice
	setscale /p y, 0,5,slice
	setaxis /w=alignmentmap left 0,min(225,5*dimsize(slice,0))
	setaxis /w=alignmentmap bottom 0,min(225,5*dimsize(slice,1))
	ModifyGraph /w=alignmentmap tick=2,mirror=1,axisOnTop=1,standoff=0,tkLblRot(left)=90
	Label /w=alignmentmap left "Size Scale [nm]";DelayUpdate
	Label /w=alignmentmap bottom "Size Scale [nm]"
	if(rev)
		make /n=3 /o tics = {1,.5,0}
	else
		make /n=3 /o tics = {0,.5,1}
	endif
	make /n=3 /t /o ticlabels = {"P3HT", "Mixture", "PCBM"}
	ColorScale /w=alignmentmap/C/N=text0/A=RB image=slice;DelayUpdate
	ColorScale /w=alignmentmap/C/N=text0 userTicks={tics,ticlabels}
	ColorScale /w=alignmentmap/C/N=text0 vert=0,side=2,tickLblRot=0
	ColorScale /w=alignmentmap/C/N=text0/X=3.00/Y=3.00 width=100,height=15
	ColorScale /w=alignmentmap/C/N=text0 lblMargin=5,tickLen=3.00
	setdatafolder foldersave
end
Window alignmentmap() : Graph
	PauseUpdate; Silent 1		// building window...
	String fldrSav0= GetDataFolder(1)
	SetDataFolder root:alignmap:
	Display /W=(270.6,59.6,518.4,308.6)/K=1  yloc vs xloc as "Slice of Alignment"
	AppendImage slice
	ModifyImage slice ctab= {0,1,yellowhot,0}
	SetDataFolder fldrSav0
	ModifyGraph margin(left)=32,margin(bottom)=33,margin(top)=1,margin(right)=1,gfSize=12
	ModifyGraph height={Plan,1,left,bottom}
	ModifyGraph mode=3
	ModifyGraph marker=42
	ModifyGraph rgb=(65535,65535,65535)
	ModifyGraph msize=0.5
	ModifyGraph mrkThick=0.1
	ModifyGraph arrowMarker(yloc)={:alignmap:arrowsyay,0.75,0,0,1}
	ModifyGraph tick=2
	ModifyGraph mirror=1
	ModifyGraph standoff=0
	ModifyGraph axOffset(left)=-4.46154,axOffset(bottom)=-0.962963
	ModifyGraph tkLblRot(left)=90
	ModifyGraph axisOnTop=1
	Label left "Size Scale [nm]"
	Label bottom "Size Scale [nm]"
	SetAxis left 0,225
	SetAxis bottom 0,225
	ColorScale/C/N=text0/A=RB/X=3.00/Y=3.00 image=slice, vert=0, side=2, tickLen=3
	ColorScale/C/N=text0 userTicks={:alignmap:tics,:alignmap:ticlabels}, lblMargin=5
EndMacro
function scatterimage(en)
	variable en
	wave scatter3dsave
	if(!waveexists(scatter3dsave))
		print "no scattering in this directory"
		return 0
	endif
	make /n=(1000,1000) /o scatterdisp
	setscale /i x,-.12,.12, scatterdisp
	setscale /i y,-.12,.12, scatterdisp
	setscale /i z,260,320,scatter3dsave
	setscale /i x,-.5,.5,scatter3dsave
	setscale /i y,-.5,.5,scatter3dsave
	scatterdisp = sqrt(x^2 +y^2) < .005 ? 0 : scatter3dsave(x)(y)(en)
	matrixfilter /n=20 gauss scatterdisp
	
	dowindow/k scatterview
	Display /W=(100,100,425,425) /k=1 /n=scatterview
	appendimage scatterdisp
	ModifyImage scatterdisp ctab= {0,500,Terrain,0}
	ColorScale/w=scatterview/C/N=text0/A=RC/X=0.00/Y=0.00 image=scatterdisp;DelayUpdate
	ColorScale/w=scatterview/C/N=text0 "Scattering Intensity"
	TextBox/C/N=text1/F=0/A=LT/X=10.00/Y=10.00/B=1/G=(65535,65535,65535) "\\F'Arial Black'\\Z20"+num2str(en)+" eV"
	ModifyGraph gfSize=16, margin(left)=78,margin(bottom)=50,margin(top)=11,margin(right)=11
	Label bottom "Momentum Transfer Qy [nm\\S-1\\M]"
	Label left "Momentum Transfer Qz [nm\\S-1\\M]"
	ModifyGraph height={Plan,1,left,bottom}
	ModifyGraph tick=2,mirror=1,minor=1,standoff=0,axRGB=(65535,65535,65535)
	ModifyImage scatterdisp ctab= {0,150,YellowHot,0}
	ColorScale/C/N=text0 axisRange={,155}
	azploter(en)
end

Window RatioGraph() : Graph
	PauseUpdate; Silent 1		// building window...
	Display /W=(400,10,300,300)
	AppendImage/T :packages:ScatterSim3D:test:ratio3dtemp
	ModifyImage ratio3dtemp ctab= {-1,1,RedWhiteGreen,0}
	ModifyGraph margin(left)=43,margin(bottom)=9,margin(top)=28,margin(right)=9,gfSize=14
	ModifyGraph log(top)=1
	ModifyGraph tick=2
	ModifyGraph mirror=2
	ModifyGraph nticks(left)=6,nticks(top)=8
	ModifyGraph minor=1
	ModifyGraph sep(left)=2
	ModifyGraph standoff=0
	ModifyGraph btLen=5
	ModifyGraph ftLen(left)=5
	ModifyGraph tlOffset=-2
	Label left "X-ray Energy [eV]"
	SetAxis/R left 275,295
	SetAxis top 0.01,1
EndMacro

macro doradiograph()
	duplicate/o ratio3dvsen, ratio3dtemp
	make/o /n=(dimsize(ratio3dtemp,0)) /o normalizewave
	normalizewave = ratio3dtemp[p](260)
	ratio3dtemp -= normalizewave[p]
	RatioGraph()
	ModifyImage ratio3dtemp ctab= {-1,1,RedWhiteBlue,0}
	smooth 5, ratio3dtemp
	ModifyGraph log(top)=1
	ModifyGraph mirror=2
	ModifyGraph nticks(left)=12
	ModifyGraph minor=1
	ModifyGraph fSize=8
	ModifyGraph standoff=0
	ModifyGraph tkLblRot(left)=90
	ModifyGraph btLen=3
	ModifyGraph tlOffset=-2
	SetAxis/A/R left
	SetAxis/R left 290,280
	SetAxis/R left 295,275
	ModifyImage ratio3dtemp ctab= {-1,1,RedWhiteGreen,0}
	SetAxis/R left 275,295
	SetAxis top 0.01,1
	ModifyGraph gfSize=16
	ModifyGraph fSize=0
	ModifyGraph nticks(left)=5
	ModifyGraph gfSize=14
	ModifyGraph tkLblRot=0
	Label left "X-ray Energy [eV]"
	ModifyGraph nticks(left)=6
	ModifyGraph tick=2,btLen=5
	ModifyGraph ftLen(left)=5
	ModifyGraph sep(left)=2
	ShowTools/A arrow
	HideTools/A
	DoWindow/C/R RatioGraph
end