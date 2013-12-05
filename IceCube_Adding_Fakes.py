#2013_11_05 Creating File Structure and Galaxy Bias

#Before you run this make sure you have the input images (both science and mask) in this directory along with the following files:
#default.conv
#default.nnw
#Pipe_sexfile.sex
#PTF_Transform_Param.param

#!!!!!!!!!!!!!!!!IMPORTANT!!!!!!!!!!!!!--- YOU NEED THE MODULES BELOW

import numpy, os, random, pylab, glob, shutil, time, subprocess
from multiprocessing import Pool
from astropy.io import fits


def file_structure():
	if not os.path.exists('Weight_Map'):
		os.makedirs('Weight_Map')
	if not os.path.exists('Catalog'):
		os.makedirs('Catalog')
	if not os.path.exists('Fake_Star_Catalog'):
		os.makedirs('Fake_Star_Catalog')
	if not os.path.exists('Fakes'):
		os.makedirs('Fakes')	
	if not os.path.exists('Galaxies'):
		os.makedirs('Galaxies')
	if not os.path.exists('Bad_Images'):
		os.makedirs('Bad_Images')



#---THIS WILL CREATE THE .fits WEIGHT_IMAGE
def weight_map(maskfile, science_image): #Taken from Create_Weight_Map.py
	'''
	weight_map(maskfile, science_image) It takes the PTF mask file and creates a Weight map
	'''
	try:
		os.remove('Weight_Map/'+science_image+'_Weight_Map.fits') #removes a file if it has the same name because astropy doesn't like overwriting existing files
		#print 'Old ', science_image, '_Weight_map.fits was removed, new file will be written'
	except OSError:
		pass
		
	hdulist_mask = fits.open(maskfile+'.fits',ignore_missing_end=True) #opens maskfile
	
	hdulist_mask[0].data=((hdulist_mask[0].data) < 3).astype(float)

	#print mask_data

	hdulist_mask.writeto('Weight_Map/'+science_image+'_Weight_Map.fits') #New fits file. Its the WEIGHT_IMAGE
	
#---THIS STEP CREATE THE .sex FILE FOR THE SExtractor routine
#---!!!!!!!!!!!!---!!!!!----!!!!!!----REDUNDANT CODE-----!!!!!-Sex_file_maker is not needed anymore---The parameters that vary on a per image basis are added into Sextract------ 
def Sex_file_maker(science_image,zeropoint,seeing,saturation,gain): #Taken from Create_sexed_catalog.py

	z=open('Pipe_sexfile.sex') #Since this is a stock file perhaps it doesn't need to be an argument in def..? #RESOLVED 29-10-13 15:05
	#v=open('temp.sex','w') #The temporary output .sex file
	j=open('Sex_File/'+science_image+'_Custom_Sex.sex','w')
	j.write(str('#Edited ')+str(time.strftime('%Y-%m-%d %H:%M:%S'))+'\n')
	#sex_file_array=[]
	for line in z:
		if line.startswith('WEIGHT_IMAGE'):
			j.write(str('WEIGHT_IMAGE')+'	'+'Weight_Map/'+str(science_image)+'_Weight_map.fits'+'\n')	
			continue		
		if line.startswith('MAG_ZEROPOINT'):			
			j.write(str('MAG_ZEROPOINT')+'	'+str(zeropoint)+'\n')
			continue	
		if line.startswith('SEEING_FWHM'):
			j.write(str('SEEING_FWHM')+'	'+str(seeing)+'\n')
			continue
		if line.startswith('SATURVAL'):			
			j.write(str('SATURVAL')+'	'+str(saturation)+'\n')
			continue
		if line.startswith('GAIN'):
			j.write(str('GAIN')+'	'+str(gain)+'\n')
			continue
		if line.startswith('CATALOG_NAME'):
			j.write(str('CATALOG_NAME')+'	'+'Catalog/'+str(science_image)+'_Catalog.cat'+'\n')
			continue
		else:
			j.write(line)
	#for g in range(len(sex_file_array)):
		#j.write(sex_file_array[g])
		#print sex_file_array[g]
	j.close() 

#---THIS STEP WILL CREATE THE .cat FILE
def Sextract(science_image,zeropoint,seeing,saturation,gain): #Taken from Create_sexed_catalog.py
	#'sex -c '+science_image+'_Custom_Sex.sex '+science_image+'.fits -PARAMETERS_NAME PTF_Transform_Param.param -FILTER_NAME default.conv'
	subprocess.call('sex -c Pipe_sexfile_CMD.sex '+science_image+'.fits -PARAMETERS_NAME PTF_Transform_Param.param -FILTER_NAME default.conv -CATALOG_NAME Catalog/'+str(science_image)+'_Catalog.cat -WEIGHT_IMAGE Weight_Map/'+str(science_image)+'_Weight_map.fits -MAG_ZEROPOINT'+'	'+str(zeropoint)+' -SEEING_FWHM '+str(seeing)+' -SATUR_LEVEL '+str(saturation)+' -GAIN '+str(gain)+' -VERBOSE_TYPE QUIET',shell=True)

def Enough_Objects(science_image):
	enough=True
	test=os.popen('wc -l Catalog/'+science_image+'_Catalog.cat').read()
	rows=test.split()
	if float(rows[0])<217:
			return False

def Detections_Hist(science_image):
	histo=[]
	this=os.getcwd()
	for i in range(1,len(os.listdir(this+'/Catalog'))):
		filename=os.listdir(this+'/Catalog')[i]
		name=os.popen('wc -l '+'Catalog/'+filename).read()
		ln=name.split()
		num=int(ln[0])	
		histo.append(num)
	pylab.hist(histo, bins=10)
	pylab.xlabel('Sources Detected in Each Image')
	pylab.ylabel('Number of Images')
	pylab.show()

	
def Selecting_Bright(science_image):


	f=open('Catalog/'+science_image+'_Catalog.cat') #CHANGE THIS SO THAT IT CAN TAKE ANY CATALOG INPUT
	fin=f.readline()
	
	xcord=[]
	ycord=[]
	flux_array=[]
	mag_array=[]
	background_array=[]
	while fin: #DEFINING BINS BASED ON MAGNITUDE OF STAR
		
		if fin.startswith('#'):
			#h.write(fin)
			#print fin
			fin=f.readline()
			continue
		ln=fin.split()
		#print ln[15]
		mag=float(ln[2]) #magnitude of star	
		x=float(ln[3])
		y=float(ln[4])
		flux=float(ln[1])	
		background=float(ln[5])

		if mag<11.0 and mag>10.0:
				if float(ln[7])<0.3: #Not Elliptical
					if x>100.0 and x<1948.0 and y>100.0 and y<3996.0: #No Edge Stars
						#h.write(fin)
						#g.write(str(x)+' '+str(y)+'\n')
						#print fin
						xcord.append(x); ycord.append(y); mag_array.append(mag); flux_array.append(flux); background_array.append(background)
						#print fin
		fin=f.readline()
	f.close()
	return xcord, ycord, mag_array, flux_array, background_array

def selecting_galaxies(science_image,):
	f=open('Catalog/'+science_image+'_Catalog.cat')
	g=open('Galaxies/'+science_image+'_Galaxy_Catalog.cat','w')
	l=open('Galaxies/'+science_image+'_Galaxy_regions.reg','w')
	m=open('Galaxies/'+science_image+'_Star_regions.reg','w')
	fin=f.readline()
	X2_array=[]
	Y2_array=[]
	FWHM_array=[]
	Ellipticity_array=[]
	class_s_array=[]
	counts=0
	while fin:
		if fin.startswith('#'):
			#h.write(fin)
			#print fin
			fin=f.readline()
			continue
		ln=fin.split()
		class_s=float(ln[9])
		gal_mag=float(ln[2])
		xcord=float(ln[3])
		ycord=float(ln[4])
		X2=float(ln[13])
		Y2=float(ln[14])
		X2_array.append(X2)
		Y2_array.append(Y2)
		
		FWHM=float(ln[16])
		FWHM_array.append(FWHM)
		elips=float(ln[7])
		Ellipticity_array.append(elips)
		class_s_array.append(class_s)
		if class_s<0.5 and gal_mag>10:
			if FWHM<15:
				if xcord>50.0 and xcord<1998.0 and ycord>50.0 and ycord<4046.0: #No Edge Galaxies
					g.write(fin)
					#g.write(str((ln[0]))+' '+str((ln[1]))+' '+str((ln[2]))+' '+str((ln[3]))+' '+str((ln[4]))+' '+str(ln[5])+' '+str((ln[6]))+' '+str((ln[7]))+' '+str((ln[8]))+' '+str((ln[9]))+' '+str((ln[10]))+' '+str((ln[11]))+' '+str(ln[12])+' '+str((ln[13]))+' '+str((ln[14]))+' '+str((ln[15]))+' '+str((ln[16]))+'\n')
					l.write(str(xcord)+' '+str(ycord)+'\n')
					counts+=1
		else:
			m.write(str(xcord)+' '+str(ycord)+'\n')		

				
		fin=f.readline()

	f.close()
	g.close()
	l.close()
	m.close()
	
	#galaxy_image(science_image, counts)


def galaxy_image(science_image, count):
	k=open('Galaxies/'+science_image+'_Galaxy_Catalog.cat')
	
	size=40.0 #How big in pixels (size*size) the square around each galaxy will be
	size2=size/2
	square=int(count**0.5)+1
	res=square*size
	

	hdulist_sci= fits.open(science_image+'.fits',ignore_missing_end=True)
	
	science_data= hdulist_sci[0].data

	try:
		os.remove('Galaxies/Primer_Galaxy_Grid.fits') #removes a file if it has the same name because astropy doesn't like overwriting existing files
		#print 'Old Galaxies Primer was removed, new file will be written'
	except OSError:
		pass

	try:
		os.remove('Galaxies/'+science_image+'Galaxy_Grid.fits') #removes a file if it has the same name because astropy doesn't like overwriting existing files
		#print 'Old '+science_image+'Galaxy_Grid.fits, new file will be written'
	except OSError:
		pass

	new=numpy.zeros(shape=(res,res))
	hdu=fits.PrimaryHDU(new)
	hdulist=fits.HDUList([hdu])
	hdulist.writeto('Galaxies/Primer_Galaxy_Grid.fits')

	hdulist_new_sci=fits.open('Galaxies/Primer_Galaxy_Grid.fits',ignore_missing_end=True)
	new_sci=hdulist_new_sci[0].data
	#print new_sci
	kin=k.readline()
	
	
	
	#while fin:
	
	for y in range(0,int(res),int(size)):
		for x in range(0,int(res),int(size)):
			if not kin:
				break
			#print '----'
			#print kin
			ln=kin.split()
			xco=float(ln[3])
			yco=float(ln[4])

			startx=int(xco-size2)
			starty=int(yco-size2)		
			finx=int(xco+size2)
			finy=int(yco+size2)

			#print float(ln[0])
			#print xco,yco, startx, finx, 'diff= ', startx-finx
			#print xco,yco, starty, finy, 'diff= ', starty-finy
			
			new_sci[y:(y+size),x:(x+size)]=science_data[starty:finy,startx:finx]
			
			
			kin=k.readline()
			

	hdulist_new_sci.writeto('Galaxies/'+science_image+'Galaxy_Grid.fits')
	k.close()
	#print square

def Scaling(science_image ,xcord, ycord, mag_array, flux_array, background_array, zpt, fake_stars, CCD_Num):
	f=open('Fake_Star_Catalog/'+science_image+'_Fake_Star_Catalog.dat','w') #--This will need to be file specific at a later date #RESOLVED 29-10-13 15:31
	#print 'lengths', len(mag_array), len(flux_array)
	ranmagarray=[]
	for i in range(0,fake_stars):
		ran_mag=random.uniform(14.0, 22.0) #The fake stars will be in this range of magnitudes
		ran_flux=10.0**((ran_mag-zpt)/(-2.5))
		ranmagarray.append(ran_mag)
		star=int(random.uniform(0,len(xcord)-1))
		
		scaling_factor=((ran_flux)/flux_array[star])
		
		newX=random.uniform(100.0,1948.0)#No edge stars, this creates a 100 pixel border
		newY=random.uniform(100.0, 3996.0)#no edge stars
		#print background_array[star]
		#'Old X', 'Old Y', 'New X', 'New Y', 'Old Mag', 'Old Flux', 'New Mag', 'New Flux', 'Scaling Factor' 		
		f.write(str(xcord[star])+' '+str(ycord[star])+' '+str(newX)+' '+str(newY)+' '+str(mag_array[star])+' '+str(flux_array[star])+' '+str(ran_mag)+' '+str(ran_flux)+' '+str(background_array[star])+' '+str(scaling_factor)+' '+str(CCD_Num)+'\n')
		
		i+=1
	#pylab.hist(ranmagarray, bins=8)
	#pylab.show()
	f.close()
#----------GALAXY BIASING----------------------------------NEW----------!!!!!!!!!!!!!!!!!!!



def add_fakes_2galaxy(science_image,boxsize):
	#This step finds the galaxies and adds fake stars to them
	h=open('Galaxies/'+science_image+'_Galaxy_Catalog.cat') #Opens the Galaxy catalog
	f=open('Fake_Star_Catalog/'+science_image+'_Fake_Star_Catalog.dat') #Opens the fake star catalog
	reg=open('Fake_Star_Catalog/'+science_image+'_Fakes_Star_Regions.reg','w') #creates region file
	hin=h.readline() #reads first line
	fin=f.readline() #reads first line

	fake_star_array=[] #prepares and array for fake stars
	gal_line_array=[] #prepares and array for galaxies
	while hin: #Adds all galxies into an array
		gal_line_array.append(hin)
		hin=h.readline()
	while fin: #adds all stars into an array
		fake_star_array.append(fin)
		fin=f.readline()
	h.close()
	f.close()

	hdulist_sci= fits.open(science_image+'.fits',ignore_missing_end=True) #THE SCIENCE DATA THAT WILL OPENED AND MODIFIED
	science_data= hdulist_sci[0].data

	#print len(fake_star_array), ' Fake Stars have been added to ', len(fake_star_array), ' Galaxies'
	#print gal_line_array[1]
	#print len(gal_line_array)
	
	j=open('Fakes/'+science_image+'_Flux_Boxes.dat','w')
	for i in range(0,len(fake_star_array)): #Will only add 200 fake stars to 200 Galaxies
		host_galaxy=gal_line_array.pop(random.randrange(0,len(gal_line_array))) #selecting a random host galaxy. Used .pop() so that the same galaxy isnt chosen twice
		source_star=fake_star_array.pop(random.randrange(0,len(fake_star_array))) #selecting a random source star. Used .pop() so that the same star isnt chosen twice
		
		#print y 
		#print 'len: ',len(gal_line_array)
		ln=host_galaxy.split()
		x=float(ln[3])
		y=float(ln[4])
		CXX=float(ln[10])
		CYY=float(ln[11])
		CXY=float(ln[12])
		R=3.0
		#----------------THIS IS A GRID USED FOR CALCULATING STAR POSITIONS !!NOT!! FOR SCALING STARS
		#Draw a large grid around the galaxy of say 20,20 pixels. Run through that grid for every x and y and if it satisfies the equation on page 32 of the sextractor manual then append it to an
		#array. Then randomly choose a coordinate and insert a fake star there.
		grid_size=20
		grid_startx=int(x-grid_size/2.0)
		grid_starty=int(y-grid_size/2.0)		
		grid_finx=int(x+grid_size/2.0)
		grid_finy=int(y+grid_size/2.0)
		#----------------
		#good_x=[]
		#good_y=[]
		good_coords=[]
		#This loops runs through the grid and finds the pixel co-ordinates that satisfy the inequality below
		for q in range(grid_starty,grid_finy):
			for w in range(grid_startx,grid_finx):
				ell_p=(CXX*((w-x)*(w-x)))+(CYY*((q-y)*(q-y)))+(CXY*((w-x)*(q-y))) #IF THIS INEQUALITY IS SATISFIED THEN A STAR CAN GO AT THIS LOCATION 
				if ell_p<=3.0:
					good_coords.append([w,q])
					#good_y.append(q)
					#good_x.append(w)
		#print i, good_coords
		fake_star_positon=random.choice(good_coords) #Chooes a random co-ordinate
		kn=source_star.split()
		sourcex=float(kn[0]) #stars current x location
		sourcey=float(kn[1]) #stars current y location
		newx=fake_star_positon[0] #where the star will go x
		newy=fake_star_positon[1] #where the star will go y

		#Saving Flux info Pre Fakes Addition

		pixel_fluxes=[]
		
		#print science_image
		#print 'newx newy sci:	', newx, newy, (science_data[newy,newx])
		pixel_fluxes.append((science_data[newy,newx]))
		for s in boxsize:
			start_box_x=newx-int(s/2)
			fin_box_x=newx+int(s/2)+1
			start_box_y=newy-int(s/2)
			fin_box_y=newy+int(s/2)+1
			#print 'Box X Dimensions: ', fin_box_x - start_box_x
			#print 'Box Y Dimensions: ', fin_box_y - start_box_y
			Box_Matrix=[[0 for x in xrange(s)] for x in xrange(s)]
			#print Box_Matrix
			countg=0
			
			for g in range (start_box_y,fin_box_y):
				counth=0
				for h in range(start_box_x,fin_box_x):
					Box_Matrix[countg][counth]=(science_data[g,h])
					#print Box_Matrix[countg][counth]
					counth+=1
					#print Box_Matrix
				countg+=1
			pixel_fluxes.append(Box_Matrix)
					
			#print 's newx newy sci:	', s, newx, newy, (science_data[newy,newx])
		j.write(str(pixel_fluxes)+'\n')
		
		reg.write(str(newx)+' '+str(newy)+'\n') #fake star region file
		
		scale_fac=float(kn[9]) #scale factor
		back=float(kn[8]) #background
		
		#---Old area to be scaled---
		startx=int(sourcex-4.0)
		starty=int(sourcey-4.0)		
		finx=int(sourcex+4.0)
		finy=int(sourcey+4.0)
		
		#---New area to have flux added---
		Nstartx=newx-4.0
		Nstarty=newy-4.0
		Nfinx=newx+4.0
		Nfiny=newy+4.0


		newdata=numpy.ones((8,8)) #Preparing a blank gird for scaled objects

		newdata[0:8,0:8]=(((science_data[starty:finy,startx:finx]))-back)*scale_fac #inserting scaled object

		science_data[Nstarty:Nfiny, Nstartx:Nfinx]= (science_data[Nstarty:Nfiny, Nstartx:Nfinx]) + newdata #Modifying the science image
		
		



	try:
		os.remove('Fakes/'+science_image+'_WITH_FAKES.fits') #removes a file if it has the same name because astropy doesn't like overwriting existing files
		#print 'Old Fake Removed'
	except OSError:
		pass
		
	hdulist_sci.writeto('Fakes/'+science_image+'_WITH_FAKES.fits') #Saving image after loop of 200 Stars is complete
	j.close()
	reg.close()	
	print 'Fake Stars Added to Galaxies in the Image: ', science_image





def Execute(run):
	#print '!!!!!!', run	
	prefix=run[0]
	suffix=run[1]
	science_image=prefix+'scie'+suffix
	#print '!!!!!!', science_image

	maskfile=prefix+'mask'+suffix
	#print  '######', maskfile

	
	#print 'Name: ', science_image

	try:
		hdulist_multi_sci=fits.open(science_image+'.fits')
		#print '++++ multi_mask assign ', science_image
	except IOError or Warning or UnboundLocalError:
		hdulist_multi_sci.close()
		
		shutil.move(science_image+'.fits','Bad_Images/')
		shutil.move(maskfile+'.fits','Bad_Images/')
		return
	try:
		multi_mask=fits.open(maskfile+'.fits')
		multi_mask.close()
	except IOError or Warning or UnboundLocalError:
		shutil.move(science_image+'.fits','Bad_Images/')
		shutil.move(maskfile+'.fits','Bad_Images/')
		return
		
	zeropoint=float(hdulist_multi_sci[0].header['IMAGEZPT'])
	seeing=float(hdulist_multi_sci[0].header['SEEING'])
	saturation=float(hdulist_multi_sci[0].header['SATURVAL'])
	gain=float(hdulist_multi_sci[0].header['GAIN'])
	CCD_Num=float(hdulist_multi_sci[0].header['CCDID'])


	fake_stars= 100 #number of fake stars per image (integer please!)

	hdulist_multi_sci.close()
	weight_map(maskfile, science_image)
	#print '@@@@@ WEIGHTMAP CREATED ', science_image 
	#Adding Weight_Map

	Sextract(science_image,zeropoint,seeing,saturation,gain)
	#print '@@@@@ sextract done'
	catsize=Enough_Objects(science_image)
	#print '@@@@@Enough_Objects Done'
	if catsize==False:
			#print science_image, 'didn\'t have enough objects detected so it was moved to Bad_Images/ and the newly created weight map, sex file and catalog have been deleted'
			shutil.move(science_image+'.fits','Bad_Images/')
			shutil.move(maskfile+'.fits','Bad_Images/')
			os.remove('Weight_Map/'+science_image+'_Weight_Map.fits')
			os.remove('Catalog/'+science_image+'_Catalog.cat')
			return
	#else: print '@@@ CATSIZE IS TRUE'
		
	#print 'Sextract done'
	x, y, mag, flux, back = Selecting_Bright(science_image)
	#print '@@@@ SELECTING BRIGHT DONE'
	#print len(x)
	if len(x)<2:
			#print 'Not enough Objects met the source star criteria in ', science_image
			shutil.move(science_image+'.fits','Bad_Images/')
			shutil.move(maskfile+'.fits','Bad_Images/')
			os.remove('Weight_Map/'+science_image+'_Weight_Map.fits')
			os.remove('Catalog/'+science_image+'_Catalog.cat')
			return
	#print '@@@@@LENGTH IS TRUE'
	Scaling(science_image, x, y, mag, flux, back, zeropoint, fake_stars, CCD_Num)
	#print '@@@@Scaling done'
	selecting_galaxies(science_image)
	boxsize=[3,5,7]
	add_fakes_2galaxy(science_image,boxsize)

#-----------------------------------RUN PIPELINE------------------------------------------

file_structure()

data_inputs=[]

#pwd=os.getcwd()
for files in glob.glob('*.fits'):
	full_file=os.path.splitext(files)[0]
	ln=full_file.split('_')
	if ln[4]=='scie': 
		data_inputs.append([(str(ln[0])+'_'+str(ln[1])+'_'+str(ln[2])+'_'+str(ln[3]))+'_','_'+str(ln[5])+'_'+str(ln[6])+'_'+str(ln[7])+'_'+str(ln[8])+'_'+str(ln[9])])
#print data_inputs
#for  run in data_inputs:
processors=4
pool=Pool(processors)

pool.map(Execute,data_inputs)

pool.close()


