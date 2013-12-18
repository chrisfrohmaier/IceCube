#2013_11_05 Creating File Structure and Galaxy Bias

#Before you run this make sure you have the input images (both science and mask) in this directory along with the following files:
#default.conv
#default.nnw
#Pipe_sexfile_CMD.sex
#PTF_Transform_Param.param

#!!!!!!!!!!!!!!!!IMPORTANT!!!!!!!!!!!!!--- YOU NEED THE MODULES BELOW

import numpy, os, random, pylab, glob, shutil, time, subprocess
from multiprocessing import Pool
from astropy.io import fits


def file_structure():
	#if not os.path.exists('Weight_Map'):
		#os.makedirs('Weight_Map')
	if not os.path.exists('Results'):
		os.makedirs('Results')
	if not os.path.exists('Results/Catalog'):
		os.makedirs('Results/Catalog')
	if not os.path.exists('Results/Fake_Star_Catalog'):
		os.makedirs('Results/Fake_Star_Catalog')
	if not os.path.exists('Results/Fakes_added'):
		os.makedirs('Results/Fakes_added')	
	if not os.path.exists('Results/Galaxies'):
		os.makedirs('Results/Galaxies')




#---THIS WILL CREATE THE .fits WEIGHT_IMAGE
def weight_map(maskfile, science_image): #Taken from Create_Weight_Map.py
	'''
	weight_map(maskfile, science_image) It takes the PTF mask file and creates a Weight map
	'''
	try:
		os.remove('Results/Weight_Map/'+science_image[1]+'_Weight_Map.fits') #removes a file if it has the same name because astropy doesn't like overwriting existing files
		#print 'Old ', science_image, '_Weight_map.fits was removed, new file will be written'
	except OSError:
		pass
	W_DEST=str('Results/Weight_Map/'+science_image[1]+'_Weight_Map.fits')	
	print 'W_DEST: ', W_DEST
	hdulist_mask = fits.open(maskfile+'.fits',ignore_missing_end=True) #opens maskfile
	
	hdulist_mask[0].data=((hdulist_mask[0].data) > 0).astype(float)

	#print mask_data

	hdulist_mask.writeto('Results/Weight_Map/'+science_image[1]+'_Weight_Map.fits') #New fits file. Its the WEIGHT_IMAGE
	
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
			j.write(str('WEIGHT_IMAGE')+'	'+'Results/Weight_Map/'+str(science_image)+'_Weight_map.fits'+'\n')	
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
			j.write(str('CATALOG_NAME')+'	'+'Results/Catalog/'+str(science_image)+'_Catalog.cat'+'\n')
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
	subprocess.call('sex -c Pipe_sexfile_CMD.sex '+science_image[0]+science_image[1]+'.fits -PARAMETERS_NAME PTF_Transform_Param.param -FILTER_NAME default.conv -CATALOG_NAME Results/Catalog/'+science_image[1]+'_Catalog.cat -WEIGHT_IMAGE '+science_image[0]+science_image[1]+'.weight.fits -MAG_ZEROPOINT'+'	'+str(zeropoint)+' -SEEING_FWHM '+str(seeing)+' -SATUR_LEVEL '+str(saturation)+' -GAIN '+str(gain)+' -VERBOSE_TYPE QUIET',shell=True)

def Enough_Objects(science_image):
	enough=True
	test=os.popen('wc -l Results/Catalog/'+science_image[1]+'_Catalog.cat').read()
	rows=test.split()
	if float(rows[0])<217:
			return False

def Detections_Hist(science_image):
	histo=[]
	this=os.getcwd()
	for i in range(1,len(os.listdir(this+'/Results/Catalog'))):
		filename=os.listdir(this+'/Results/Catalog')[i]
		name=os.popen('wc -l '+'Results/Catalog/'+filename).read()
		ln=name.split()
		num=int(ln[0])	
		histo.append(num)
	pylab.hist(histo, bins=10)
	pylab.xlabel('Sources Detected in Each Image')
	pylab.ylabel('Number of Images')
	pylab.show() #Do Not USE!!

	
def Selecting_Bright(science_image):


	f=open('Results/Catalog/'+science_image[1]+'_Catalog.cat') #CHANGE THIS SO THAT IT CAN TAKE ANY CATALOG INPUT
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

		if mag<14 and mag>11.0:
				if float(ln[7])<0.3: #Not Elliptical
					if float(ln[9])>0.8: #Considered a good star
						#if x>100.0 and x<1948.0 and y>100.0 and y<3996.0: #No Edge Stars
						if int(ln[8])==0:
							#print 'Good Mag!!!!!!!!!!!!!!!!', mag
							#h.write(fin)
							#g.write(str(x)+' '+str(y)+'\n')
							#print fin
							xcord.append(x); ycord.append(y); mag_array.append(mag); flux_array.append(flux); background_array.append(background)
							#print fin
		fin=f.readline()
	f.close()
	#print 'Mag Arrays', mag_array
	return xcord, ycord, mag_array, flux_array, background_array

def selecting_galaxies(science_image,):
	f=open('Results/Catalog/'+science_image[1]+'_Catalog.cat')
	g=open('Results/Galaxies/'+science_image[1]+'_Galaxy_Catalog.cat','w')
	l=open('Results/Galaxies/'+science_image[1]+'_Galaxy_regions.reg','w')
	m=open('Results/Galaxies/'+science_image[1]+'_Star_regions.reg','w')
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
		'''
	pylab.subplot(2,2,1)
	pylab.xlabel('FWHM (Pixels)')
	pylab.ylabel('Number of Objects')
	pylab.title('FWHM Histogram')
	pylab.hist(FWHM_array,bins=200)
	#pylab.savefig('Galaxy_FWHM_Histogram')
	
	pylab.subplot(2,2,2)
	pylab.ylim([0,6])
	pylab.xlim([0,6])	
	pylab.xlabel('X2')
	pylab.ylabel('Y2')
	pylab.title('2nd Moments')
	pylab.scatter(X2_array,Y2_array,s=8,color='r') #Temp addition

	pylab.subplot(2,2,3)
	pylab.xlabel('ELLIPTICITY')
	pylab.ylabel('Number')
	pylab.hist(Ellipticity_array,bins=20)
	#pylab.savefig('2nd_Order')

	pylab.subplot(2,2,4)
	pylab.xlabel('Class 0=Galaxy 1=Star')
	pylab.ylabel('Number')
	pylab.hist(class_s_array,bins=100)
	
	pylab.suptitle(science_image[1]+'_Characteristics.png')

	pylab.savefig(science_image[1]+'_Characteristics.png',dpi=300)
	#pylab.show()#Temp addition
	'''
	f.close()
	g.close()
	l.close()
	m.close()
	pylab.close()
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
		os.remove('Results/Galaxies/Primer_Galaxy_Grid.fits') #removes a file if it has the same name because astropy doesn't like overwriting existing files
		#print 'Old Galaxies Primer was removed, new file will be written'
	except OSError:
		pass

	try:
		os.remove('Results/Galaxies/'+science_image+'Galaxy_Grid.fits') #removes a file if it has the same name because astropy doesn't like overwriting existing files
		#print 'Old '+science_image+'Galaxy_Grid.fits, new file will be written'
	except OSError:
		pass

	new=numpy.zeros(shape=(res,res))
	hdu=fits.PrimaryHDU(new)
	hdulist=fits.HDUList([hdu])
	hdulist.writeto('Results/Galaxies/Primer_Galaxy_Grid.fits')

	hdulist_new_sci=fits.open('Results/Galaxies/Primer_Galaxy_Grid.fits',ignore_missing_end=True)
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
			

	hdulist_new_sci.writeto('Results/Galaxies/'+science_image+'Galaxy_Grid.fits')
	k.close()
	#print square #DO NOT USE UNTIL MODIFIED

def Scaling(science_image ,xcord, ycord, mag_array, flux_array, background_array, zpt, fake_stars, CCD_Num):
	#f=open('Fake_Star_Catalog/'+science_image[1]+'_Fake_Star_Catalog.dat','w') #--This will need to be file specific at a later date #RESOLVED 29-10-13 15:31
	#print 'lengths', len(mag_array), len(flux_array)
	ranmagarray=[]
	xcord_star=[]
	ycord_star=[]
	newx_star=[]
	newy_star=[]
	mag_array_star=[]
	flux_array_star=[]
	ran_mag_star=[]
	ran_flux_star=[]
	background_array_star=[]
	scaling_factor_star=[]
	CCD_Num_star=[]
	for i in range(0,fake_stars):
		ran_mag=random.uniform(15.0, 23.0) #The fake stars will be in this range of magnitudes
		ran_flux=10.0**((ran_mag-zpt)/(-2.5))
		ranmagarray.append(ran_mag)
		star=int(random.uniform(0,len(xcord)-1))
		
		scaling_factor=((ran_flux)/flux_array[star])
		
		newX=random.uniform(100.0,1948.0)#No edge stars, this creates a 100 pixel border
		newY=random.uniform(100.0, 3996.0)#no edge stars
		#print background_array[star]
		#'Old X', 'Old Y', 'New X', 'New Y', 'Old Mag', 'Old Flux', 'New Mag', 'New Flux', 'Scaling Factor' 		
		#f.write(str(xcord[star])+' '+str(ycord[star])+' '+str(newX)+' '+str(newY)+' '+str(mag_array[star])+' '+str(flux_array[star])+' '+str(ran_mag)+' '+str(ran_flux)+' '+str(background_array[star])+' '+str(scaling_factor)+' '+str(CCD_Num)+'\n')
		xcord_star.append(xcord[star]); ycord_star.append(ycord[star]); newx_star.append(newX); newy_star.append(newY); mag_array_star.append(mag_array[star]); flux_array_star.append(flux_array[star]); ran_mag_star.append(ran_mag); ran_flux_star.append(ran_flux); background_array_star.append(background_array[star]); scaling_factor_star.append(scaling_factor); CCD_Num_star.append(CCD_Num)
		i+=1
	return xcord_star, ycord_star, newx_star, newy_star, mag_array_star, flux_array_star, ran_mag_star, ran_flux_star, background_array_star, scaling_factor_star, CCD_Num_star
	#pylab.hist(ranmagarray, bins=8)
	#pylab.show()
	#f.close()
#----------GALAXY BIASING----------------------------------NEW----------!!!!!!!!!!!!!!!!!!!



def add_fakes_2galaxy(science_image,boxsize, xcord_star, ycord_star, newx_star, newy_star, mag_array_star, flux_array_star, ran_mag_star, ran_flux_star, background_array_star, scaling_factor_star, CCD_Num_star):
	#This step finds the galaxies and adds fake stars to them
	h=open('Results/Galaxies/'+science_image[1]+'_Galaxy_Catalog.cat') #Opens the Galaxy catalog
	f=open('Results/Fake_Star_Catalog/'+science_image[1]+'_Fake_Star_Catalog.dat','w') #Opens the fake star catalog
	reg=open('Results/Fake_Star_Catalog/'+science_image[1]+'_Fakes_Star_Regions.reg','w') #creates region file
	hin=h.readline() #reads first line
	
	fake_star_array=[] #prepares and array for fake stars
	
	for stars in range(0,len(xcord_star)):
		fake_star_array.append(stars) #creates an array [0,1,2,3, etc]

	gal_line_array=[] #prepares and array for galaxies
	while hin: #Adds all galxies into an array
		gal_line_array.append(hin)
		hin=h.readline()
	#while fin: #adds all stars into an array
	#	fake_star_array.append(fin)
	#	fin=f.readline()
	h.close()
	#f.close()

	hdulist_sci= fits.open(science_image[0]+science_image[1]+'.fits',ignore_missing_end=True) #THE SCIENCE DATA THAT WILL OPENED AND MODIFIED
	science_data= hdulist_sci[0].data

	resy=science_data.shape[0]
	resx=science_data.shape[1]

	#print len(fake_star_array), ' Fake Stars have been added to ', len(fake_star_array), ' Galaxies'
	#print gal_line_array[1]
	#print 'Number of Fakes to Be added: ', len(xcord_star)
	num_of_gal_fakes=0
	j=open('Results/Fakes_added/'+science_image[1]+'_Flux_Boxes.dat','w')
	galaxy_mask=numpy.ones((resy,resx),dtype=bool)
	for i in range(0,len(xcord_star)): #Will only add n fake stars to n Galaxies
		#host_galaxy=gal_line_array.pop(random.randrange(0,len(gal_line_array))) #selecting a random host galaxy. Used .pop() so that the same galaxy isnt chosen twice
		
		source_star=fake_star_array.pop(random.randrange(0,len(fake_star_array))) #selecting a random source star. Used .pop() so that the same star isnt chosen twice
			
		
		#print y 
		#print 'len: ',len(gal_line_array)
		#ln=host_galaxy.split()
		#x=float(ln[3])
		#y=float(ln[4])
		

		
		#print 'Lenth of Possible Galaxies: ', len(gal_line_array)
		while len(gal_line_array)>0: #and num_of_gal_fakes<len(xcord_star):
			#host_galaxy=random.choice(gal_line_array)
			#print 'Host Galaxy: ', host_galaxy
			host_galaxy=gal_line_array.pop(random.randrange(0,len(gal_line_array))) #selecting a random host galaxy. Used .pop() so that the same galaxy isnt chosen twice
			ln=host_galaxy.split()
			x=float(ln[3])
			y=float(ln[4])
			
			if galaxy_mask[y,x]==False:
				#print 'Cant Go there'
				continue
			else:

				r=40 #Radius of Mask
				ym,xm = numpy.ogrid[-y:resy-y, -x:resx-x] #Some clever numpy stuff
				mask = xm*xm + ym*ym <= r*r
				galaxy_mask[mask]=False
				
				'Doing Usual Business'
				#
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
				#kn=source_star.split()
				sourcex=xcord_star[source_star] #stars current x location
				sourcey=ycord_star[source_star] #stars current y location
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
				
				scale_fac=scaling_factor_star[source_star] #scale factor
				back=background_array_star[source_star] #background
				
				#---Old area to be scaled---
				startx=int(sourcex-8.0)
				starty=int(sourcey-8.0)		
				finx=int(sourcex+8.0)
				finy=int(sourcey+8.0)
				
				#---New area to have flux added---
				Nstartx=newx-8.0
				Nstarty=newy-8.0
				Nfinx=newx+8.0
				Nfiny=newy+8.0


				newdata=numpy.ones((16,16)) #Preparing a blank gird for scaled objects

				newdata[0:16,0:16]=(((science_data[starty:finy,startx:finx]))-back)*scale_fac #inserting scaled object

				science_data[Nstarty:Nfiny, Nstartx:Nfinx]= (science_data[Nstarty:Nfiny, Nstartx:Nfinx]) + newdata #Modifying the science image
				
				f.write(str(xcord_star[source_star])+' '+str(ycord_star[source_star])+' '+str(newx)+' '+str(newy)+' '+str(mag_array_star[source_star])+' '+str(flux_array_star[source_star])+' '+str(ran_mag_star[source_star])+' '+str(ran_flux_star[source_star])+' '+str(background_array_star[source_star])+' '+str(scaling_factor_star[source_star])+' '+str(CCD_Num_star[source_star])+'\n')
				
				num_of_gal_fakes+=1
				break

		

	try:
		os.remove('Results/Fakes_added/'+science_image[1]+'_WITH_FAKES.fits') #removes a file if it has the same name because astropy doesn't like overwriting existing files
		#print 'Old Fake Removed'
	except OSError:
		pass
		
	hdulist_sci.writeto('Results/Fakes_added/'+science_image[1]+'_WITH_FAKES.fits') #Saving image after loop of 200 Stars is complete
	j.close()
	reg.close()	
	f.close()
	print num_of_gal_fakes, 'fake Stars Added to Galaxies in the Image: ', science_image[1]

	





def Execute(run):
	#print '!!!!!!', run	
	science_image=run
	#print '@@@@@', run[0], run[1]
	#print '!!!!!',science_image[0], science_image[1]

	sci_fil=science_image[0]+science_image[1]+'.fits'
	#maskfile=science_image[0]+science_image[1]+'.weight'
	#print  '######', science_image[1],'.weight'
	print sci_fil
	#print maskfile

	
	#print 'Name: ', science_image

	try:
		hdulist_multi_sci=fits.open(science_image[0]+science_image[1]+'.fits')
		#print '++++ multi_mask assign ', science_image
		
	except IOError or Warning or UnboundLocalError:
		bad_images=open('Results/Bad_Images.dat','a')
		bad_images.write(str(science_image[0])+str(science_image[1])+str('.fits')+' '+str('Reason: Astropy Could not Open the .fits file')+'\n')
		bad_images.close()
		print 'Cant open Science'
		
		
		#shutil.move(science_image+'.fits','Bad_Images/')
		#shutil.move(maskfile+'.fits','Bad_Images/')
		return
	#try:
	#	multi_mask=fits.open(maskfile+'.fits')
	#	multi_mask.close()
	#except IOError or Warning or UnboundLocalError:
	#	'Cant Open Mask'
	#	bad_images=open('Bad_Images.dat','a')
	#	bad_images.write(str(science_image[0])+str(science_image[1])+str('.fits')+'\n')
	#	bad_images.close()
	#	#shutil.move(science_image+'.fits','Bad_Images/')
	#	#shutil.move(maskfile+'.fits','Bad_Images/')
	#	return
		
	zeropoint=float(hdulist_multi_sci[0].header['UB1_ZP'])
	seeing=float(hdulist_multi_sci[0].header['SEEING'])
	saturation=55000.0 #float(hdulist_multi_sci[0].header['SATURATE'])
	gain=float(hdulist_multi_sci[0].header['GAIN'])
	CCD_Num=float(hdulist_multi_sci[0].header['CCDID'])


	fake_stars= 150 #number of fake stars per image (integer please!)

	hdulist_multi_sci.close()
	#weight_map(maskfile, science_image)
	#print '@@@@@ WEIGHTMAP CREATED ', science_image 
	#Adding Weight_Map

	Sextract(science_image,zeropoint,seeing,saturation,gain)
	#print '@@@@@ sextract done'
	catsize=Enough_Objects(science_image)
	#print '@@@@@Enough_Objects Done'
	if catsize==False:
			print science_image, 'didn\'t have enough objects detected so it was moved to Results/Bad_Images/ and the newly created weight map, sex file and catalog have been deleted'
			bad_images=open('Results/Bad_Images.dat','a')
			bad_images.write(str(science_image[0])+str(science_image[1])+str('.fits')+' '+str('Reason: Sextractor did not detect enough objects (<200)')+'\n')
			#shutil.move(science_image+'.fits','Bad_Images/')
			#shutil.move(maskfile+'.fits','Bad_Images/')
			#os.remove('Weight_Map/'+science_image[1]+'_Weight_Map.fits')
			os.remove('Results/Catalog/'+science_image[1]+'_Catalog.cat')
			return
	#else: print '@@@ CATSIZE IS TRUE'
		
	#print 'Sextract done'
	x, y, mag, flux, back = Selecting_Bright(science_image)
	#print '@@@@ SELECTING BRIGHT DONE'
	#print len(x)
	if len(x)<5:
		print 'Not enough Objects met the source star criteria in ', science_image[1]
		bad_images=open('Results/Bad_Images.dat','a')
		bad_images.write(str(science_image[0])+str(science_image[1])+str('.fits')+' '+str('Reason: Less than 5 stars met the source criteria to be a fake')+'\n')
			#shutil.move(science_image+'.fits','Bad_Images/')
			#shutil.move(maskfile+'.fits','Bad_Images/')
			#os.remove('Weight_Map/'+science_image[1]+'_Weight_Map.fits')
		os.remove('Results/Catalog/'+science_image[1]+'_Catalog.cat')
		return
	#print '@@@@@LENGTH IS TRUE'
	xcord_star, ycord_star, newx_star, newy_star, mag_array_star, flux_array_star, ran_mag_star, ran_flux_star, background_array_star, scaling_factor_star, CCD_Num_star=Scaling(science_image, x, y, mag, flux, back, zeropoint, fake_stars, CCD_Num)
	#print '@@@@Scaling done'
	selecting_galaxies(science_image)
	boxsize=[3,5,7]
	add_fakes_2galaxy(science_image,boxsize, xcord_star, ycord_star, newx_star, newy_star, mag_array_star, flux_array_star, ran_mag_star, ran_flux_star, background_array_star, scaling_factor_star, CCD_Num_star)

#-----------------------------------RUN PIPELINE------------------------------------------

file_structure()



all_fits=[] #Establishing an array to find the files
#path=[]
#fnames=[]
for dirpath,dirname,filenames in os.walk(os.path.abspath('../fakes')): #Traverses through a directory tree from 'refs' (refs is a folder in your CWD)
	#print dirname
	#if dirname=='Results':
	#	print '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@Results skipped'
	#	continue
	for file in filenames:
		fileex=os.path.splitext(file)[-1] #Splits the file name, [-1] means it will look at the extension
		if fileex== '.fits': #wanted all .fits files 
			#print dirpath, '$$$', file 
			#complete=os.path.join(dirpath,file) #joins the path and file name (path from 'refs')
			#print complete
			all_fits.append([dirpath, file])
#print all_fits
science_fits=[]
for i in range(len(all_fits)):
	#fname=all_fits[1]
	ln=all_fits[i]
	fname=ln[1].split('.')
	#print fname
	
	
	if fname[-2]=='w':
	
		science_fits.append([ln[0]+str('/'), (os.path.splitext(ln[1])[0])])

#print science_fits

bad_images=open('Results/Bad_Images.dat','w')
bad_images.close()


processors=4
pool=Pool(processors)

pool.map(Execute,science_fits)

pool.close()


