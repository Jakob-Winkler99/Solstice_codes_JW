#
# This is an example script to simulate a central receiver system (CRS)
# using Solstice via solsticepy
#
import solsticepy
from solsticepy.master import Master
import numpy as np
import os

#==================================================
# INPUT PARAMETERS

# whether run a new case (True) or load pre-existing input.yaml inputs (False):
new_case=True 

# the sun
# =========

DNI = 1000 # W/m2
sunshape = 'pillbox'
half_angle_deg = 0.2664
sun = solsticepy.Sun(dni=DNI, sunshape=sunshape, half_angle_deg=half_angle_deg)

# S3. sun position
from solsticepy.cal_sun import *
sun_tubes=SunPosition()
day=sun_tubes.days(21, 'Jun')
latitude=37.4425
delta=sun_tubes.declination(day)
daytime, sunrise=sun_tubes.solarhour(delta, latitude)
omega=sunrise+15.*2. # solar noon
theta=sun_tubes.zenith(latitude, delta, omega)
phi=sun_tubes.azimuth(latitude, theta, delta, omega)

azimuth, elevation=sun_tubes.convert_convention(tool='solstice', azimuth=phi, zenith=theta)

num_rays=10000000
#
# the field
# ==========
# F1.Layout
layoutfile='./PS10_layout.csv'
hst_field=True # simulate the full heliostat field (if False: give single heliostat information)

# F2. Heliostat
hst_w=12.84 # m
hst_h=9.45# m
rho_refl=0.88 # mirror reflectivity
slope_error=2.9e-3 # radians
# F3. Tower
tower_h=100.5 # tower height
tower_r=9 # tower radius
#
# the receiver
# ============
# R1. shape
receiver='tubular_cylinder' # 'tubular_cylinder' or 'stl' or ' tubular_cylinder
# R2. Size
gap=0
cylinder_radius=0.05 #radius, m
num_tubes=10 # number of tubes
reciever_width=(num_tubes*(2*cylinder_radius))+((num_tubes-1)*(0*(2*cylinder_radius)))
cylinder_height=reciever_width
rec_w=reciever_width
rec_h=cylinder_height

# R3. tilt angle
tilt=0.  # deg
# R4. position
loc_x=0. # m
loc_y=0 # m
loc_z=108 # m
# R5. Abosrptivity
rec_abs=0.9 # normal value = 0.9

if receiver=='flat':
    # receiver mesh, for binning the flux distribution
    rec_mesh=100
elif receiver=='stl':
    stlfile='./demo_plane.stl'
elif receiver=='tubular_cylinder':
    cylinder_stacks=num_tubes
    cylinder_slices=12

# set the folder name for saving the output files
# False: an automatically generated name  
# or
# string: name of the user defined folder
# (Note: you must set this folder name if you want to use new_case==False)
userdefinedfolder=False


# NO NEED TO CHANGE THE CONTENT BELOW
# ===============================================================
# the ray-tracing scene will be generated 
# based on the settings above

if hst_field:
	# extract the heliostat positions from the loaded CSV file
    layout=np.loadtxt(layoutfile, delimiter=',', skiprows=2)
    hst_pos=layout[:,:3]
    hst_foc=layout[:,3] # F2.5
    hst_aims=layout[:,4:] # F4.
    one_heliostat=False
else:
    one_heliostat=True

if new_case==False:
	assert os.path.isdir(userdefinedfolder)

if userdefinedfolder:
    casefolder=userdefinedfolder
else:
	# define a unique case folder for the user
    snum = 0
    suffix = ""
    while 1:
        import datetime,time
        dt = datetime.datetime.now()
        ds = dt.strftime("%a-%H-%M")
        casefolder = os.path.join(os.getcwd(),'case-%s-%s%s'%(os.path.basename(__file__),ds,suffix))
        if os.path.exists(casefolder):
            snum+=1
            suffix = "-%d"%(snum,)
            if snum > 200:
                raise RuntimeError("Some problem with creating casefolder")
        else:
            # good, we have a new case dir
            break

if receiver =='flat':
    rec_param=np.r_[rec_w, rec_h, num_tubes, rec_mesh, loc_x, loc_y, loc_z, tilt]
elif receiver =='stl':
    rec_param=np.array([rec_w, rec_h, stlfile, loc_x, loc_y, loc_z, tilt])
elif receiver =='tubular_cylinder':
    rec_param=np.r_[num_tubes, cylinder_height, cylinder_radius, reciever_width, cylinder_slices, cylinder_stacks, loc_x, loc_y, loc_z, tilt, tower_h, gap]

master=Master(casedir=casefolder)
outfile_yaml = master.in_case(folder=casefolder, fn='input.yaml')
outfile_recv = master.in_case(folder=casefolder, fn='input-rcv.yaml')

if new_case:
	# generate the YAML file from the input parameters specified above
    solsticepy.gen_yaml(sun, hst_pos, hst_foc, hst_aims,hst_w, hst_h
		, rho_refl, slope_error, receiver, rec_param, rec_abs
		, outfile_yaml=outfile_yaml, outfile_recv=outfile_recv
		, hemisphere='North',tower_h=tower_h, tower_r=tower_r, spectral=False
		, medium=0, one_heliostat=one_heliostat)

#   they were in the above command

# run Solstice using the generate inputs, and run all required post-processing
master.run(azimuth, elevation, num_rays, rho_refl,sun.dni, folder=casefolder, gen_vtk=True)

# annual solution (see instructions)
#master.run_annual(nd=5, nh=5, latitude=latitude, num_rays=num_rays, num_hst=len(hst_pos),rho_mirror=rho_refl, dni=DNI)

