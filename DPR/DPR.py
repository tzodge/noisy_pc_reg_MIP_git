
# for i in range(1,6):

i = 1
# data = 'bunnyLow1000_pts_180_deg_log_{}'.format(i)
# path = '/home/biorobotics/Desktop/tejas/gurobiCodesPython/tejas5_21_jan_clean_up/faceNeglect/datasets_large_number_model_pts'+ '/'
# print data,"data"

data = "david_face1000_noisy_log_1"
path = "/home/biorobotics/Desktop/tejas/gurobiCodesPython/datasets"+"/"
# ####### SWITCHES


set_start_sol_switch = 0 ### 1 sets the initial seed value using ICP/GICP 

callback_switch = 1 #### on or off

from phimax_as_correspondence_modified import main_func  


####### constants
Ns_sampled = 50 ######## fewer sensor points for optimization 
num_sampled_sens_points_ICP = 500
approx_sampled_model_points_ICP = 200
model_tol_band = 20


main_func(data, \
	path, \
	set_start_sol_switch,\
	callback_switch,\
	Ns_sampled,\
	num_sampled_sens_points_ICP,\
	approx_sampled_model_points_ICP,\
	model_tol_band )

### saving in file