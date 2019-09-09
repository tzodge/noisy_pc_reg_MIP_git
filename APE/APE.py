
set_start_sol_switch = 0
callback_switch = 1



St_as_M_switch = 0  #### sets transformed sensor points as model points



outlier_detection_switch = 2

if outlier_detection_switch == 1:
	from alpha_with_outlier_detection import main_func

elif outlier_detection_switch ==2 :	
	from  alpha_S_to_M_without_Sigma import main_func


data = "david_face1000_noisy_log_1"
path = "/home/biorobotics/Desktop/tejas/gurobiCodesPython/datasets"+"/"

i=1 
# data = 'bunny1000_pts_180_deg_log_{}'.format(i)
# path = '/home/biorobotics/Desktop/tejas/gurobiCodesPython/datasets/partial'+ '/'

print "data being used :", data

Ns_sampled = 20 #number of points used for pose estimation
num_sampled_sens_points_ICP = 1000
approx_sampled_model_points = 500


main_func(data, \
	path, \
	set_start_sol_switch,\
	callback_switch,\
	Ns_sampled,\
	num_sampled_sens_points_ICP,\
	approx_sampled_model_points,\
	)

