import numpy as np 
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle
from scipy.spatial import KDTree
import scipy
from icp_simple import icp_test
import xlwt
import openpyxl		
import transforms3d #https://matthew-brett.github.io/transforms3d/reference/transforms3d.euler.html#eulerfuncs

from gurobipy import *




######### functions
def Bt_Nsx6_to_3x3_UD (B_nx6):
	B_return = np.zeros((B_nx6.shape[0] , 3,3))
	for ii in range(0,B_nx6.shape[0]):
		temp = np.zeros((3,3))
		temp[0,0] 			= B_nx6[ii,0] 
		temp[0,1],temp[1,0] = B_nx6[ii,1], 0 
		temp[0,2],temp[2,0] = B_nx6[ii,2], 0 
		temp[1,1] 			= B_nx6[ii,3] 
		temp[1,2],temp[2,1] = B_nx6[ii,4], 0 
		temp[2,2] 			= B_nx6[ii,5] 
		B_return[ii] = temp
		

	return B_return  ######## B_return nx3x3

def Bt_Nsx6_to_3x3 (B_nx6):
	B_return = np.zeros((B_nx6.shape[0] , 3,3))
	temp = np.zeros((3,3))
	temp[0,0] 			= B_nx6[ii,0] 
	temp[0,1],temp[1,0] = B_nx6[ii,1],  B_nx6[ii,1] 
	temp[0,2],temp[2,0] = B_nx6[ii,2], B_nx6[ii,2] 
	temp[1,1] 			= B_nx6[ii,3] 
	temp[1,2],temp[2,1] = B_nx6[ii,4],  B_nx6[ii,4] 
	temp[2,2] 			= B_nx6[ii,5] 
	B_return[ii] = temp

	return B_return  ######## B_return nx3x3



def R_T_to_4X4 ((R,T)):
	T = np.reshape(T,(3,1))
	X = np.concatenate((R,T), axis = 1)
	X = np.concatenate((X,np.array([[0, 0, 0, 1]])), axis = 0)
	return X


def _4X4_to_R_T ((Transf)):
	R = Transf [0:3, 0:3]
	T = Transf [0:3, 3]
	return R,np.reshape(T,(3,1))




def array_to_list (arr):
	output_arr = arr.flatten()
	output_list = output_arr.tolist() 
	return output_list, len(output_list)


def find_Hb (transformed_sens_points):
	Hb = np.zeros((transformed_sens_points.shape[0], model_tol))
	for ii in range (transformed_sens_points.shape[0]):
		d_min = 1e8
		for jj in range(model_tol):
			d_new = np.sum(abs(B_combined[Q[ii,jj]].dot(transformed_sens_points[ii,:] - M[:,Q[ii,jj]])))
			if d_new < d_min:
				d_min = d_new	
				Hb[ii,:] = 0
				Hb[ii,jj] = 1
	return Hb


def find_lam(R):
	lam = np.zeros((num_partitions_SOS2,3,3))
	q = np.linspace(-1,1,num_partitions_SOS2)
	for ii in range(3):
		for jj in range(3):
			for kk in range (num_partitions_SOS2-1):
				if q[kk] <= R[ii,jj] and  R[ii,jj] <= q[kk+1]:
					lam[kk][ii,jj] = (R[ii,jj]-q[kk+1])/(q[kk]-q[kk+1])
					lam[kk+1][ii,jj] = 1 - (R[ii,jj]-q[kk+1])/(q[kk]-q[kk+1])
	return lam


def find_w (lam):
	w = np.zeros((3,3))
	for ii in range(3):
		for jj in range(3):
			sum1 = 0
			for kk in range (num_partitions_SOS2):
				sum1 = sum1 + lam[kk,ii,jj]*q[kk]*q[kk]
			w[ii,jj] = sum1
	return w



def find_beta(R,T,S,M):
	Ns_temp = S.shape[1]

	beta = np.zeros((Ns_temp,model_tol,3))
	for ii in range(Ns_temp):
		for jj in range(model_tol):
			beta[ii,jj,:] = abs(B_combined[Q[ii,jj]].dot(-M[:,Q[ii,jj]] + R.dot(S[:,ii]) + T.reshape(3,)))

	return beta		
			


def find_phi(beta,H):
	Ns_temp = beta.shape[0]
	
	phi = np.zeros(( Ns_temp,model_tol ))
	
	for jj in range(model_tol):
		for ii in range(Ns_temp):
			phi[ii,jj] = np.sum(beta[ii,jj,:]) - (1 - H[ii,jj])*bigM
			if phi[ii,jj] < 0:
				phi[ii,jj] = phiMax 


	return phi	



def find_all_opt_variables(sens_points, model_points, model_tree, R, T):
	### here R and T are from sensor to model M = RS + T
	### sens_points shape nx3
	### model_points shape nx3
	T = T.reshape(3,1)
	transformed_sens_points = (np.matmul(R,sens_points.T) + T)
	H_start = find_Hb(transformed_sens_points.T)
	lam_start  = find_lam(R)
	w_start = find_w(lam_start)
	beta_start = find_beta(R,T,sens_points.T,model_points.T)
	phi_start = find_phi(beta_start, H_start) 
	return H_start, lam_start, w_start, phi_start, beta_start 


def add_McCormick(m,x,x_lb,x_ub,y,y_lb,y_ub,v,temp0,temp1,temp2,temp3):
	m.addConstr( temp0 == (x_lb*y + y_ub*x - x_lb*y_ub))
	m.addConstr( temp1 == (x_ub*y + y_lb*x - x_ub*y_lb))
	m.addConstr( temp2 == (x_lb*y + y_lb*x - x_lb*y_lb))
	m.addConstr( temp3 == (x_ub*y + y_ub*x - x_ub*y_ub))
	m.addConstr(v >= temp2)
	m.addConstr(v >= temp3)
	m.addConstr(v <= temp1)
	m.addConstr(v <= temp0)
	return 0


def chol_sigma_inv (SigmaS_6x1,SigmaMsampled_6x1):
	SigmaSi_3x3	= Bt_Nsx6_to_3x3 (SigmaS_6x1.reshape(1,6))
	SigmaMsampled_3x3	=Bt_Nsx6_to_3x3 (SigmaMsampled_6x1)
	C_inv = np.linalg.inv(SigmaSi_3x3+SigmaMsampled_3x3)
	B_cov = np.linalg.cholesky(C_inv).T
	return B_cov


def main_func(data, \
	path, \
	set_start_sol_switch,\
	callback_switch,\
	Ns_sampled,\
	num_sampled_sens_points_ICP,\
	approx_sampled_model_points,\
	model_tol_1,\
	TimeLimit,\
	corr):

	S_given = np.loadtxt(path +data +'/S.txt',delimiter = ',') ##file is nx3
	M_given = np.loadtxt(path +data +'/Msampled.txt',delimiter = ',')##file is nx3
	gt_given = np.loadtxt(path +data +'/gt.txt',delimiter = ',') ##file is 4x4 
	B_combined_given = np.loadtxt(path +data +'/B_combined.txt',delimiter = ',') ##file is NfxNm

	global S,M,gt,num_partitions_SOS2,q,phiMax,bigM,B_combined,Q,model_tol 
	model_tol = model_tol_1
	bigM = 1e+4
	phiMax = 1000
	IntFeasTol =1e-9
	OptimalityTol = 1e-9

	S = S_given.T  # S 3xNs  sensor points
	M = M_given.T  # M 3XNm  Model points
	gt = gt_given  # 4x4 

	correspondenceIdx = corr
	gt_inv = np.linalg.inv(gt)
	B_combined_Nsx6 = B_combined_given    # Nmx6 
	B_combined = Bt_Nsx6_to_3x3_UD (B_combined_Nsx6)  # B_combined = Nmx3x3
	
	print S.shape, "S.shape"
	print M.shape, "M.shape"
	print gt, "gt"

	Ns = S.shape[1]
	Nm = M.shape[1]
	num_partitions_SOS2 = 50

	interval_sampled_model_points = Nm/approx_sampled_model_points
	M_sampled = M[:,0::interval_sampled_model_points] ######### remove this later if not required

	Nm_sampled = M_sampled.shape[1]
	print M.shape, "M.shape" 
	print M_sampled.shape,  "M_sampled.shape" 


	tree_M = KDTree(M.T) 
	tree_M_sampled = KDTree(M_sampled.T)



#########################################################################
############# creating the band matrix ##################################
#########################################################################

	Q = np.zeros((Ns_sampled,model_tol),dtype=int)
	for ii in range(Ns_sampled):
		Q[ii,:] = tree_M.query(M[:,correspondenceIdx[ii]],model_tol)[1]
	# print Q,"Q"
	# print correspondenceIdx
	# print "............................>>>>>>>>>>...>>>>>>>....>>>>.>>>...>>..>>>........"
#########################################################################




	m = Model("tejas5.2")

	###### gurobi parameters
	# m.Params.BestObjStop = phi_minimum.sum()/Ns_sampled
	m.Params.TimeLimit = TimeLimit
	m.Params.OptimalityTol = OptimalityTol
	m.Params.IntFeasTol =  IntFeasTol 
	# number_of_solutions = 100 
	# m.Params.PoolSolutions = number_of_solutions

	errorMat = gt
	angleError = transforms3d.axangles.mat2axangle(errorMat)[1]*180/np.pi
	positionError = errorMat[0:3,3]
	print ".........................."
	print "error before optimisation"
	print np.linalg.norm(angleError), "angleError"
	print np.linalg.norm(positionError), "positionError" 


	plot_switch = 0
	if plot_switch == 1:
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		plt.title('before optimisation')
		ax.scatter(M_sampled[0,:],M_sampled[1,:],M_sampled[2,:],s = 50, c= 'red',alpha = 0.5, label = 'M_sampled')
		ax.scatter(S[0,:],S[1,:],S[2,:], s = 20, c= 'blue', label ='S_given', alpha = 1)
		plt.legend()
		plt.show()
	 



	####### Gurobi variables
	T = m.addVars(3,1 , lb=-GRB.INFINITY, ub = GRB.INFINITY)
	R = m.addVars(3,3 , lb=-1, ub = 1, name='rotationMatrix')
	w = m.addVars(3,3 , lb=0, ub = 1)
	q = np.linspace(-1,1,num_partitions_SOS2)
	lam = m.addVars(num_partitions_SOS2,3,3, lb = 0, name='lam') #lambda
	phi =  m.addVars(Ns_sampled,model_tol , lb=0, ub = GRB.INFINITY)
	Hb_sampled = m.addVars(Ns_sampled, model_tol, vtype=GRB.BINARY )
	beta = m.addVars(Ns_sampled, model_tol,3 , lb = 0, ub = GRB.INFINITY)
	 


	m.update()


	if set_start_sol_switch == 1:

		m.Params.StartNodeLimit = 2000000000-2

		
		# bunnyLow6000_pts_180_deg_log_1_results
		# R_start = Tejas1_results["R"]
		# T_start = Tejas1_results["T"]

		# Hb_sampled_start, lam_start, w_start, phi_start, beta_start 	=\
		# find_all_opt_variables(S[:,0:Ns_sampled].T, M_sampled[:,0:Nm_sampled].T, tree_M_sampled, R_start, T_start)
		# print R_start,"R_start"
		# print T_start,'T_start'
		# print np.linalg.det(R_start.dot(R_start.T)), "RR' det"
		# print ""
		for ii in range (3):
			T[ii,0].start = gt[ii,3]
			for jj in range(3):
				R[ii,jj].start = gt[ii,jj]
				# a= 8
		# Hb_sampled.start = Hb_sampled_start
		# for ii in range (Ns_sampled):
		# 	for jj in range (Nm_sampled):
		# 		Hb_sampled[ii,jj].start =Hb_sampled_start[ii,jj] 

		# lam.start = lam_start
		# for kk in range(num_partitions_SOS2):
		# 	for ii in range(3):
		# 		for jj in range(3):
		# 			lam[kk,ii,jj].start = lam_start[kk,ii,jj] 

		# for ii in range(3):
		# 	for jj in range(3):
		# 		w[ii,jj].start = w_start[ii,jj] 

		# for ii in range(Ns_sampled):
		# 	for jj in range(Nm_sampled):
		# 		for kk in range(3):
		# 			beta[ii,jj,kk].start = beta_start[ii,jj,kk]
				
		# for ii in range(Ns_sampled):
		# 	for jj in range(Nm_sampled):
		# 		phi[ii,jj].start = phi_start[ii,jj] 
######################################################################################################################################


################################################
##########	setting R values  ##################
################################################

	# for i in range(3):
	# 	for j in range(3):
	# 		if i==j:
	# 			m.addConstr(R[i,j] ==1, name='rotConstr2u{},{}'.format(i,j))
	# 		else:
	# 			m.addConstr(R[i,j] ==0, name='rotConstr2u{},{}'.format(i,j))
	# ##set T value
	# for l in range(3):
	# 	m.addConstr(T[l,0] == 0)	

	# gt_inv = np.linalg.inv(gt)
	# for i in range(3):
	# 	for j in range(1):
	# 		m.addConstr(R[i,j] ==gt[i,j], name='rotConstr2u{},{}'.format(i,j))

	# ##set T value
	# for l in range(3):
	# 	m.addConstr(T[l,0] == gt[l,3])	




	###################


	####################################################################3
	######### constraints

	# t0 =  time.time() 

###### relaxed rotation constraints
	# # #Constraints on rotation matrix

	# # sum1mation R(i)^2 <= 1
	for ii in range(3):
		sum1 = 0
		for jj in range(3):
			sum1 = sum1 + R[jj,ii]*R[jj,ii]
		m.addConstr((sum1 <= 1 ) , name='rotConstr{},{}'.format(ii,jj))


	for ii in range(3):	
		for jj in range(ii+1,3):
			sum1 = 0
			for kk in range (0,3):
				sum1 = sum1 + (R[kk,ii] + R[kk,jj])*(R[kk,ii] + R[kk,jj])
			m.addConstr((sum1 <= 2 ) , name='rotConstr1{},{}'.format(ii,jj))

	for ii in range(3):
		for jj in range(ii+1,3):
			sum1 = 0
			for kk in range (0,3):
				sum1 = sum1 + (R[kk,ii] - R[kk,jj])*(R[kk,ii] - R[kk,jj])
			m.addConstr((sum1 <= 2 ) , name='rotConstr1{},{}'.format(ii,jj))


	for tmp1 in [-1,1]:
		for tmp2 in [-1,1]:
			sum1 = 0
			for ii in range(0,3):
				sum1 = sum1 + (R[ii,0] + tmp1* R[ii,1] + tmp2 * R[ii,2])*(R[ii,0] + tmp1* R[ii,1] + tmp2 * R[ii,2])
			m.addConstr((sum1 <= 3) , name='rotConstr1{}'.format(ii))

######################################
	v = m.addVars(1,6,lb= -1, ub =1)
	temp_var = m.addVars(4,6,lb = -3, ub =3)
	R_ub = 1.0
	R_lb = -1.0

 	for i in range(6):
 		m.addConstr(temp_var[0,i] >= -1)
 		m.addConstr(temp_var[1,i] >= -1)
 		m.addConstr(temp_var[2,i] <=  1)
 		m.addConstr(temp_var[3,i] <=  1)


 	add_McCormick(m,R[0,2],R_lb,R_ub,R[1,1],R_lb,R_ub,v[0,0],temp_var[0,0],temp_var[1,0],temp_var[2,0],temp_var[3,0])
 	add_McCormick(m,R[0,1],R_lb,R_ub,R[1,2],R_lb,R_ub,v[0,1],temp_var[0,1],temp_var[1,1],temp_var[2,1],temp_var[3,1])
 	add_McCormick(m,R[1,0],R_lb,R_ub,R[0,2],R_lb,R_ub,v[0,2],temp_var[0,2],temp_var[1,2],temp_var[2,2],temp_var[3,2])
 	add_McCormick(m,R[0,0],R_lb,R_ub,R[1,2],R_lb,R_ub,v[0,3],temp_var[0,3],temp_var[1,3],temp_var[2,3],temp_var[3,3])
 	add_McCormick(m,R[1,0],R_lb,R_ub,R[0,1],R_lb,R_ub,v[0,4],temp_var[0,4],temp_var[1,4],temp_var[2,4],temp_var[3,4])
 	add_McCormick(m,R[0,0],R_lb,R_ub,R[1,1],R_lb,R_ub,v[0,5],temp_var[0,5],temp_var[1,5],temp_var[2,5],temp_var[3,5])

	
 	m.addConstr(-v[0,0] + v[0,1] == R[2,0])
 	m.addConstr( v[0,2] - v[0,3] == R[2,1])
 	m.addConstr(-v[0,4] + v[0,5] == R[2,2])

#######################################


	# ######## W constraints
	sum1 = 0
	for ii in range(3):
		sum1 = 0
		for jj in range(3):
			sum1 = sum1 + w[ii,jj]
		m.addConstr(( sum1 >= 1) , name = 'w' )



	# ##SOS constraints  
	# #### R and W constraints
	for ii in range (3):
		for jj in range(3):
			sum1 = 0
			sum2 = 0
			sum3 = 0
			for kk in range(num_partitions_SOS2):
				sum1 = sum1 + lam[kk,ii,jj]
				sum2 = sum2 + lam[kk,ii,jj]*q[kk]
				sum3 = sum3 + lam[kk,ii,jj]*q[kk]*q[kk]
			m.addConstr(sum1 == 1 )
			m.addConstr(R[ii,jj] == sum2)
			m.addConstr(w[ii,jj] == sum3)


	#########SOS2 constrain
	lam_list = [0]*num_partitions_SOS2
	weight_list =  [0]*num_partitions_SOS2
	for ii in range (3):
		for jj in range(3):
			for k in range(num_partitions_SOS2):	
				lam_list[k] = lam[k,ii,jj]
				weight_list[k] = k+1
			m.addSOS(GRB.SOS_TYPE2, lam_list, weight_list)
			





	for jj in range(model_tol): 	
		for kk in range(3):
			for ii in range(Ns_sampled):
				sum1 = 0
				for ll in range(3):
					sum1 = sum1 + B_combined[Q[ii,jj],kk,ll]*M[ll,Q[ii,jj]]    
				
				sum2 = 0
				for ll in range(3):
					sum4 = 0
					for mm in range(3):
						sum4 = sum4 + R[ll,mm]*S[mm,ii]
					sum2 = sum2 + B_combined[Q[ii,jj],kk,ll]*sum4   	

				sum3 = 0
				for ll in range(3):
					sum3 = sum3 + B_combined[Q[ii,jj],kk,ll]*T[ll,0] 	
 
 				m.addConstr((beta[ii,jj,kk] >=  ( sum1 - sum2 -sum3)  ), name="betaplus,{},{}".format(ii,jj) )
				m.addConstr((beta[ii,jj,kk] >= -( sum1 - sum2 -sum3)  ), name="betaminus,{},{}".format(ii,jj) ) 
				

	for ii in range(Ns_sampled):
		for jj in range(model_tol):
			sum6 = 0
			for kk in range(3):
				sum6 = sum6 + (beta[ii,jj,kk])
			m.addConstr(phi[ii,jj] >= sum6 - bigM*(1-Hb_sampled[ii,jj])) 
			m.addConstr(phi[ii,jj] >= phiMax*(1-Hb_sampled[ii,jj]))			
                

	for ii in range(Ns_sampled):
		sum7 = 0
		for jj in range(model_tol):
			sum7 = sum7 + Hb_sampled[ii,jj]
		m.addConstr(sum7 == 1)


	m.setObjective(phi.sum() - abs(float(model_tol-1)*(Ns_sampled)*phiMax) , GRB.MINIMIZE)


	# print "time taken = " ,time.time() - t0 
	global sol_repeat_counter 	
	global solflag 	
	solflag  = 0
	sol_repeat_counter = 0

	def callbackStoM(m, where):
		global sol_repeat_counter 	
		global solflag 
		# print "aa gaya"
		t_callback_start = time.time()
		if 	sol_repeat_counter < 100 :
			# print "aa gaya"

			temp_T = np.zeros((1,3))
			temp_R = np.zeros((3,3))
			# R_final = np.zeros((3,3))  
			# T_final = np.zeros((1,3))

			if where == GRB.Callback.MIPNODE:

				sol_at_node =  m.cbGetNodeRel(m._vars)
				l = 0
				start_count = 0
				for i in range(start_count , start_count +3): ###Extracting R and T from current node
					temp_T[0][l] = sol_at_node[i]
					temp_R[0][l] = sol_at_node[i+3]
					temp_R[1][l] = sol_at_node[i+6]
					temp_R[2][l] = sol_at_node[i+9]
					l = l+1
				
				###### SVD for valid transformation 		
				U_svd, S_svd, Vt_svd = np.linalg.svd(temp_R)
				R_init_ICP = np.dot(U_svd,Vt_svd)
				if np.linalg.det(R_init_ICP)<0:
					Vt_svd [2,:] *= -1
					R_init_ICP = np.dot(U_svd,Vt_svd)
				

				T_init_ICP = temp_T.T
			
				############# ICP applied on num_sampled_sens_points_ICP sensor Points
				transformed_for_ICP = np.matmul(R_init_ICP, S[:,0:num_sampled_sens_points_ICP]) + T_init_ICP 
				
				# R_ICP_output, T_ICP_output  = icp_test(V.T, transformed_for_ICP.T, tree_M, SigmaS, F, points_per_face)
				# ICP_or_GICP_switch = 2
				# ICP_triangle_proj_switch = 0 
				# R_ICP_output, T_ICP_output = icp_test( V.T,transformed_for_ICP.T, M.T, tree_M, M_sampled.T, tree_M_sampled,SigmaS,F,points_per_face,ICP_triangle_proj_switch_callback,ICP_or_GICP_switch_callback,Q)
				# t = time.time()
				R_ICP_output, T_ICP_output = icp_test( transformed_for_ICP.T, M_sampled.T, tree_M_sampled)
				 
									 

				
				total_transf_after_ICP_S_to_M = np.matmul(R_T_to_4X4((R_ICP_output, T_ICP_output)), R_T_to_4X4((R_init_ICP,T_init_ICP))) 				
 
				
				T_after_ICP_S_to_M = np.reshape(total_transf_after_ICP_S_to_M[0:3, 3],(3,1))
				R_after_ICP_S_to_M = np.reshape(total_transf_after_ICP_S_to_M[0:3,0:3],(3,3))
				
 
				
				Hb_sampled_after_ICP, lam_after_ICP, w_after_ICP, phi_after_ICP, beta_after_ICP =\
					find_all_opt_variables( S[:,0:Ns_sampled].T, M.T, tree_M, total_transf_after_ICP_S_to_M[0:3,0:3], total_transf_after_ICP_S_to_M[0:3,3])
				 

				list_T_after_ICP_S_to_M, len_list_T_after_ICP_S_to_M = array_to_list(T_after_ICP_S_to_M)
				list_R_after_ICP_S_to_M, len_list_R_after_ICP_S_to_M = array_to_list(R_after_ICP_S_to_M)
				list_Hb_sampled_after_ICP, len_list_Hb_sampled_after_ICP = array_to_list(Hb_sampled_after_ICP)
				list_w_after_ICP, len_list_w_after_ICP = array_to_list(w_after_ICP)
				list_phi_after_ICP, len_list_phi_after_ICP = array_to_list(phi_after_ICP)
				list_beta_after_ICP, len_list_beta_after_ICP = array_to_list(beta_after_ICP)
				list_lam_after_ICP, len_list_lam_after_ICP = array_to_list(lam_after_ICP)

				###T
				start_count = 0 ## counting for translation vector elements # 3*Ns for Alpha, Ns for Phi
				m.cbSetSolution(m._vars[start_count : start_count + len_list_T_after_ICP_S_to_M ], list_T_after_ICP_S_to_M )
				###R
				start_count = 3 
				m.cbSetSolution(m._vars[start_count : start_count + len_list_R_after_ICP_S_to_M ], list_R_after_ICP_S_to_M )
				###lambda
				start_count = 3 + 3*3 + 3*3
				m.cbSetSolution(m._vars[start_count : start_count + len_list_lam_after_ICP], list_lam_after_ICP)
				###phi
				# start_count = 3 + 3*3 + 3*3 + 3*3*num_partitions_SOS2
				# m.cbSetSolution(m._vars[start_count : start_count + len_list_phi_after_ICP ], list_phi_after_ICP)
				# ###Hb_sampled
				start_count = 3 + 3*3 + 3*3 + 3*3*num_partitions_SOS2 + Ns_sampled*Nm_sampled
				m.cbSetSolution(m._vars[start_count : start_count + len_list_Hb_sampled_after_ICP ], list_Hb_sampled_after_ICP )
				# # ####beta
				# start_count =  3 + 3*3 + 3*3 + 3*3*num_partitions_SOS2 + Ns_sampled*Nm_sampled*2  
				# m.cbSetSolution(m._vars[start_count : start_count + len_list_beta_after_ICP], list_beta_after_ICP)
			 				
				objval = m.cbUseSolution()
				
 				print phi_after_ICP.sum() - abs(float(model_tol-1)*(Ns_sampled)*phiMax) , "obj funct value using heuristics"
 				if  objval == 1e+100 :							
					sol_repeat_counter = sol_repeat_counter + 1 

				else:
					sol_repeat_counter = 0
				print objval
				print 'sol_repeat_counter: ', sol_repeat_counter
				print time.time()-t_callback_start,"time for callback"
				print "................................................................................................................"


	m.update()
	m._vars = m.getVars()
	if callback_switch == 1:
		m.optimize(callbackStoM)
	else:
		m.optimize()


	######
 
	number_of_solutions = 100 
	m.Params.PoolSolutions = number_of_solutions
	m.Params.IntFeasTol = 1e-5
 	# m.Params.SubMIPNodes = 2000


	FivePointEstimation = np.empty((number_of_solutions,4))
	# wb = openpyxl.Workbook()
	# ws = wb.active
 

	t5_2_ang_err = []
	t5_2_tra_err = []
	t5_2_tf = []
	t5_2_corr_python = []
	t5_2_corr_matlab = []
	t5_2_obj_val = m.objVal

	for iii in range(number_of_solutions):
		Bucket = np.zeros((m.SolCount,4,4))
		if iii < m.SolCount:
			m.Params.SolutionNumber = iii 
			print ""
			print ""
			print "bucket number = ", iii

			Hb_sampled_before_bucket = np.zeros((Ns_sampled,model_tol))
			for k in range (0,Ns_sampled):
				for j in range(0,model_tol):
					Hb_sampled_before_bucket[k,j]  =  Hb_sampled[k,j].Xn

			R_bucket_input_S_to_M = np.zeros((3,3))
			for k in range (0,3):
				for j in range (0,3):
					R_bucket_input_S_to_M[k,j] =  R[k,j].Xn

			T_bucket_input_S_to_M = np.zeros((3,1))
			for k in range(0,3):
				T_bucket_input_S_to_M[k,0] = T[k,0].Xn

			# phi_before_bucket = np.zeros((Ns_sampled,model_tol))	
			# for i in range(Ns_sampled):
			# 	for j in range(model_tol):
			# 		phi_before_bucket[i,j] = phi[i,j].Xn
			
			# beta_before_bucket = np.zeros((Ns_sampled,model_tol,3))
			# for i in range(Ns_sampled):
			# 	for j in range(model_tol):
			# 		for k in range(3):
			# 			beta_before_bucket[i,j,k] = beta[i,j,k].Xn	



			U_svd, S_svd, Vt_svd = np.linalg.svd(R_bucket_input_S_to_M)
			R_bucket_input_valid = np.dot(U_svd,Vt_svd)
		
			if np.linalg.det(R_bucket_input_valid)<0:
				Vt_svd [2,:] *= -1
				R_bucket_input_valid = np.dot(U_svd, Vt_svd)



			corr_mat_ids = np.zeros((Ns_sampled,1),dtype = int)
			for i in range(Ns_sampled):
				corr_mat_ids[i] = Q[i,np.where((Hb_sampled_before_bucket[i,:] <=1.2) & (Hb_sampled_before_bucket[i,:] >= 0.9))]
			
			# if iii == 0:
			error_mat = R_T_to_4X4((R_bucket_input_valid,T_bucket_input_S_to_M)).dot(np.linalg.inv(gt))
			angleError = np.asarray(transforms3d.axangles.mat2axangle(error_mat[0:3,0:3]))[1]*180/np.pi
			positionError = error_mat[0:3,3]
			print np.linalg.norm(angleError),"angleError_valid"
			print np.linalg.norm(positionError),"positionError_valid"


			t5_2_ang_err.append(np.linalg.norm(angleError))
			t5_2_tra_err.append(np.linalg.norm(positionError))
			t5_2_tf.append(R_T_to_4X4((R_bucket_input_valid,T_bucket_input_S_to_M)))
			t5_2_corr_python.append(corr_mat_ids)
			t5_2_corr_matlab.append(corr_mat_ids+1)
 		
 			print ""
			# print R_bucket_input_S_to_M, "invalid R"
			print R_bucket_input_valid,"valid"
			print T_bucket_input_S_to_M,"T_bucket_input" 
			print (corr_mat_ids+1).T,"corr_mat_ids, 1 already added"
			print "...."
 			print ""

		else:
			break

	# wb = openpyxl.load_workbook("tejas5_test_21_jan.xlsx")
	# ws = wb.active
	# ws.append(bucket_list_ang_err)

	# wb.save("tejas5_test_21_jan.xlsx")

	return 	t5_2_ang_err,\
			t5_2_tra_err,\
			t5_2_tf,\
			t5_2_corr_python,\
			t5_2_corr_matlab,\
			t5_2_obj_val

