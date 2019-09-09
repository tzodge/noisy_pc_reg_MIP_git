from gurobipy import *
import numpy as np 
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from icp_simple import icp_test
from bucket_refinement import bucket_refinement
from scipy.spatial import KDTree
import transforms3d #https://matthew-brett.github.io/transforms3d/reference/transforms3d.euler.html#eulerfuncs
import openpyxl
 


def find_Cb (transformed_sens_points, num_model_points ,tree_of_model ):
	### transformed_sens_points Nsx3
	### returns Cb matrix of dim [num(transformed_sens_points) X num(tree_of_model)]
	Cb = np.zeros((transformed_sens_points.shape[0], num_model_points))
	neighbor_idx = tree_of_model.query(transformed_sens_points)[1]
	for ii in range (transformed_sens_points.shape[0]):
		Cb[ii,neighbor_idx[ii]] = 1
	return Cb,neighbor_idx


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



def find_alpha(R,T,S,M,C):  
## dimensions --> S 3xn,
## dimensions --> T 3x1, 
## dimensions --> C S.shape[1]xM.shape[1], 
	alpha = np.zeros((3,S.shape[1]))
	T = np.reshape(T,(3,))
	for ii in range (S.shape[1]):
		# alpha[:,ii] = abs(B[ii].dot(S[:,ii] - T - R.dot(M).dot(C[ii,:].T)))
		alpha[:,ii] = abs(R.dot(S[:,ii]) + T - M.dot(C[ii,:].T))
			
	return alpha


def find_phi(alpha):
	phi = alpha.sum(axis = 0)
	return np.reshape(phi,(1,alpha.shape[1]))	


def find_all_opt_variables(sens_points, model_points, model_tree, R, T):
	
	transformed_sens_points = (np.matmul(R,sens_points.T)+T)
	Cb_start,_ = find_Cb(transformed_sens_points.T, model_points.shape[0], model_tree)
	lam_start  = find_lam(R)
	w_start = find_w(lam_start)
	# print sens_points.T , "sens_points.T given to find_alpha"
	alpha_start = find_alpha(R,T,sens_points.T, model_points.T,Cb_start)
	phi_start = find_phi(alpha_start) 

	return Cb_start, lam_start, w_start, alpha_start, phi_start 

def R_T_to_4X4 ((R,T)):
	T= np.reshape(T,(3,1))
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

# def save2Ddata(f,var, name ):
# 	num_rows = var.shape[0]	
# 	num_cols = var.shape[1]
# 	f.write(name+"\n")
# 	f.write("{}\n".format(num_rows))
# 	for i in range(num_rows):
# 		for j in range(num_cols):
# 			if j == num_cols-1:	
# 				f.write("{}\n".format(var[i,j] ))
# 			else:
# 				f.write("{} ".format(var[i,j] ))
	
# def save1Ddata(f,var):
# 	num_rows = var.shape[0]	
# 	num_cols = var.shape[1]
# 	for i in range(num_rows):
# 		for j in range(num_cols):
# 			f.write("{} ".format(var[i,j]))
# 	f.write("		")
			
def transf_mat_to_euler_6x1(transf_mat):
	angle_sxyz = np.asarray(transforms3d.euler.mat2euler(transf_mat))
	transl_xyz = np.array(transf_mat[0:3,3])
	return np.append(angle_sxyz,transl_xyz).reshape(1,6)

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



def main_func(data, \
	path, \
	set_start_sol_switch,\
	callback_switch,\
	Ns_sampled,\
	num_sampled_sens_points_ICP,\
	approx_sampled_model_points,\
	approx_sampled_model_points_ICP,\
	TimeLimit):
	t_0 = time.time()
	
	global S,V,M,gt,B,SigmaS,F,num_partitions_SOS2,q

	S_given = np.loadtxt( path+data +'/S.txt',delimiter = ',') # S_given Ns_given X 3	
	M_given = np.loadtxt(path +data +'/Msampled.txt',delimiter = ',') #M_given Nm_given X 3
	gt = np.loadtxt(path +data +'/gt.txt',delimiter = ',') # gt 4 X 4

	M = M_given.T
	S = S_given.T
	
	interval_sampled_model_points = M.shape[1]/approx_sampled_model_points
	M_sampled = M[:,0::interval_sampled_model_points] ### shape 3 x Nm_sampled  ####[0::N] selects every Nth point  
	
	interval_sampled_model_points_ICP = M.shape[1]/approx_sampled_model_points
	M_sampled_ICP = M[:,0::interval_sampled_model_points_ICP] ### shape 3 x Nm_sampled  ####[0::N] selects every Nth point  
	


	tree_M_sampled = KDTree(M_sampled.T)
	tree_M_sampled_ICP = KDTree(M_sampled_ICP.T)
	tree_M = KDTree(M.T)


	# Model
	m = Model("GurobiTest2")

	Ns = S.shape[1] #Number of sensor points
	Nm = M.shape[1]

	Nm_sampled = M_sampled.shape[1]
	num_partitions_SOS2 = 50 #for SOS2 constraint
	print "Nm_sampled :" ,Nm_sampled
	print "Ns_sampled :", Ns_sampled
	print "Nm :", Nm



	### Don't disturb the order of the variables 
	########
	alpha = m.addVars(3,Ns_sampled, lb=0, ub = GRB.INFINITY)
	phi =  m.addVars(1,Ns_sampled , lb=0, ub = GRB.INFINITY)
	T = m.addVars(3,1 , lb=-GRB.INFINITY, ub = GRB.INFINITY)
	R = m.addVars(3,3 , lb=-1, ub = 1, name='rotationMatrix')
	w = m.addVars(3,3 , lb=0, ub = 1)
	Cb_sampled = m.addVars(Ns_sampled, Nm_sampled, vtype=GRB.BINARY )
	lam = m.addVars(num_partitions_SOS2,3,3, lb = 0) #lambda
	q = np.linspace(-1,1,num_partitions_SOS2)
	BigM = 1000
	phiMax = 1000
	m.update()
	#Constraints 



	####set R value ### Can be used for debugging purposes
	# for i in range(3):
	# 	for j in range(3):
	# 		if i==j:
	# 			m.addConstr(R[i,j] ==1, name='rotConstr2u{},{}'.format(i,j))
	# 		else:
	# 			m.addConstr(R[i,j] ==0, name='rotConstr2u{},{}'.format(i,j))
	# ##set T value
	# for l in range(3):
	# 	m.addConstr(T[l,0] == 0)	

	# for ii in range(0,3):
	# 	for jj in range(0,3):
	# 		m.addConstr(R[ii,jj] ==CorrectTransf[ii,jj], name='rotConstr2u{},{}'.format(ii,jj))

	# for ii in range (0,3):
	# 	m.addConstr(T[ii,0] ==CorrectTransf[ii,3], name='rotConstr2u{},{}'.format(ii,jj))

	#############


	# ###Adding initial guess using heuristics

	R_start, T_start = icp_test(S[:,0:num_sampled_sens_points_ICP].T, M_sampled_ICP.T, tree_M_sampled_ICP)
	m.Params.StartNodeLimit = 2000000000-2
	m.Params.TimeLimit = TimeLimit
	m.Params.OutputFlag = 1
	m.Params.OptimalityTol = 1e-9
	m.Params.IntFeasTol = 1e-9
	


	if set_start_sol_switch == 1:
		Cb_sampled_start, lam_start, w_start, alpha_start, phi_start  = find_all_opt_variables(S[:,0:Ns_sampled].T, M_sampled.T, tree_M_sampled, R_start, T_start)
		for ii in range (3):
			T[ii,0].start = T_start[ii,0]
			for jj in range(3):
				R[ii,jj].start = R_start[ii,jj]

		# Cb_sampled.start = Cb_sampled_start
		for ii in range (Ns_sampled):
			for jj in range (Nm_sampled):
				Cb_sampled[ii,jj].start =Cb_sampled_start[ii,jj] 

		# lam.start = lam_start
		for kk in range(num_partitions_SOS2):
			for ii in range(3):
				for jj in range(3):
					lam[kk,ii,jj].start = lam_start[kk,ii,jj] 

		# w.start = w_start
		for ii in range(3):
			for jj in range(3):
				w[ii,jj].start = w_start[ii,jj] 

		# print alpha_start, "alpha_start "
		for ii in range(3):
			for jj in range(Ns_sampled):
				alpha[ii,jj].start = alpha_start[ii,jj] 

		# print phi_start, "phi_start "
		for jj in range(Ns_sampled):
			phi[0,jj].start = phi_start[0,jj] 


		print phi_start.sum()/Ns_sampled ,"phi_start.sum()/Ns_sampled "

 








#############################################################################
################## Rotation matrix constraints ##############################
#############################################################################

	for i in range(3):
		sum1 = 0
		for j in range(3):
			sum1 = sum1 + R[j,i]*R[j,i]
		m.addConstr((sum1 <= 1 ) , name='rotConstr{},{}'.format(i,j))


	for i in range(3):	
		for j in range(i+1,3):
			sum1 = 0
			for k in range (0,3):
				sum1 = sum1 + (R[k,i] + R[k,j])*(R[k,i] + R[k,j])
			m.addConstr((sum1 <= 2 ) , name='rotConstr1{},{}'.format(i,j))


	for i in range(3):
		
		for j in range(i+1,3):
			sum1 = 0
			for k in range (0,3):
				sum1 = sum1 + (R[k,i] - R[k,j])*(R[k,i] - R[k,j])
			m.addConstr((sum1 <= 2 ) , name='rotConstr2{},{}'.format(i,j))


	for tmp1 in [-1,1]:
		for tmp2 in [-1,1]:
			sum1 = 0
			for i in range(0,3):
				sum1 = sum1 + (R[i,0] + tmp1* R[i,1] + tmp2 * R[i,2])*(R[i,0] + tmp1* R[i,1] + tmp2 * R[i,2])
			m.addConstr((sum1 <= 3) , name='rotConstr3{}'.format(i))


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




	######### W constraints
	sum1 = 0
	for i in range(3):
		sum1 = 0
		for j in range(3):
			sum1 = sum1 + w[i,j]
		m.addConstr(( sum1 >= 1) , name = 'w' )


	#### R and W constraints
	for i in range (3):
		for j in range(3):
			sum1 = 0
			sum2 = 0
			sum3 = 0
			for k in range(num_partitions_SOS2):
				sum1 = sum1 + lam[k,i,j]
				sum2 = sum2 + lam[k,i,j]*q[k]
				sum3 = sum3 + lam[k,i,j]*q[k]*q[k]
			m.addConstr(sum1 == 1 )
			m.addConstr(R[i,j] == sum2)
			m.addConstr(w[i,j] == sum3)

 
		#########SOS2 constrain
	lam_list = [0]*num_partitions_SOS2
	weight_list =  [0]*num_partitions_SOS2
	for ii in range (3):
		for jj in range(3):
			for k in range(num_partitions_SOS2):	
				lam_list[k] = lam[k,ii,jj]
				weight_list[k] = k+1
			m.addSOS(GRB.SOS_TYPE2, lam_list, weight_list)
					
##############################################################################



	# summation Cij + oi (Correspondance matrix with outliers)
	for ii in range (Ns_sampled): #try later with quicksum1 
		sum1 = 0
		for jj in range (Nm_sampled):
			sum1 = sum1 + Cb_sampled[ii,jj]  
		m.addConstr((sum1 == 1 ) , name='Combination Matrix with Outliers')



	for ii in range (Ns_sampled): #try later with quicksum1 		  
		m.addConstr((phi[0,ii] >=  0 ) , name='Phi')



	for l in range (3):	 
		for i in range (Ns_sampled):
			sum1 = 0
			sum3 = 0
			for k in range (3):
				sum1 = sum1 + R[l,k]*S[k,i]
			for j in range (Nm_sampled):
				sum3 = sum3 + M_sampled[l,j]*Cb_sampled[i,j]
			m.addConstr((alpha[l,i] >= (sum1 + T[l,0] - sum3 ) ), name="alphaplus({},{})".format(l,i) )
			m.addConstr((alpha[l,i] >= -(sum1 + T[l,0] - sum3 )), name="alphaminus({},{})".format(l,i) )
		

			
			
	for i in range(Ns_sampled):
		sum1 = 0
		for l in range(3):
			sum1 = sum1 + alpha[l,i]
		m.addConstr(phi[0,i] >= sum1)
	m.setObjective((phi.sum()/Ns_sampled) , GRB.MINIMIZE)
 


	global SolRepeatCounter 	
	global objval2

	SolRepeatCounter = 0
	objval2 = 0.0
	def mycallback(m, where):
		global SolRepeatCounter 	
		global objval2
		if 	SolRepeatCounter < 100:
			temp_T_ICP = np.zeros((1,3))
			temp_R_ICP = np.zeros((3,3))
			R_total_after_ICP = np.zeros((3,3))
			T_total_after_ICP = np.zeros((1,3))

			# transformed_for_ICP = S[:,0:Ns_sampled]		
			

			if where == GRB.Callback.MIPNODE:

				solAtNode =  m.cbGetNodeRel(m._vars)
				solAtNode_After_ICP = solAtNode 
				
				l = 0
				for i in range(3*Ns_sampled + 1*Ns_sampled, 3*Ns_sampled + 1*Ns_sampled +3): ###Extracting R and T from current node
					temp_T_ICP[0][l] = solAtNode[i]
					temp_R_ICP[0][l] = solAtNode[i+3]
					temp_R_ICP[1][l] = solAtNode[i+6]
					temp_R_ICP[2][l] = solAtNode[i+9]
					l = l+1

				###### SVD for valid transformation 		
				U_svd, S_svd, Vt_svd = np.linalg.svd(temp_R_ICP)
				R_init_ICP = np.dot(Vt_svd.T, U_svd.T)
				
				if np.linalg.det(R_init_ICP)<0:
					Vt_svd [2,:] *= -1
					R_init_ICP = np.dot(Vt_svd.T, U_svd.T)
				
				T_init_ICP = temp_T_ICP

				transformed_for_ICP = np.matmul(R_init_ICP,S[:,0:num_sampled_sens_points_ICP]) + np.transpose(T_init_ICP) 
				
				t = time.time()

				
				 
				# R_after_ICP, T_after_ICP = icp_test( V.T,transformed_for_ICP.T, M.T, tree_M, M_sampled.T, tree_M_sampled,SigmaS,F,points_per_face,ICP_triangle_proj_switch_callback,ICP_or_GICP_switch_callback)
				R_after_ICP,T_after_ICP = icp_test(transformed_for_ICP.T, M_sampled_ICP.T, tree_M_sampled_ICP)
				elapsed = time.time() - t
				print "time for ICP :" ,elapsed
				

				#########Processing After ICP
				total_transf_after_ICP = np.matmul(R_T_to_4X4 ((R_after_ICP, T_after_ICP)), R_T_to_4X4((R_init_ICP,T_init_ICP.T)))
				
				
				R_total_after_ICP = total_transf_after_ICP[0:3,0:3]
				T_total_after_ICP = total_transf_after_ICP[0:3,3]
				
				Cb_sampled_after_ICP, lam_after_ICP, w_after_ICP, alpha_after_ICP, phi_after_ICP = find_all_opt_variables(S[:,0:Ns_sampled].T, M_sampled.T, tree_M_sampled, R_total_after_ICP, T_total_after_ICP.reshape(3,1))

				print "phi_after_ICP.sum()/Ns_sampled", phi_after_ICP.sum()/Ns_sampled  

				list_T_after_ICP, len_list_T_after_ICP = array_to_list(T_total_after_ICP)
				list_R_after_ICP, len_list_R_after_ICP = array_to_list(R_total_after_ICP)
				list_Cb_sampled_after_ICP, len_list_Cb_sampled_after_ICP = array_to_list(Cb_sampled_after_ICP)
				list_w_after_ICP, len_list_w_after_ICP = array_to_list(w_after_ICP)
				list_phi_after_ICP, len_list_phi_after_ICP = array_to_list(phi_after_ICP)
				list_alpha_after_ICP, len_list_alpha_after_ICP = array_to_list(alpha_after_ICP)
				list_lam_after_ICP, len_list_lam_after_ICP = array_to_list(lam_after_ICP)

				start_count = 4*Ns_sampled ## counting for translation vector elements # 3*Ns for Alpha, Ns for Phi
				m.cbSetSolution(m._vars[start_count : start_count + len_list_T_after_ICP ], list_T_after_ICP )
				start_count = 4*Ns_sampled + 3 
				m.cbSetSolution(m._vars[start_count : start_count + len_list_R_after_ICP ], list_R_after_ICP )
				start_count = 4*Ns_sampled +3 +9 +9 
				m.cbSetSolution(m._vars[start_count : start_count + len_list_Cb_sampled_after_ICP ], list_Cb_sampled_after_ICP )
				start_count = 0
				m.cbSetSolution(m._vars[start_count : start_count + len_list_alpha_after_ICP ], list_alpha_after_ICP )
				start_count = 3*Ns_sampled
				m.cbSetSolution(m._vars[start_count : start_count + len_list_phi_after_ICP ], list_phi_after_ICP)
				start_count = 4*Ns_sampled +3 +9 +9 + Ns_sampled*Nm_sampled
				 				
				objval = m.cbUseSolution()
				
				elapsed = time.time() - t

				if  objval == 1e+100 :							
					SolRepeatCounter = SolRepeatCounter + 1 

				else:
					SolRepeatCounter = 0
				print "objval", objval  
				print 'SolRepeatCounter: ', SolRepeatCounter
				print "................"



	m._vars = m.getVars()
	if callback_switch == 1:
		m.optimize(mycallback)
	else: 
		m.optimize()





	##############Number of solutions 
	number_of_solutions = 100 
	m.Params.PoolSolutions = number_of_solutions



	t1_ang_err = []
	t1_tra_err = []
	t1_tf = []
	t1_corr_python = []
	t1_corr_matlab = []
	t1_obj_val = m.objVal
	# wb = xlwt.Workbook()
	# ws = wb.add_sheet("bucket_error",cell_overwrite_ok=True)
	# ws.write(0,0,"Tejas1")
	# ws.write(1,0,"bucket number")
	# ws.write(1,1,"Error in deg")
	print "data = ", data
	for iii in range(number_of_solutions):
		Bucket = np.zeros((m.SolCount,4,4))
		if iii < m.SolCount:
			m.Params.SolutionNumber = iii 
			
			# Cb_sampled_sol = np.zeros((Ns_sampled,Nm_sampled))
			# for k in range (0,Ns_sampled):
			# 	for j in range(0,Nm_sampled):
			# 		Cb_sampled_sol[k,j]  =  Cb_sampled[k,j].Xn

			R_bucket_input = np.zeros((3,3))
			for k in range (0,3):
				for j in range (0,3):
					R_bucket_input[k,j] =  R[k,j].Xn

			T_bucket_input = np.zeros((3,1))
			for k in range(0,3):
				T_bucket_input[k,0] = T[k,0].Xn


			U_svd, S_svd, Vt_svd = np.linalg.svd(R_bucket_input)
			R_bucket_input_valid = np.dot(U_svd,Vt_svd)
		
			if np.linalg.det(R_bucket_input_valid)<0:
				Vt_svd [2,:] *= -1
				R_bucket_input_valid = np.dot(U_svd, Vt_svd)

			# Tejas1_corr_ids = np.zeros((Ns_sampled,1))
			# for i in range(Ns_sampled):
			# 	Tejas1_corr_ids[i,0] = np.where((Cb_sampled_sol[i,:] <= 1.1) & (Cb_sampled_sol[i,:] >= 0.9))[0][0] 
			
			
			transformed_sens_points = R_bucket_input_valid.dot(S)+T_bucket_input.reshape(3,1)
			_,corr_mat_ids = find_Cb (transformed_sens_points.T, Nm ,tree_M )


			T1_input_valid = R_T_to_4X4 ((R_bucket_input_valid,T_bucket_input))
			error_mat = T1_input_valid.dot(np.linalg.inv(gt))
			angle_error = abs(transforms3d.axangles.mat2axangle(error_mat[0:3,0:3])[1]*180/np.pi)
			print " bucket number = ", iii
			print np.linalg.norm(angle_error),"angle_error_norm"
			position_error = error_mat[0:3,3]


################################### bucket refinement for Tejas1 ################################
			# tf_after_bucket = bucket_refinement(R_bucket_input,T_bucket_input, S, M, tree_M )
			# error_mat = tf_after_bucket.dot(np.linalg.inv(gt))
			# angle_error = np.asarray(transforms3d.euler.mat2euler(error_mat))*180/np.pi
			# position_error = error_mat[0:3,3]
			# print np.linalg.norm(angle_error),"angle_error after bucket"
			# print ".."
################################### bucket refinement for Tejas1 ################################




			# if iii == 0:
 
			print corr_mat_ids[0:Ns_sampled] + 1, "Correspondance for matlab (1 alredy added)" 
			print R_bucket_input, "invalid R"
			print R_bucket_input_valid, "valid R"
			print T_bucket_input,"T_bucket_input"
			print ""
			print ""
			print ""
			
			t1_ang_err.append(np.linalg.norm(angle_error))
			t1_tra_err.append(np.linalg.norm(position_error))
			t1_tf.append(T1_input_valid)
			t1_corr_python.append(corr_mat_ids)
			t1_corr_matlab.append(corr_mat_ids+1)
###################### saving inputs for Tejas5
			# if iii == 0:	

			# 	path_name_save = "/home/biorobotics/Desktop/tejas/gurobiCodesPython/tejas5/faceNeglect/Tejas1_after_bucket/"  
			# 	filename_save = data + "_T1_corr.txt"
			# 	print path_name_save+filename_save, "path_name_save+filename_save"

			# 	np.savetxt(path_name_save+filename_save,Tejas1_corr_ids,fmt="%i")

			# 	filename_save = data + "_results"
			# 	np.savez(path_name_save+filename_save, R = tf_after_bucket[0:3,0:3], T = tf_after_bucket[0:3,3].reshape(3,1))
###################### saving inputs for Tejas5


		else:
			break	
	# wb = openpyxl.load_workbook("tejas1_results_21_jan.xlsx")
	# ws = wb.active
	# ws.append(t1_ang_err)

	# wb.save("tejas1_results_21_jan.xlsx")
	# print "yaha aaya?"
	# f.close()
	# f2.close()

	print "\n\n\n\n"
	return t1_ang_err,\
			t1_tra_err,\
			t1_tf,\
			t1_corr_python,\
			t1_corr_matlab,\
			t1_obj_val

