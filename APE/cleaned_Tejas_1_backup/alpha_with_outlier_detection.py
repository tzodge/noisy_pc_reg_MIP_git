from gurobipy import *
import numpy as np 
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pickle
# from sklearn.neighbors import KDTree
from scipy.spatial import KDTree
import transforms3d #https://matthew-brett.github.io/transforms3d/reference/transforms3d.euler.html#eulerfuncs

from icp_with_outlier_handling import icp_test
#http://www.gurobi.com/documentation/8.0/refman/py_model_addconstrs.html i !=j , i!=-j+1 

#matrix functions



def find_Cb (transformed_sens_points, num_model_points ,tree_of_model ):
	### transformed_sens_points Nsx3
	### returns Cb matrix of dim [num(transformed_sens_points) X num(tree_of_model)]
	Cb = np.zeros((transformed_sens_points.shape[0], num_model_points))
	neighbor_idx = tree_of_model.query(transformed_sens_points)[1]
	neighbor_dist = tree_of_model.query(transformed_sens_points)[0]
	for ii in range (transformed_sens_points.shape[0]):
		if neighbor_dist[ii] <= phiMax:
			Cb[ii,neighbor_idx[ii]] = 1
	return Cb


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



def find_alpha(R,T,S,M,C,o):  
## dimensions --> S 3xn,
## dimensions --> T 3x1, 
## dimensions --> C S.shape[1]xM.shape[1], 
	alpha = np.zeros((3,S.shape[1]))
	T = np.reshape(T,(3,))
	for ii in range (S.shape[1]):
		# alpha[:,ii] = abs(B[ii].dot(S[:,ii] - T - R.dot(M).dot(C[ii,:].T)))
		# print C[ii,:].T.shape, "inside find_alpha"
		# print R.dot(S[:,ii] + T 
		alpha[:,ii] = abs(R.dot(S[:,ii]) + T - M.dot(C[ii,:].T)) -BigM*o[0,ii]
			
	return alpha


def find_phi(alpha):
	phi = alpha.sum(axis = 0)
	for ii in range(phi.shape[0]):
		if phi[ii] <= 0:
			phi[ii] = phiMax
	# print phi , "phi during callback"		
	return np.reshape(phi,(1,alpha.shape[1]))	

def find_o(C):	
	o = np.ones((1,C.shape[0])) - C.sum(axis = 1)
	print o, "Outliers"
	return  o



def find_all_opt_variables(sens_points, model_points, model_tree, R, T):
	
	transformed_sens_points = (np.matmul(R,sens_points.T)+T)
	# print transformed_sens_points.T.shape,"transformed_sens_points.T.shape......................"
	# print model_points.shape[0], "model_points.shape[0] ......"
	Cb_start = find_Cb(transformed_sens_points.T, model_points.shape[0], model_tree)
	o_start = find_o(Cb_start) 
	lam_start  = find_lam(R)
	w_start = find_w(lam_start)
	# print sens_points.T , "sens_points.T given to find_alpha"
	alpha_start = find_alpha(R,T,sens_points.T, model_points.T,Cb_start,o_start)
	phi_start = find_phi(alpha_start)
	# print phi_start.shape,  "phi_start.shape"

	return Cb_start, lam_start, w_start, alpha_start, phi_start,o_start 

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

def save2Ddata(f,var, name ):
	num_rows = var.shape[0]	
	num_cols = var.shape[1]
	f.write(name+"\n")
	f.write("{}\n".format(num_rows))
	for i in range(num_rows):
		for j in range(num_cols):
			if j == num_cols-1:	
				f.write("{}\n".format(var[i,j] ))
			else:
				f.write("{} ".format(var[i,j] ))
	
def save1Ddata(f,var):
	num_rows = var.shape[0]	
	num_cols = var.shape[1]
	for i in range(num_rows):
		for j in range(num_cols):
			f.write("{} ".format(var[i,j]))
	f.write("		")
			
def transf_mat_to_euler_6x1(transf_mat):
	angle_sxyz = np.asarray(transforms3d.euler.mat2euler(transf_mat))
	transl_xyz = np.array(transf_mat[0:3,3])
	return np.append(angle_sxyz,transl_xyz).reshape(1,6)


def main_func(model_data, \
	model_path, \
	sensor_data, \
	sensor_path, \
	set_start_sol_switch,\
	callback_switch,\
	ICP_or_GICP_switch_bucket,\
	ICP_triangle_proj_switch_bucket,\
	ICP_or_GICP_switch_callback,\
	ICP_triangle_proj_switch_callback,\
	Ns_sampled,\
	num_sampled_sens_points_ICP,\
	approx_sampled_model_points,\
	partial_sens_pnt_switch,\
	St_as_M_switch):
	
	global S,V,M,gt,B,SigmaS,F,num_partitions_SOS2,q, phiMax, BigM

	outl_fract = 0.05
	if partial_sens_pnt_switch == 1 :

		S_given = np.loadtxt( sensor_path+ sensor_data +'/S.txt',delimiter = ',') # S_given Ns_given X 3
		St_given = np.loadtxt( sensor_path+ sensor_data +'/St.txt',delimiter = ',') # St_given Ns_given X 3
		M_given = np.loadtxt(model_path + model_data +'/Msampled.txt',delimiter = ',') #M_given Nm_given X 3
		gt = np.loadtxt( sensor_path + sensor_data +'/gt.txt',delimiter = ',') # gt 4 X 4

		V_given = np.loadtxt(model_path + model_data +'/M.txt',delimiter = ',') ## Nv x 3
		B_given = np.loadtxt( sensor_path + sensor_data +'/B.txt',delimiter = ',') ## Ns_given x 6
		SigmaS_given = np.loadtxt( sensor_path + sensor_data +'/SigmaS.txt',delimiter = ',') ## Ns_given x 6
		F_given = np.loadtxt(model_path + model_data +'/F.txt',delimiter = ',') ## Nf x Nm


	else: 
		data = model_data
		path = model_path
		S_given = np.loadtxt( path+data +'/S.txt',delimiter = ',') # S_given Ns_given X 3
		St_given = np.loadtxt( path+data +'/St.txt',delimiter = ',') # St_given Ns_given X 3
		M_given = np.loadtxt(path +data +'/Msampled.txt',delimiter = ',') #M_given Nm_given X 3
		gt = np.loadtxt(path +data +'/gt.txt',delimiter = ',') # gt 4 X 4

		V_given = np.loadtxt(path +data +'/M.txt',delimiter = ',') ## Nv x 3
		B_given = np.loadtxt(path +data +'/B.txt',delimiter = ',') ## Ns_given x 6
		SigmaS_given = np.loadtxt(path +data +'/SigmaS.txt',delimiter = ',') ## Ns_given x 6
		F_given = np.loadtxt(path +data +'/F.txt',delimiter = ',') ## Nf x Nm

		

		
	
	SigmaS = SigmaS_given
	M = M_given.T
	F = F_given
	V = V_given.T  # V 3xNv	 Vertices	
	S = S_given.T
	St = St_given.T
	print S.shape, "S.shape"
	print M.shape, "M.shape"
	print F.shape, "F.shape"
	print V_given.shape, "V_given.shape"

###########  testing purposes
	
	M = St_given.T

############
	####  we have sampled from total M to save the testing time
	#### you can just set the divider to 1 instead of 300 and the code would take all the model points 
	interval_sampled_model_points = M.shape[1]/approx_sampled_model_points
	M_sampled = M[:,0::interval_sampled_model_points] ### shape 3 x Nm_sampled  ####[0::N] selects every Nth point  
	print M_sampled.shape, "M_sampled.shape"
	# M_icp = np.transpose(M_sampled)
	tree_M_sampled = KDTree(M_sampled.T)
	tree_M = KDTree(M.T)

	f = open("output.txt",'w+' )
	save2Ddata(f,gt,"gt")
	f.write("\n\n")

	f2 = open("output2.txt",'w+' )

	f2.write("gt		")
	save1Ddata(f2,transf_mat_to_euler_6x1(gt))
	f2.write("\n")

	# Model
	m = Model("GurobiTest2")

	Ns = S.shape[1] #Number of sensor points
	Nm = M.shape[1]
	Nf = F.shape[0]

	Nm_sampled = M_sampled.shape[1]
	num_partitions_SOS2 = 10 #for SOS2 constraint
	points_per_face = Nm/Nf

	num_outl_pts = int(outl_fract* Ns)
	
	print "num_outl_pts", num_outl_pts  
	print "points_per_face = Nm/Nf  ", points_per_face 
	print "Nm_sampled :" ,Nm_sampled
	print "Ns_sampled :", Ns_sampled

	###### Adding Outliers
	



	add_synthentic_outlier_switch = 0
	if add_synthentic_outlier_switch ==1: 	
		for i in range(num_outl_pts):
			a = np.random.randint(Ns)
			S[:,a] = (5*np.random.random((3,1)) - 2.5*np.ones((3,1))).reshape(3,)  
			print a,"", "a"


	plot_switch = 1
	if plot_switch == 1:
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		plt.title('before refinement')
		ax.scatter(M_sampled[0,:],M_sampled[1,:],M_sampled[2,:],s = 20, c= 'red',alpha = 0.5, label = 'M_sampled')
		ax.scatter(S[0,:],S[1,:],S[2,:], s = 20, c= 'blue', label ='S_given', alpha = 0.1)
		ax.scatter(S[0,0:Ns_sampled],S[1,0:Ns_sampled],S[2,0:Ns_sampled], s = 50, c= 'black', label ='S_sampled', alpha = 1)
		plt.legend()
		plt.show()

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
	o = m.addVars(1,Ns , vtype=GRB.BINARY )
	BigM = 1000
	phiMax = 0.15
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

	R_start, T_start = icp_test(V.T,S[:,0:num_sampled_sens_points_ICP].T, M.T, tree_M, M_sampled.T, tree_M_sampled,SigmaS,F,points_per_face,ICP_triangle_proj_switch_bucket,ICP_or_GICP_switch_bucket)
	# print S[:,0:num_sampled_sens_points_ICP].T,"S[:,0:num_sampled_sens_points_ICP].T"
	m.Params.StartNodeLimit = 2000000000-2
	m.Params.TimeLimit = 10*60
	m.Params.OutputFlag = 1
	# m.Params.MIPFocus=3 
	Cb_sampled_start, lam_start, w_start, alpha_start, phi_start, o_start  = find_all_opt_variables(S[:,0:Ns_sampled].T, M_sampled.T, tree_M_sampled, R_start, T_start)

##### Providing a warm start 
	if set_start_sol_switch == 1:
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

	

#####################################################################
##Constraints
####################################################################3

	# Ns_sampled
	# #Constraints on rotation matrix

	# sum1mation R(i)^2 <= 1
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



	######### W constraints
	sum1 = 0
	for i in range(3):
		sum1 = 0
		for j in range(3):
			sum1 = sum1 + w[i,j]
		m.addConstr(( sum1 >= 1) , name = 'w' )


	###SOS constraints  http://www.gurobi.com/documentation/8.0/refman/constraints.html
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


	#SOS2 constrain
	for i in range (3):
		for j in range(3):
			# m.addSOS(GRB.SOS_TYPE2, [lam[i,j,0],lam[i,j,1],lam[i,j,2],lam[i,j,3],lam[i,j,4],lam[i,j,5]] , [1,2,3,4,5,6])
			m.addSOS(GRB.SOS_TYPE2, [lam[0,i,j],lam[1,i,j],lam[2,i,j],lam[3,i,j],lam[4,i,j],lam[5,i,j],lam[6,i,j],lam[7,i,j],lam[8,i,j],lam[9,i,j]] , [1,2,3,4,5,6,7,8,9,10])
			



	# summation Cij + oi (Correspondance matrix with outliers)
	for ii in range (Ns_sampled): #try later with quicksum1 
		sum1 = 0
		for jj in range (Nm_sampled):
			sum1 = sum1 + Cb_sampled[ii,jj]  
		m.addConstr((sum1 + o[0,ii] == 1 ) , name='Combination Matrix with Outliers')
		# m.addConstr((sum1 == 1 ) , name='Combination Matrix with Outliers')



	for ii in range (Ns_sampled): #try later with quicksum1 		  
		m.addConstr((phi[0,ii] >=  phiMax *o[0,ii] ) , name='Phi')
		# m.addConstr((phi[0,ii] >=  0 ) , name='Phi')



	for l in range (3):	 
		for i in range (Ns_sampled):
			sum1 = 0
			sum3 = 0
			for k in range (3):
				sum1 = sum1 + R[l,k]*S[k,i]
			for j in range (Nm_sampled):
				sum3 = sum3 + M_sampled[l,j]*Cb_sampled[i,j]
			m.addConstr((alpha[l,i] >= (sum1 + T[l,0] - sum3 )- BigM*o[0,i]    ), name="fdf" )
			m.addConstr((alpha[l,i] >= -(sum1 + T[l,0] - sum3 )- BigM*o[0,i]  ), name="fdf" )
			# m.addConstr((alpha[l,i] >= (sum1 + T[l,0] - sum3 ) ), name="alphaplus({},{})".format(l,i) )
			# m.addConstr((alpha[l,i] >= -(sum1 + T[l,0] - sum3 )), name="alphaminus({},{})".format(l,i) )
		

			
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
				start_count = 4*Ns_sampled +3 +9 +9 +Ns_sampled*Nm_sampled +num_partitions_SOS2*3*3 
				o_before_ICP = solAtNode[start_count:start_count+Ns_sampled]
				

				l = 0
				for i in range(3*Ns_sampled + 1*Ns_sampled, 3*Ns_sampled + 1*Ns_sampled +3): ###Extracting R and T from current node
					temp_T_ICP[0][l] = solAtNode[i]
					temp_R_ICP[0][l] = solAtNode[i+3]
					temp_R_ICP[1][l] = solAtNode[i+6]
					temp_R_ICP[2][l] = solAtNode[i+9]
					l = l+1

				# print temp_R_ICP, "temp_R_ICP"
				# print temp_T_ICP, "temp_T_ICP"
				
				###### SVD for valid transformation 		
				U_svd, S_svd, Vt_svd = np.linalg.svd(temp_R_ICP)
				R_init_ICP = np.dot(Vt_svd.T, U_svd.T)
				
				if np.linalg.det(R_init_ICP)<0:
					Vt_svd [2,:] *= -1
					R_init_ICP = np.dot(Vt_svd.T, U_svd.T)
				
				T_init_ICP = temp_T_ICP

				transformed_for_ICP = np.matmul(R_init_ICP,S[:,0:num_sampled_sens_points_ICP]) + np.transpose(T_init_ICP) 
				
				t = time.time()

				
				 
				R_after_ICP, T_after_ICP = icp_test( V.T,transformed_for_ICP.T, M.T, tree_M, M_sampled.T, tree_M_sampled,SigmaS,F,points_per_face,ICP_triangle_proj_switch_callback,ICP_or_GICP_switch_callback)
				elapsed = time.time() - t
				print "time for ICP :" ,elapsed
				

				#########Processing After ICP
				# Transformed_After_ICP = np.matmul(R_after_ICP,transformed_for_ICP[:,0:Ns]) + T_after_ICP 
				total_transf_after_ICP = np.matmul(R_T_to_4X4 ((R_after_ICP, T_after_ICP)), R_T_to_4X4((R_init_ICP,T_init_ICP.T)))
				
				# print errorAfterICP((total_transf_after_ICP[0:3,0:3], total_transf_after_ICP[0:3,3])) , "errorAfterICP"
				
				R_total_after_ICP = total_transf_after_ICP[0:3,0:3]
				T_total_after_ICP = total_transf_after_ICP[0:3,3]
				
				Cb_sampled_after_ICP, lam_after_ICP, w_after_ICP, alpha_after_ICP, phi_after_ICP, o_after_ICP = find_all_opt_variables(S[:,0:Ns_sampled].T, M_sampled.T, tree_M_sampled, R_total_after_ICP, T_total_after_ICP.reshape(3,1))

				print "phi_after_ICP.sum()/Ns_sampled", phi_after_ICP.sum()/Ns_sampled  

				list_T_after_ICP, len_list_T_after_ICP = array_to_list(T_total_after_ICP)
				list_R_after_ICP, len_list_R_after_ICP = array_to_list(R_total_after_ICP)
				list_Cb_sampled_after_ICP, len_list_Cb_sampled_after_ICP = array_to_list(Cb_sampled_after_ICP)
				list_w_after_ICP, len_list_w_after_ICP = array_to_list(w_after_ICP)
				list_phi_after_ICP, len_list_phi_after_ICP = array_to_list(phi_after_ICP)
				list_alpha_after_ICP, len_list_alpha_after_ICP = array_to_list(alpha_after_ICP)
				list_lam_after_ICP, len_list_lam_after_ICP = array_to_list(lam_after_ICP)
				list_o_after_ICP, len_list_o_after_ICP = array_to_list (o_after_ICP) 
				# print list_o_after_ICP, "list_o_after_ICP"

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
				start_count = 4*Ns_sampled +3 +9 +9 +Ns_sampled*Nm_sampled +num_partitions_SOS2*3*3 
				m.cbSetSolution(m._vars[start_count : start_count + len_list_o_after_ICP ], list_o_after_ICP)

				# start_count = 4*Ns_sampled +3 +9 +9 + Ns_sampled*Nm_sampled
				# m.cbSetSolution(m._vars[start_count : start_count + len_list_lam_after_ICP], list_lam_after_ICP)


				# m.cbSetSolution(m._vars[4*Ns_sampled:4*Ns_sampled+12],solAtNode_After_ICP[4*Ns_sampled:4*Ns_sampled+12] )
				
				objval = m.cbUseSolution()
				
				elapsed = time.time() - t
				# ############# Termination Criteria

				if  objval == 1e+100 :							
					SolRepeatCounter = SolRepeatCounter + 1 

				else:
					SolRepeatCounter = 0
				print "objval", objval  
				print 'SolRepeatCounter: ', SolRepeatCounter
				print "................"




					
			

			######plotting	
			# TransformedAfterICPAll = np.matmul(R_after_ICP,transformed_fCP_All) + T_after_ICP
			
			# fig = plt.figure()
			# ax = fig.add_subplot(111, projection='3d')

			# ax.scatter(M_sampled[0,:],M_sampled[1,:],M_sampled[2,:], s = 5, c= 'red')
			# ax.scatter(transformed_for_ICP[0,:],transformed_for_ICP[1,:],transformed_for_ICP[2,:], s = 50, c= 'orange')
			# ax.scatter(TransformedAfterICPAll[0,:],TransformedAfterICPAll[1,:],TransformedAfterICPAll[2,:], s = 15, c= 'blue')
			# ax.scatter(S[0,:],S[1,:],S[2,:], s = 5, c = 'green')
			# ax.scatter(St[0,:],St[1,:],St[2,:], s = 5, c = 'green')
			
			# plt.show()





	### setting minimum objective value, setting parameters
	# S[:,0:Ns_sampled].T, M_sampled.T, tree_M_sampled, R_total_after_ICP, T_total_after_ICP.reshape(3,1)
	# _, _, _, _, phi_minimum  = find_all_opt_variables(S.T, M.T,  KDTree(M.T), gt[0:3,0:3], gt[0:3,3].reshape(3,1))
	# print  phi_minimum.sum()/Ns, " phi_minimum.sum()/Ns-------------------------------------"
	# m.Params.BestObjStop = phi_minimum.sum()/Ns
	# m.Params.TimeLimit = 500



	m._vars = m.getVars()
	if callback_switch == 1:
		m.optimize(mycallback)
	else: 
		m.optimize()




	R_bucket_input = np.zeros((3,3))
	for ii in range (0,3):
		for jj in range (0,3):
			R_bucket_input[ii,jj] =  R[ii,jj].X
	# print R_bucket_input
	# print ""

	T_bucket_input = np.zeros((3,1))
	for ii in range(0,3):
		T_bucket_input[ii,0] = T[ii,0].X

	# print T_bucket_input


	# o_sol = np.zeros((1,Ns))
	# for ii in range(0,Ns):
	# 	o_sol[0,ii] = o[0,ii].X

	# print "DDDDDDDDDDDDDDDDDDDDD"
	# print o_sol 

	# modified =  np.matmul(R_bucket_input,S) + T_bucket_input


	# fig = plt.figure()
	# ax = fig.add_subplot(111, projection='3d')

	# ax.scatter(M_sampled[0,:],M_sampled[1,:],M_sampled[2,:],s = 50, c= 'red')
	# ax.scatter(S[0,:],S[1,:],S[2,:], s = 10, c= 'blue')
	# ax.scatter(modified[0,:],modified[1,:],modified[2,:],s=10, c ='green')

	# plt.show()


	##############Number of solutions 
	number_of_solutions = 100 
	m.Params.PoolSolutions = number_of_solutions



	FivePointEstimation = np.empty((number_of_solutions,4))

	for iii in range(number_of_solutions):
		Bucket = np.zeros((m.SolCount,4,4))
		if iii < m.SolCount:
			m.Params.SolutionNumber = iii 
			
			Cb_sampled_sol = np.zeros((Ns_sampled,Nm_sampled))
			for k in range (0,Ns_sampled):
				for j in range(0,Nm_sampled):
					Cb_sampled_sol[k,j]  =  Cb_sampled[k,j].Xn

			o_sampled_sol = np.zeros((Ns_sampled,1))
			for k in range(0,Ns_sampled):
				o_sampled_sol[k,0] = o[0,k].Xn

			phi_sampled_sol = np.zeros((Ns_sampled,1))
			for k in range(0,Ns_sampled):
				phi_sampled_sol[k,0] = phi[0,k].Xn


			R_bucket_input = np.zeros((3,3))
			for k in range (0,3):
				for j in range (0,3):
					R_bucket_input[k,j] =  R[k,j].Xn

			T_bucket_input = np.zeros((3,1))
			for k in range(0,3):
				T_bucket_input[k,0] = T[k,0].Xn


			print o_sampled_sol.T, "o_sampled_sol "
			# print phi_sampled_sol.T, " phi_sampled_sol...." 

			U_svd, S_svd, Vt_svd = np.linalg.svd(R_bucket_input)
			R_bucket_input_valid = np.dot(U_svd,Vt_svd)
		
			if np.linalg.det(R_bucket_input_valid)<0:
				Vt_svd [2,:] *= -1
				R_bucket_input_valid = np.dot(U_svd, Vt_svd)


			#################################	
			error_mat = np.matmul(np.linalg.inv(R_T_to_4X4 ((R_bucket_input,T_bucket_input))), gt)
			angle_error = np.asarray(transforms3d.euler.mat2euler(error_mat))*180/np.pi
			position_error = error_mat[0:3,3]
			print "before bucket refinement"
			norm_angle_error_before = np.linalg.norm(angle_error)
			norm_position_error_before =  np.linalg.norm(position_error)	
			print "angle_error", norm_angle_error_before
			print "position_error", norm_position_error_before
			print "."
			modified =  np.matmul(R_bucket_input_valid,S) + T_bucket_input
			

			plot_switch = 1
			if plot_switch == 1:
				fig = plt.figure()
				ax = fig.add_subplot(111, projection='3d')
				plt.title('before refinement')
				ax.scatter(M_sampled[0,:],M_sampled[1,:],M_sampled[2,:],s = 10, c= 'red', label = 'M_sampled')
				ax.scatter(S[0,:],S[1,:],S[2,:], s = 010, c= 'blue', label ='S_given')
				ax.scatter(modified[0,:],modified[1,:],modified[2,:],s=010, c ='green', label = 'modified')
				# ax.scatter(S[0,0:Ns_sampled],S[1,0:Ns_sampled],S[2,0:Ns_sampled],s=100, c ='yellow', label = 'modified')
				
				plt.legend()
				plt.show()
			print " .."

			########################################
			bucket_refinement_switch = 0
			if bucket_refinement_switch == 1:	
				import bucket_refinement
		        
				total_transf_after_bucket_refinement = bucket_refinement.bucket_refinement(R_bucket_input_valid,T_bucket_input,Ns_sampled,V,S,M_sampled,tree_M_sampled,M,tree_M,gt,num_sampled_sens_points_ICP, SigmaS, F, points_per_face,ICP_or_GICP_switch_bucket,ICP_triangle_proj_switch_bucket)
				Bucket[iii] = total_transf_after_bucket_refinement

				error_mat = np.matmul(np.linalg.inv(R_T_to_4X4 ((total_transf_after_bucket_refinement[0:3,0:3],total_transf_after_bucket_refinement[0:3,3]))), gt)
				angle_error = np.asarray(transforms3d.euler.mat2euler(error_mat))*180/np.pi
				position_error = error_mat[0:3,3]
				print "after bucket refinement"
				norm_angle_error = np.linalg.norm(angle_error)
				norm_position_error =  np.linalg.norm(position_error) 	
				print "angle_error", norm_angle_error
				print "position_error", norm_position_error
				print ".."
				print ""


				save2Ddata(f,Bucket[iii],"Bucket[{}]".format(iii))
				f.write("Error in orientation before bucket refinement : " +"{}".format(norm_angle_error_before) + "\n")		
				f.write("Error in orientation after bucket refinement : " +"{}".format(norm_angle_error) + "\n")		
				f.write("\n\n")

				f2.write("{}		".format(iii))
				save1Ddata(f2,transf_mat_to_euler_6x1(Bucket[iii]))
				f2.write("{}		".format(norm_position_error_before))
				f2.write("{}		".format(norm_position_error))


				# print  Bucket[iii][0:3,0:3], "Bucket[iii][0:3,0:3]", "......"
	# find_all_opt_variables(sens_points, model_points, model_tree, R, T):
				
				_, _, _, _, phi_iii,_  = find_all_opt_variables(S.T,M.T,tree_M,Bucket[iii][0:3,0:3],Bucket[iii][0:3,3].reshape(3,1))
				f2.write("{}\n".format(phi_iii.sum()/Ns))


		else:
			break	

	f.close()
	f2.close()

	print "\n\n\n\n"
