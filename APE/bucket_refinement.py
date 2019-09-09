import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def R_T_to_4X4 ((R,T)):
	T= np.reshape(T,(3,1))
	X = np.concatenate((R,T), axis = 1)
	X = np.concatenate((X,np.array([[0, 0, 0, 1]])), axis = 0)
	return X

def _4X4_to_R_T ((Transf)):
	R = Transf [0:3, 0:3]
	T = Transf [0:3, 3]
	return R,np.reshape(T,(3,1))


def bucket_refinement(R_bucket_input,T_bucket_input, S, M,tree_M ):	
	from icp_simple import icp_test

	transformed_sens_points_before_bucket = R_bucket_input.dot(S) + T_bucket_input.reshape(3,1)	
	R_after_ICP,T_after_ICP = icp_test(transformed_sens_points_before_bucket.T, M.T, tree_M)
	total_transf_after_bucket_refinement = (np.matmul(R_T_to_4X4((R_after_ICP, T_after_ICP)), R_T_to_4X4((R_bucket_input,T_bucket_input)))) 				
	transformed_after_bucket_refinement =  np.matmul(R_after_ICP,transformed_sens_points_before_bucket) + np.reshape(T_after_ICP,(3,1))

########### Plotting After bucket refinement
	
	plot_switch = 0
	if plot_switch ==1:
		fig = plt.figure()
		ax = fig.add_subplot(111, projection='3d')
		plt.title('After bucket refinement')

		ax.scatter(S[0,:],S[1,:],S[2,:], s = 20, c= 'blue' , label = 'S',alpha = 0.3)
		ax.scatter(transformed_after_bucket_refinement[0,:],transformed_after_bucket_refinement[1,:],transformed_after_bucket_refinement[2,:],s=20, c ='green', label = 'transformed_after_bucket_refinement', alpha = 0.3)
		ax.scatter(transformed_after_bucket_refinement[0,0:Ns_sampled],transformed_after_bucket_refinement[1,0:Ns_sampled],transformed_after_bucket_refinement[2,0:Ns_sampled],s=20, c ='black', label = 'transformed_sampled', alpha = 0.9)
		plt.legend()
		plt.show()

	return total_transf_after_bucket_refinement