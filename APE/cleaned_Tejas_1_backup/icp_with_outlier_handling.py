""" ICP algorithm

    References:
    (ICP)
    [1] Paul J. Besl and Neil D. McKay,
        "A method for registration of 3-D shapes",
        PAMI Vol. 14, Issue 2, pp. 239-256, 1992.
    (SVD)
    [2] K. S. Arun, T. S. Huang and S. D. Blostein,
        "Least-Squares Fitting of Two 3-D Point Sets",
        PAMI Vol. 9, Issue 5, pp.698--700, 1987
"""
import numpy as np
from scipy.spatial import KDTree
from GICP_gaussNewton import GICP_gaussnewton
# import  NearestNeighbour as NN 
import math
import time
import transforms3d

####switches
# ICP_or_GICP_switch_1 = 1
    ######### ICP_or_GICP_switch_1 = 1 for ICP 
    ######### ICP_or_GICP_switch_1 = 1 for GICP
# switch_for_triang_proj = 1
    ######### flag = 1 
    ######### flag = 2 


def _icp_find_rigid_transform(p_from, p_target, SigmaS,flag):
    ######### flag = 1 for SVD 
    ######### flag = 2 for Gaussian
    A, B = np.copy(p_from), np.copy(p_target)

    if flag == 1:

        centroid_A = np.mean(A, axis=0)
        centroid_B = np.mean(B, axis=0)

        A -= centroid_A
        B -= centroid_B

        H = np.dot(A.T, B)
        U, S, Vt = np.linalg.svd(H)
        R = np.dot(Vt.T, U.T)

        # special reflection case
        if np.linalg.det(R) < 0:
            Vt[2,:] *= -1
            R = np.dot(Vt.T, U.T)

        t = np.dot(-R, centroid_A) + centroid_B

    elif flag ==2:

        R,t1 = GICP_gaussnewton(A,B,SigmaS)
        t = t1[0]
    return R, t


def _icp_Rt_to_matrix(R, t):
    # matrix M = [R, t; 0, 1]
    Rt = np.concatenate((R, np.expand_dims(t.T, axis=-1)), axis=1)
    a = np.concatenate((np.zeros_like(t), np.ones(1)))
    M = np.concatenate((Rt, np.expand_dims(a, axis=0)), axis=0)
    return M

class ICP_triangle_proj:
    """ Estimate a rigid-body transform g such that:
        p0 = g.p1
    """
    def __init__(self, p0, p1,tree_M, SigmaS,F,points_per_face):
        self.p0 = p0
        self.p1 = p1
        self.nearest = tree_M
        self.g_series = None
        self.F = F
        self.SigmaS = SigmaS
        self.points_per_face = points_per_face
    def compute(self, max_iter):
        ftol = 1.0e-7
        dim_k = self.p0.shape[1]
        g = np.eye(dim_k + 1, dtype=self.p0.dtype)
        p = np.copy(self.p1)

        Nm = self.F.shape[1]
        # print Nm
        # = p1.shape[0]

        self.g_series = np.zeros((max_iter + 1, dim_k + 1, dim_k + 1), dtype=g.dtype)
        self.g_series[0, :, :] = g
        
        p_all_sens = np.copy(self.p1)
        num_sens_pts_without_outl = int(p_all_sens.shape[0]*(1-outl_fract_icp)) 
        num_outl = int(p_all_sens.shape[0]*outl_fract_icp) 

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        
        targets=np.zeros((num_sens_pts_without_outl,3))
        itr = -1
        # print "inside compute"
        # self.nearest.query(p[0,:])[1]
        import  NearestNeighbour as NN
        for itr in range(max_iter):

            # fig = plt.figure()
            # ax = Axes3D(fig)
            # ax = fig.add_subplot(111, projection='3d')

            # neighbor_idx = self.nearest.query(p)[1]
            search_result = self.nearest.query(p_all_sens)
            neighbor_idx = search_result[1]
            neighbor_dist = search_result[0]
            sorted_ind = np.argsort(neighbor_dist)
            p = np.delete(p_all_sens, sorted_ind[-num_outl:],axis = 0)
            neighbor_idx = np.delete(neighbor_idx, sorted_ind[-num_outl:],axis = 0)
            

            for ii in range(p.shape[0]):
                
                ind = [0,0,0]
                tt=0;
                for kk in range(Nm):
                    if self.F[int(float(neighbor_idx[ii])/self.points_per_face),kk] == 1:
                        
                        ind[tt] = kk
                        tt = tt+1
                    if tt == 3:
                        break
    

                targets[ii,:]= NN.ProjectedOnTriangle(p[ii,:], self.p0[ind[0],:], self.p0[ind[1],:], self.p0[ind[2],:] )[0]
                
                A0 = np.linalg.norm( np.cross(self.p0[ind[2],:]-self.p0[ind[0],:] , self.p0[ind[1],:]-self.p0[ind[0],:] ) )
                A1 = np.linalg.norm( np.cross(self.p0[ind[0],:]-targets[ii,:] , self.p0[ind[1],:]-targets[ii,:] ) ) 
                A2 = np.linalg.norm( np.cross(self.p0[ind[0],:]-targets[ii,:] , self.p0[ind[2],:]-targets[ii,:] ) ) 
                A3 = np.linalg.norm( np.cross(self.p0[ind[1],:]-targets[ii,:] , self.p0[ind[2],:]-targets[ii,:] ) ) 

                if abs(A0-(A1-A2-A3)) >= 0.001:
                    targets[ii,:] = (self.p0[ind[0],:]+ self.p0[ind[1],:]+ self.p0[ind[2],:])/3


                # ax.scatter(targets[ii,0], targets[ii,1], targets[ii,2],"o", color="black") 

            R, t = _icp_find_rigid_transform(p, targets, self.SigmaS, ICP_or_GICP_switch_1)

            
            new_p = np.dot(R, p.T).T + t
            # p_all_sens = np.copy(new_p)
            # new_p = np.dot(R, p.T).T + t
            error = np.sum(np.abs(p - new_p)) 
            if np.sum(np.abs(p - new_p)) < ftol:
                break


            new_p = np.dot(R, p_all_sens.T).T + t
            p_all_sens = np.copy(new_p)                
            dg = _icp_Rt_to_matrix(R, t)
            new_g = np.dot(dg, g)
            g = np.copy(new_g)

            self.g_series[itr + 1, :, :] = g
            
            # ax.set_label("x - axis")
            # ax.set_label("y - axis")
            # ax.set_label("z - axis")
            # ax.plot(self.p0[:,0], self.p0[:,1], self.p0[:,2], "o", color="red", ms=4, mew=0.5)
            # ax.plot(p[:,0], p[:,1], p[:,2], "o", color="blue", ms=8, mew=0)   
            # # ax.plot(targets[:,0],targets[:,1],targets[:,2],"o",color="green", ms=4, mew=0)
            # plt.show()        


        self.g_series[(itr+1):, :, :] = g

        return g, p, (itr + 1), error, neighbor_idx



class ICP:
    """ Estimate a rigid-body transform g such that:
        p0 = g.p1
    """
    def __init__(self, p0, p1,tree_M_sampled,SigmaS):
        self.p0 = p0
        self.p1 = p1
        # self.nearest = KDTree_M_sampled(self.p0)
        self.nearest = tree_M_sampled
        self.SigmaS = SigmaS

        self.g_series = None

    def compute(self, max_iter):
        ftol = 1.0e-7
        dim_k = self.p0.shape[1]
        g = np.eye(dim_k + 1, dtype=self.p0.dtype)
        p_all_sens = np.copy(self.p1)
        num_sens_pts_without_outl = int(p_all_sens.shape[0]*(1-outl_fract_icp)) 
        num_outl = int(p_all_sens.shape[0]*outl_fract_icp) 
        self.g_series = np.zeros((max_iter + 1, dim_k + 1, dim_k + 1), dtype=g.dtype)
        self.g_series[0, :, :] = g

        itr = -1
        for itr in range(max_iter):
            search_result = self.nearest.query(p_all_sens)
            neighbor_idx = search_result[1]
            neighbor_dist = search_result[0]
            sorted_ind = np.argsort(neighbor_dist)
            p = np.delete(p_all_sens, sorted_ind[-num_outl:],axis = 0)
            neighbor_idx = np.delete(neighbor_idx, sorted_ind[-num_outl:],axis = 0)
            
            targets = self.p0[neighbor_idx]
            R, t = _icp_find_rigid_transform(p, targets, self.SigmaS, ICP_or_GICP_switch_1)
           
            new_p = np.dot(R, p.T).T + t
     
            error = np.sum(np.abs(p - new_p)) 
            if np.sum(np.abs(p - new_p)) < ftol:
                break
            new_p = np.dot(R, p_all_sens.T).T + t
            
            p_all_sens = np.copy(new_p)
            dg = _icp_Rt_to_matrix(R, t)
            new_g = np.dot(dg, g)
            g = np.copy(new_g)
            self.g_series[itr + 1, :, :] = g


            show_plot_switch = 00 
            if show_plot_switch == 1:                       
                # print targets.shape, "targets.shape"  
                import matplotlib.pyplot as plt
                from mpl_toolkits.mplot3d import Axes3D

                fig = plt.figure()
                ax = Axes3D(fig)
                ax = fig.add_subplot(111, projection='3d')
                ax.set_label("x - axis")
                ax.set_label("y - axis")
                ax.set_label("z - axis")
                ax.plot(self.p0[:,0], self.p0[:,1], self.p0[:,2], "o", color="red", ms=4, mew=0.5)
                ax.plot(p_all_sens[:,0], p_all_sens[:,1], p_all_sens[:,2], "o", color="blue", ms=8, mew=0)   
                ax.plot(targets[:,0],targets[:,1],targets[:,2],"o",color="green", ms=4, mew=0)
                plt.show()     
                # time.sleep(1)   


        self.g_series[(itr+1):, :, :] = g


        return g, p, (itr + 1), error, neighbor_idx



#### please note that the S over here, is not total niumber of sensor points
#### it is the number  sensor points sampled for ICP 
def icp_test(V,S,M,tree_M, M_sampled, tree_M_sampled,SigmaS,F,points_per_face,ICP_triangle_proj_switch,ICP_or_GICP_switch):
    #### V = 3 x Nv
    #### S = 3 x sampled_sensor_points_for_icp
    #### F = Nf x Nm

    #### V = [R,t;0,1] S
    #### gt = [matrix 4x4] S

    from math import sin, cos
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    import NearestNeighbour
    global ICP_or_GICP_switch_1, outl_fract_icp 
    

    show_plot_switch = 1
    max_iter = 30


    outl_fract_icp = 0.1    

    ICP_or_GICP_switch_1 = ICP_or_GICP_switch

    if ICP_triangle_proj_switch == 1:
        icp = ICP_triangle_proj(V, S, tree_M, SigmaS, F, points_per_face)
    else: 
        icp = ICP(M_sampled, S, tree_M_sampled,SigmaS)
    matrix, points, itr,error, neighbor_idx = icp.compute(max_iter)



    show_plot_switch = 0
    error_calc_switch = 0
      
    if show_plot_switch == 1:
    ##### plotting after ICP_triangle_proj
        # print "elapsed time ", time.time() -t0
        rot = matrix[0:3,0:3]
        tra = np.reshape(np.transpose(matrix[0:3,3]),(3,1))
        transformed = np.matmul(rot, np.transpose(S)) + tra    

        if error_calc_switch ==1:
            error_mat =  matrix.dot(np.linalg.inv(gt_given))
            angle_error = np.asarray(transforms3d.euler.mat2euler(error_mat))*180/np.pi
            print "orientation before ICP", np.linalg.norm(np.asarray(transforms3d.euler.mat2euler(gt_given))*180/np.pi)
            print "orientation after ICP",np.linalg.norm(angle_error)
        fig = plt.figure()
        ax = Axes3D(fig)
        ax = fig.add_subplot(111, projection='3d')
        ax.set_label("x - axis")
        ax.set_label("y - axis")
        ax.set_label("z - axis")
        ax.plot(V[:,0], V[:,1], V[:,2], "o", color="red", ms=4, mew=0.5)
        ax.plot(S[:,0], S[:,1], S[:,2], "o", color="blue", ms=1, mew=0)   
        ax.plot(transformed[0,:],transformed[1,:],transformed[2,:],"o",color="green", ms=4, mew=0)
        plt.show()        



    print itr, "itr" 
    return matrix[0:3,0:3] , np.reshape(matrix[0:3,3],(3,1))  





# if __name__ == '__main__':
#     # print "aa gaya"
#     # icp_test(A,B,tree_M,SigmaS,F)
#     icp_test(V,S,tree_M,SigmaS,F,points_per_face)






##########EOF

######for testing purposes


# data = 'dragon_20_deg'
# path = '/home/biorobotics/Desktop/tejas/Downloads/tempdatasets/noise_free'+ '/'

# ####### switches
# outlier_switch = 0
# M_sampled_switch_global_optimiser = 1 ########### if 1, uses M_sampled in global optimizer 


# # # ########

# S_given = np.loadtxt(path +data +'/S.txt',delimiter = ',')
# V_given = np.loadtxt(path +data +'/M.txt',delimiter = ',')
# M_given = np.loadtxt(path +data +'/Msampled.txt',delimiter = ',')
# F_given = np.loadtxt(path +data +'/F.txt',delimiter = ',')
# gt_given = np.loadtxt(path +data +'/gt.txt',delimiter = ',')
# B_given = np.loadtxt(path +data +'/B.txt',delimiter = ',')
# SigmaS_given = np.loadtxt(path +data +'/SigmaS.txt',delimiter = ',')

# print F_given.shape, " F_given.shape" 

# V = V_given.T
# S = S_given[0:500,:].T


# model_data = 'bunny_partial_view_1'
# model_path = '/home/biorobotics/Desktop/tejas/gurobiCodesPython/datasets'+ '/'

# sensor_data = 'bunny_partial_view_1'
# sensor_path = '/home/biorobotics/Desktop/tejas/gurobiCodesPython/datasets'+ '/'


# partial_sens_pnt_switch = 1
# if partial_sens_pnt_switch == 1 :

#     S_given = np.loadtxt( sensor_path+ sensor_data +'/S.txt',delimiter = ',') # S_given Ns_given X 3
#     St_given = np.loadtxt( sensor_path+ sensor_data +'/St.txt',delimiter = ',') # St_given Ns_given X 3
#     M_given = np.loadtxt(model_path + model_data +'/Msampled.txt',delimiter = ',') #M_given Nm_given X 3
#     gt = np.loadtxt( sensor_path + sensor_data +'/gt.txt',delimiter = ',') # gt 4 X 4

#     V_given = np.loadtxt(model_path + model_data +'/M.txt',delimiter = ',') ## Nv x 3
#     B_given = np.loadtxt( sensor_path + sensor_data +'/B.txt',delimiter = ',') ## Ns_given x 6
#     SigmaS_given = np.loadtxt( sensor_path + sensor_data +'/SigmaS.txt',delimiter = ',') ## Ns_given x 6
#     F_given = np.loadtxt(model_path + model_data +'/F.txt',delimiter = ',') ## Nf x Nm


# else: 
#     data = model_data
#     path = model_path
#     S_given = np.loadtxt( path+data +'/S.txt',delimiter = ',') # S_given Ns_given X 3
#     St_given = np.loadtxt( path+data +'/St.txt',delimiter = ',') # St_given Ns_given X 3
#     M_given = np.loadtxt(path +data +'/Msampled.txt',delimiter = ',') #M_given Nm_given X 3
#     gt = np.loadtxt(path +data +'/gt.txt',delimiter = ',') # gt 4 X 4

#     V_given = np.loadtxt(path +data +'/M.txt',delimiter = ',') ## Nv x 3
#     B_given = np.loadtxt(path +data +'/B.txt',delimiter = ',') ## Ns_given x 6
#     SigmaS_given = np.loadtxt(path +data +'/SigmaS.txt',delimiter = ',') ## Ns_given x 6
#     F_given = np.loadtxt(path +data +'/F.txt',delimiter = ',') ## Nf x Nm



# V = V_given.T
# S = S_given[0:500,:].T





#####adding outlier for testing purpose
# for i in range (50):
#     S[0,i] = 1*np.random.random() - 0.5
#######################

# M = M_given.T
# tree_M = KDTree(M.T)
# SigmaS = SigmaS_given
# F = F_given
# Ns = S.shape[1]
# Nm = M.shape[1]
# Nv = V.shape[1]

# Nf = F.shape[0]
# num_sampled_model_points = Nm/300
# M_sampled = M[:,0::num_sampled_model_points]



# ICP_or_GICP_switch = 1 ## 1 for ICP, 2 for GICP
# ICP_triangle_proj_switch = 1 ## 1 to use triangle projection 0 to not

# show_plot_switch = 1 

# points_per_face = int(M.shape[1]/F.shape[0])
# print points_per_face, " points_per_face"
# t0 = time.time()

# tree_M = KDTree(M.T)
# tree_M_sampled = KDTree(M_sampled.T)

# icp_test (V.T,S.T,M.T,tree_M, M_sampled.T, tree_M_sampled,SigmaS,F,points_per_face,ICP_triangle_proj_switch,ICP_or_GICP_switch)
