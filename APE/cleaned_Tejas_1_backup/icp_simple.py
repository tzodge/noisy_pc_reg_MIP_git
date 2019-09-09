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
import time
import transforms3d


def _icp_find_rigid_transform(p_from, p_target):

    A, B = np.copy(p_from), np.copy(p_target)


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
    return R, t


class ICP:
    """ Estimate a rigid-body transform g such that:
        p0 = g.p1
        p0 = model
        p1 = sensor
    """ 
    def __init__(self, p0, p1,tree_M_sampled):

        self.p0 = p0
        self.p1 = p1
        self.nearest = tree_M_sampled
        self.g_series = None

    def compute(self, max_iter):
        ftol = 1.0e-4
        dim_k = self.p0.shape[1]
        g = np.eye(dim_k + 1, dtype=self.p0.dtype)
        p = np.copy(self.p1)

        self.g_series = np.zeros((max_iter + 1, dim_k + 1, dim_k + 1), dtype=g.dtype)
        self.g_series[0, :, :] = g

        itr = -1
        for itr in range(max_iter):

            neighbor_idx = self.nearest.query(p)[1]
            targets = self.p0[neighbor_idx]
            R, t = _icp_find_rigid_transform(p, targets)
            new_p = np.dot(R, p.T).T + t


            if np.linalg.det(R-np.eye(3)) < ftol and np.linalg.norm(t) < ftol:
                break



            p = np.copy(new_p)
            # dg = _icp_Rt_to_matrix(R, t)
            dg = np.eye(4)
            dg[0:3,0:3] = R
            dg[0:3,3] = t
            new_g = np.dot(dg, g)
            g = np.copy(new_g)
            self.g_series[itr + 1, :, :] = g


        self.g_series[(itr+1):, :, :] = g


        return g, (itr + 1)


def icp_test(S,M_sampled,tree_M_sampled):
    #### S = 3 x sampled_sensor_points_for_icp

    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D

    icp = ICP(M_sampled, S, tree_M_sampled)
    matrix, itr= icp.compute(20)



    plot_switch = 0
    if plot_switch == 1:
        
        # error_mat =  matrix.dot(np.linalg.inv(gt_given))
        # angle_error = np.asarray(transforms3d.euler.mat2euler(error_mat))*180/np.pi
        # print "orientation before ICP", np.linalg.norm(np.asarray(transforms3d.euler.mat2euler(gt_given))*180/np.pi)
        # print  "orientation after ICP",np.linalg.norm(angle_error)


    
        tra = np.reshape(np.transpose(matrix[0:3,3]),(3,1))
        transformed = np.matmul(matrix[0:3,0:3], np.transpose(S)) + tra    
        fig = plt.figure()
        ax = Axes3D(fig)
        ax = fig.add_subplot(111, projection='3d')
        ax.set_label("x - axis")
        ax.set_label("y - axis")
        ax.set_label("z - axis")
        ax.plot(M_sampled[:,0], M_sampled[:,1], M_sampled[:,2], "o", color="red", ms=4, mew=0.5, label="M")
        ax.plot(transformed[0,:],transformed[1,:],transformed[2,:],"o",color="green", ms=4, mew=0, label = "transformed")
        ax.plot(S[:,0],S[:,1],S[:,2],"o",color="blue", ms=4, mew=0,label="S")
        plt.legend()
        plt.show()        



    print itr, "itr" 
    return matrix[0:3,0:3] , np.reshape(matrix[0:3,3],(3,1))  


###########EOF






#########for testing purposes
if __name__ ==  '__main__':
    
    data = 'bunny1000_pts_10_deg'
    path = '/home/biorobotics/Desktop/tejas/gurobiCodesPython/tejas5_21_jan/faceNeglect/datasets'+ '/'

    S_given = np.loadtxt(path +data +'/S.txt',delimiter = ',')
    M_given = np.loadtxt(path +data +'/Msampled.txt',delimiter = ',')
    gt_given = np.loadtxt(path +data +'/gt.txt',delimiter = ',')

    S = S_given[0:500,:].T
    M = M_given.T

    Ns = S.shape[1]
    Nm = M.shape[1]

    num_sampled_model_points = Nm/300
    M_sampled = M[:,0::num_sampled_model_points]

    tree_M_sampled = KDTree(M_sampled.T)

    t0 = time.time()

    icp_test (S.T,M_sampled.T,tree_M_sampled )
    ## input to ICP
    ## S = Ns*3
    ## M = Nm*3
    ## tree = KDTree of M

    print time.time() - t0, "time for icp "
    print gt_given,"gt"
     

