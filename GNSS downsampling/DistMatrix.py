###########################################################################
#                                                                         #
#                             Lin,Li-Chieh                                #
#                      Earth and Planetary Sciences                       #
#                  University of California, Riverside                    #
#                              2023.08.08                                 #
#                                                                         #
#                                                                         #
# Calculate distances between each input points. (e.g. GNSS stations)     #
#                                                                         #
# Input:                                                                  #
# 1. InMat: MxN matrix with local coordinates. First two columns must be  #
# coordinates (X,Y).                                                      #                                           
#                                                                         #
# Output:                                                                 #
# 1. OutMat: Distance matrix of each input points                         #
###########################################################################
import numpy as np
def DistMatrix (InMat):
    OutMat = np.zeros((InMat.shape[0],InMat.shape[0]))
    for i in range(InMat.shape[0]):
        tmpMat = InMat[:,(0,1)]
        tmprow = tmpMat.shape[0]
        Now = tmpMat[i,:]
        NowMat = np.ones((tmprow,1))*Now
        dist = np.sqrt(np.sum((NowMat-tmpMat)**2, axis=1, keepdims=True))
        disttmp = dist*np.ones((1,tmprow))
        OutMat[i,:] = disttmp[:,i]
    return OutMat
