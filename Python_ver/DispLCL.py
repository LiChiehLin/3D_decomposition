###########################################################################
#                                                                         #
#                             Lin,Li-Chieh                                #
#                      Earth and Planetary Sciences                       #
#                  University of California, Riverside                    #
#                                                                         #
#                                                                         #
# The Python version of codes dveloped in MatLab                          #
# Function names are identical to MatLab's                                #
# Detailed description and usage please refer to GitHub page              #
# https://github.com/LiChiehLin/3D_decomposition                          #
#                                                                         #
#                                                                         #
# Include this .py file in your path and import it as a regular package   #
# Example:                                                                #
#    import DispLCL                                                       #
#    # Use Dsample function, then type                                    #
#    [Output1, Output2] = DispLCL.Dsample(input_variable)                 #
#                                                                         #
# Please report error to:                                                 #
# llin148@ucr.edu                                                         #
###########################################################################
import numpy as np

#--------------Down sampling of GNSS stations (main Dsample)--------------#
# Updates:
    # 2023.08.08 DistMatrix
    # 2023.08.24 Dsample

def Dsample (InMat,Dist):
    # Call DistMatrix to generate distance matrix
    DistMat = DistMatrix(InMat)
    
    OutMat = InMat
    NearSta = np.sum(DistMat<Dist,1)-1
    MaxStaCount = np.max(NearSta)
    if MaxStaCount != 0:
        MaxStaInd = np.asarray(np.where(NearSta == MaxStaCount))
        if MaxStaInd.shape[1] > 1:
            Sta = DistMat[MaxStaInd]
            ndist = Sta*(Sta<Dist)
            ndistAvg = np.sum(ndist[0,:,:],axis=1)/MaxStaCount
            ndistAvgMinInd = np.argmin(ndistAvg)
            OmitInd = MaxStaInd[:,ndistAvgMinInd]
            OutMat = np.delete(OutMat,OmitInd,axis=0)
        else:
            OmitInd = MaxStaInd
            OutMat = np.delete(OutMat,OmitInd,axis=0)
            
    return OutMat, OmitInd
    
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






    
    
    