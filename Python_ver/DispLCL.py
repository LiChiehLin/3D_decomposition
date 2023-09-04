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
import math
import warnings

# Updates:
    # 2023.08.08 DistMatrix
    # 2023.08.24 Dsample
    # 2023.08.26 Deno
        # Dependencies:
            # histbin
            # sub2ind
    # 2023.09.04 nearneighbor
        # Dependencies:
            # ind2sub
#--------------Down sampling of GNSS stations (main Dsample)--------------#
def Dsample (InMat,Dist):
    # Call DistMatrix to generate distance matrix
    DistMat = DistMatrix(InMat)
    
    OutMat = InMat.copy()
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

#----------------------------Denoise (Deno)-------------------------------#
def Deno(NoiseMat,ws):
    R = NoiseMat.shape[0]
    C = NoiseMat.shape[1]
    Grd = NoiseMat.copy()
    
    # Determine the edge of the matrix with respect to the window size
    Edge = math.floor(ws/2)
    # Construct boundary row and column
    RowFirst = []
    RowEnd = []
    ColFirst = []
    ColEnd = []
    for i in range(Edge):
        RowFirst.append(i)
        RowEnd.append((R-1) - i)
        ColFirst.append(i)
        ColEnd.append((C-1) - i)
        
    EdgeRow = np.concatenate((RowFirst,RowEnd[::-1]),axis=None)
    EdgeCol = np.concatenate((ColFirst,ColEnd[::-1]),axis=None)

    # Determine the relative noise-free window
    DiffAvg = []
    for r in range(Edge,R-Edge):
        for c in range(Edge,C-Edge):
            IndexR = list(range(r-Edge,r+Edge+1))
            IndexC = list(range(c-Edge,c+Edge+1))
            
            CenterPix = Grd[r,c]
            MatCal = Grd[IndexR[0]:IndexR[-1]+1,IndexC[0]:IndexC[-1]+1].flatten(order='F')
            MatCal = np.delete(MatCal, np.isnan(MatCal))
            Diff = np.abs(CenterPix - MatCal)
            with warnings.catch_warnings():
                # This warning is expected, so it's fine
                warnings.filterwarnings(action='ignore', message='Mean of empty slice')
                DiffAvg.append(np.mean(Diff))
            
            
    DiffAvg = np.asarray(DiffAvg)
    DiffAvg = np.delete(DiffAvg,np.isnan(DiffAvg))
    
    # Compute CDF of DiffAvg to determine the Tolerance
    Bins = histbin(DiffAvg)
    Counts = np.histogram(DiffAvg,bins=Bins)[0]
    Bins = Bins[1:]
    
    CDF = np.cumsum(Counts)/np.sum(Counts)
        
    # Tolerance is the differnce of 95%
    Index = np.abs(CDF - 0.95)
    Index = np.argmin(Index)
    Tol = Bins[Index]
    
    out = np.zeros(Grd.shape)
    outind = []
    # Deal with inner pixels first
    for r in range(Edge,R-Edge):
        for c in range(Edge,C-Edge):
            # Window that the Pix is centered on
            IndexR = list(range(r-Edge,r+Edge+1))
            IndexC = list(range(c-Edge,c+Edge+1))
            CenterPix = Grd[r,c]
            MatCal = Grd[IndexR[0]:IndexR[-1]+1,IndexC[0]:IndexC[-1]+1].flatten(order='F')
            MatCal = np.delete(MatCal, np.isnan(MatCal))
            Diff = np.abs(CenterPix - MatCal)
            DiffAvg = np.mean(Diff.flatten(order='F'))
            if DiffAvg > Tol:
                # Difference larger than tolerance
                CenterPix = float('nan')
                out[r,c] = CenterPix
                outind.append(sub2ind(out.shape,r,c))
            else:
                # Difference smaller than tolerance
                out[r,c] = CenterPix
                    
    # Deal with boundary pixels
    for r in range(R):
        if np.any(r == EdgeRow):
            for c in range(C):
                if np.asarray(np.where(r == EdgeRow))[0] <= (len(EdgeRow)/2):
                    IndexR = np.asarray(list(range(r-Edge,r+Edge+1)))
                    IndexR = np.delete(IndexR, IndexR < 0)
                    if np.asarray(np.where(c == EdgeCol))[0] <= (len(EdgeCol)/2):
                        # Upperleft corner
                        IndexC = np.asarray(list(range(c-Edge,c+Edge+1)))
                        IndexC = np.delete(IndexC, IndexC < 0)
                        CenterPix = Grd[r,c]
                        MatCal = Grd[IndexR[0]:IndexR[-1]+1,IndexC[0]:IndexC[-1]+1].flatten(order='F')
                        MatCal = np.delete(MatCal, np.isnan(MatCal))
                        Diff = np.abs(CenterPix - MatCal)
                        DiffAvg = np.mean(Diff.flatten(order='F'))
                        if DiffAvg > Tol:
                            # Difference larger than tolerance
                            CenterPix = float('nan')
                            out[r,c] = CenterPix
                            outind.append(sub2ind(out.shape,r,c))
                        else:
                            # Difference smaller than tolerance
                            out[r,c] = CenterPix
                            
                    elif np.asarray(np.where(c == EdgeCol))[0] > (len(EdgeCol)/2):
                        # Upperright corner
                        IndexC = np.asarray(list(range(c-Edge,c+Edge+1)))
                        IndexC = np.delete(IndexC, IndexC >= C)
                        CenterPix = Grd[r,c]
                        MatCal = Grd[IndexR[0]:IndexR[-1]+1,IndexC[0]:IndexC[-1]+1].flatten(order='F')
                        MatCal = np.delete(MatCal, np.isnan(MatCal))
                        Diff = np.abs(CenterPix - MatCal)
                        DiffAvg = np.mean(Diff.flatten(order='F'))
                        if DiffAvg > Tol:
                            # Difference larger than tolerance
                            CenterPix = float('nan')
                            out[r,c] = CenterPix
                            outind.append(sub2ind(out.shape,r,c))
                        else:
                            # Difference smaller than tolerance
                            out[r,c] = CenterPix
                            
                    else:
                        # Upper boundary
                        IndexC = np.asarray(list(range(c-Edge,c+Edge+1)))
                        CenterPix = Grd[r,c]
                        MatCal = Grd[IndexR[0]:IndexR[-1]+1,IndexC[0]:IndexC[-1]+1].flatten(order='F')
                        MatCal = np.delete(MatCal, np.isnan(MatCal))
                        Diff = np.abs(CenterPix - MatCal)
                        DiffAvg = np.mean(Diff.flatten(order='F'))
                        if DiffAvg > Tol:
                            # Difference larger than tolerance
                            CenterPix = float('nan')
                            out[r,c] = CenterPix
                            outind.append(sub2ind(out.shape,r,c))
                        else:
                            # Difference smaller than tolerance
                            out[r,c] = CenterPix
                            
                else:
                    IndexR = np.asarray(list(range(r-Edge,r+Edge+1)))
                    IndexR = np.delete(IndexR, IndexR >= R)
                    if np.asarray(np.where(c == EdgeCol))[0] <= (len(EdgeCol)/2):
                        # Lowerleft corner
                        IndexC = np.asarray(list(range(c-Edge,c+Edge+1)))
                        IndexC = np.delete(IndexC, IndexC < 0)
                        CenterPix = Grd[r,c]
                        MatCal = Grd[IndexR[0]:IndexR[-1]+1,IndexC[0]:IndexC[-1]+1].flatten(order='F')
                        MatCal = np.delete(MatCal, np.isnan(MatCal))
                        Diff = np.abs(CenterPix - MatCal)
                        DiffAvg = np.mean(Diff.flatten(order='F'))
                        if DiffAvg > Tol:
                            # Difference larger than tolerance
                            CenterPix = float('nan')
                            out[r,c] = CenterPix
                            outind.append(sub2ind(out.shape,r,c))
                        else:
                            # Difference smaller than tolerance
                            out[r,c] = CenterPix
                            
                    elif np.asarray(np.where(c == EdgeCol)) > (len(EdgeCol)/2):
                        # Lowerright corner
                        IndexC = np.asarray(list(range(c-Edge,c+Edge+1)))
                        IndexC = np.delete(IndexC, IndexC >= C)
                        CenterPix = Grd[r,c]
                        MatCal = Grd[IndexR[0]:IndexR[-1]+1,IndexC[0]:IndexC[-1]+1].flatten(order='F')
                        MatCal = np.delete(MatCal, np.isnan(MatCal))
                        Diff = np.abs(CenterPix - MatCal)
                        DiffAvg = np.mean(Diff.flatten(order='F'))
                        if DiffAvg > Tol:
                            # Difference larger than tolerance
                            CenterPix = float('nan')
                            out[r,c] = CenterPix
                            outind.append(sub2ind(out.shape,r,c))
                        else:
                            # Difference smaller than tolerance
                            out[r,c] = CenterPix
                            
                    else:
                        # Lower boundary
                        IndexC = np.asarray(list(range(c-Edge,c+Edge+1)))
                        CenterPix = Grd[r,c]
                        MatCal = Grd[IndexR[0]:IndexR[-1]+1,IndexC[0]:IndexC[-1]+1].flatten(order='F')
                        MatCal = np.delete(MatCal, np.isnan(MatCal))
                        Diff = np.abs(CenterPix - MatCal)
                        DiffAvg = np.mean(Diff.flatten(order='F'))
                        if DiffAvg > Tol:
                            # Difference larger than tolerance
                            CenterPix = float('nan')
                            out[r,c] = CenterPix
                            outind.append(sub2ind(out.shape,r,c))
                        else:
                            # Difference smaller than tolerance
                            out[r,c] = CenterPix
                            
        else:
            for c in EdgeCol:
                IndexR  = np.asarray(list(range(r-Edge,r+Edge+1)))
                if np.asarray(np.where(c == EdgeCol)) <= (len(EdgeCol)/2):
                    # Left
                    IndexC = np.asarray(list(range(c-Edge,c+Edge+1)))
                    IndexC = np.delete(IndexC, IndexC < 0)
                    CenterPix = Grd[r,c]
                    MatCal = Grd[IndexR[0]:IndexR[-1]+1,IndexC[0]:IndexC[-1]+1].flatten(order='F')
                    MatCal = np.delete(MatCal, np.isnan(MatCal))
                    Diff = np.abs(CenterPix - MatCal)
                    DiffAvg = np.mean(Diff.flatten(order='F'))
                    if DiffAvg > Tol:
                        # Difference larger than tolerance
                        CenterPix = float('nan')
                        out[r,c] = CenterPix
                        outind.append(sub2ind(out.shape,r,c))
                    else:
                        # Difference smaller than tolerance
                        out[r,c] = CenterPix
                        
                else:
                     # Right
                     IndexC = np.asarray(list(range(c-Edge,c+Edge+1)))
                     IndexC = np.delete(IndexC, IndexC >=C )
                     CenterPix = Grd[r,c]
                     MatCal = Grd[IndexR[0]:IndexR[-1]+1,IndexC[0]:IndexC[-1]+1].flatten(order='F')
                     MatCal = np.delete(MatCal, np.isnan(MatCal))
                     Diff = np.abs(CenterPix - MatCal)
                     DiffAvg = np.mean(Diff.flatten(order='F'))
                     if DiffAvg > Tol:
                         # Difference larger than tolerance
                         CenterPix = float('nan')
                         out[r,c] = CenterPix
                         outind.append(sub2ind(out.shape,r,c))
                     else:
                         # Difference smaller than tolerance
                         out[r,c] = CenterPix
                         
    DenoiseMat = out
    NInd = outind
    return DenoiseMat, NInd
                 
def histbin(data):
    # Calculate histogram bins edges and pass to np.histogram
    # This binning method is referenced from MatLab and this discussion:
    # https://au.mathworks.com/matlabcentral/answers/396614-how-is-the-number-of-bins-chosen-with-the-auto-binning-algorithm-in-histcounts
    
    rawBinWidth = 3.5*np.std(data)/len(data)**(1/3)
    xmin = np.min(data)
    xmax = np.max(data)
    
    powOfTen = 10**np.floor(np.log10(rawBinWidth))
    relSize = rawBinWidth/powOfTen
    if relSize < 1.5:
        binWidth = 1*powOfTen
    elif relSize < 2.5:
        binWidth = 2*powOfTen
    elif relSize < 4:
        binWidth = 3*powOfTen
    elif relSize < 7.5:
        binWidth = 5*powOfTen
    else:
        binWidth = 10*powOfTen
        
    leftEdge = np.min([binWidth*np.floor(xmin/binWidth), xmin])
    nbinsActual = np.max([1, np.ceil((xmax - leftEdge))/binWidth])
    rightEdge = np.max([leftEdge + nbinsActual*binWidth, xmax])
    bins = np.linspace(int(leftEdge),int(rightEdge),int(nbinsActual)+1)
    return bins
        
def sub2ind(array_shape,rows,cols):
    # Modified from:
    # https://blog.csdn.net/weixin_35193147/article/details/108940367
    ind = rows + cols*array_shape[0]
    ind = int(ind)
    if ind < 0:
        ind = -1
        print("Index Value ERROR!")
    elif ind >= array_shape[0]*array_shape[1]:
        ind = -1
        print("Index Value OVERFLOW!")
    else:
         return ind
     
        
#----------------------------Nearneighbor-------------------------------#
def nearneighbor(InMat,Radius,Weight):
    # Find total missing values and indices
    missind = np.where(np.isnan(InMat.flatten(order='F')))[0]
    missrow,misscol = ind2sub(InMat.shape,missind)
    miss = np.sum(np.isnan(InMat.flatten(order='F')))
    print('Missing values:',str(miss))
    
    out = InMat.copy()
    ##### Search neighbors for each missing values
    for i in range(len(missind)):
        print('Working on',str(i),end='\n')
        # Make sure search box is within matrix
        #Matind = missind[i]
        row = missrow[i]
        rtmp = np.asarray(range(row-Radius,row+Radius+1))
        col = misscol[i]
        ctmp = np.asarray(range(col-Radius,col+Radius+1))
        rowtmp = rtmp[np.all([(rtmp >= 0),rtmp < InMat.shape[0]],axis=0)] # rows to be searched
        coltmp = ctmp[np.all([(ctmp >= 0),ctmp < InMat.shape[1]],axis=0)] # cols to be searched
        
        ##### Start searching neighbors
        # Get neighboring values
        rowind = rowtmp*np.ones([len(coltmp),1])
        rowind = rowind.flatten(order='F')
        colind = coltmp*np.ones([len(rowtmp),1])
        colind = colind.flatten(order='C')
        
        linind = []
        for j in range(len(rowind)):
            linind.append(sub2ind(InMat.shape,rowind[j],colind[j]))
        
        ##### No weighting
        if Weight == 0:
            neival = InMat.flatten(order='F')[linind]
            with warnings.catch_warnings():
                # This warning is expected, so it's fine
                warnings.filterwarnings(action='ignore', message='Mean of empty slice')
                out[missrow[i],misscol[i]] = np.nanmean(neival)
            
        elif Weight == 1:
            ##### Include weighting
            Smseq = np.asarray(range(1,Radius+1,1))
            rlow = row - Smseq
            rlow[rlow < 0] = 0
            rup = row + Smseq
            rup[rup >= np.max(rowind)] = np.max(rowind)
            clow = col - Smseq
            clow[clow < 0] = 0
            cup = col + Smseq
            cup[cup >= np.max(colind)] = np.max(colind)
            

            rowindtmp = np.tile(rowind,(len(rlow),1)).T
            rlowtmp = np.tile(rlow,(len(rowind),1))
            ruptmp = np.tile(rup,(len(rowind),1))
            Smrow = np.all([rowindtmp >= rlowtmp, rowindtmp <= ruptmp],axis=0)
            
            colindtmp = np.tile(colind,(len(clow),1)).T
            clowtmp = np.tile(clow,(len(colind),1))
            cuptmp = np.tile(cup,(len(colind),1))
            Smcol = np.all([colindtmp >= clowtmp, colindtmp <= cuptmp],axis=0)
            
            Smlinind = np.tile(linind,(Radius,1)).T
            Smbool = np.all([Smrow,Smcol],axis=0)
            Smind = Smlinind.flatten(order='F')[Smbool.flatten(order='F')]
            Smneival = InMat.flatten(order='F')[Smind]
            
            with warnings.catch_warnings():
                # This warning is expected, so it's fine
                warnings.filterwarnings(action='ignore', message='Mean of empty slice')
                out[missrow[i],misscol[i]] = np.nanmean(Smneival)
                
        else:
             return -1
             
    return out
        
        
        
def ind2sub(array_shape, ind):
    # Modified from:
    # https://gist.github.com/mizunototori/9b26f8b5bdcccabdadb460da7dcfffc1
    cols = (ind.astype("int32") // array_shape[0])
    rows = (ind.astype("int32") % array_shape[0])
    return (rows, cols)


    
    
    
        
        
        
        
        












    
        
