# Perform Bit-Plane slicing on LOS velocity/displacement field   
Bit-plane slicing is a technique widely used in **Digital Image Processing** field for recognizing which bit plane contributes the most to form the overall image. It is also used for image compression. Typicaally, the image can be well reconstructed using bit planes from `5` to `8`. Bit plane 8 (**Most Significant Bit, MSB**) contains the overall pattern and the subsequent planes carry more finer details (Could be seen in the following example). The low bit planes could potentially be used to denoise the LOS velocity. 

In geophysics, some deformation patterns/traits might be separated out using **Bit-Plane slicing** to isolate the signal we want and to exclude the irrelevant parts. This process includes **3 steps** and **3 Matlab routines** as shown below:  

---
### BPslice.m
Perform the Bit-Plane slicing on the converted ***"Gray"*** scale image
##### Input variable:  
  * GrayMatrix: matrix. An MxN matrix whose values range from 0 to 255.

##### Output variable:
  * BitPlanes: cell. A 8x1 cell which contains the bit planes. 1 is the `Least Significant Bit (LSB)` 8 is the `Most Significant Bit (MSB)`
  * Reconstruct: matrix. An MxN matrix of the reconstructed image from the sliced bit planes

---
### BPrecon.m
Reconstruct the original matrix based on the bit planes in interest. This would require visual inspection to determine what bit planes to use
##### Input variable:  
  * OrigMatrix: matrix. An MxN matrix that was used in `BPslice.m` or the original LOS velocity/displacement
  * BitPlanes: cell. The output ***BitPlanes*** from `BPslice.m`
  * BitPlaneInd: vector. Contains the bit planes to use

##### Output variable:
  * Reconstruct: matrix. The reconstructed matrix based on the bit planes in interest

---
### BPshow.m
A simple function to show the result of Bit-Plane slicing. It takes time for Matlab to show the figure if the input matrix is large (~20 seconds). 

---


### Example:
First, convert the LOS velocity/displacement into [0~255]
```matlab
[~,~,z] = grdread2('filename');
GrayScale = round(255.*rescale(z));
```
![Example](https://github.com/LiChiehLin/3D_decomposition/blob/f658d2f0e529261f8a8cb85f07ad8e193cc32595/Figures/GrayScale_convert_Example.png)

Second, perform Bit-Plane Slicing using `BPslice.m`. You can also view the slicing using `BPshow.m`
```matlab
% Bit-Plane Slicing
[BitPlanes,~] = BPslice(Gray);
% Show the sliced result
BPshow(BitPlanes)
```
![Example](https://github.com/LiChiehLin/3D_decomposition/blob/f658d2f0e529261f8a8cb85f07ad8e193cc32595/Figures/Bit-plane_slice_Example.png)

Lastly, use `BPrecon.m` to reconstruct the LOS velocity. Do a visual inspect on what planes you want to use, or there might exist some more delicate and smarter way. In this case, I chose bit plane `5` and `7` to reconstruct the LOS velocity.
```matlab
BPInd = [5 7];
VeloRecon = BPrecon(z,BitPlanes,BPInd);
```
![Example](https://github.com/LiChiehLin/3D_decomposition/blob/f658d2f0e529261f8a8cb85f07ad8e193cc32595/Figures/Reconstrcuted_Example.png)

### Note that, the first order tectonic activity (fault creep) is separated out using `bit plane 5` and `bit plane 7`.


