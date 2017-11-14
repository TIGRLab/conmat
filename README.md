dti_conn
--------

Computes probabilistic connectivity matricies using deterministic and probabilistic methods on diffusion weighted imaging (DTI) data. Relies on [CAMINO](http://camino.cs.ucl.ac.uk/index.php?n=Main.ManPages), and optionally works on the outputs of [BEDPOSTX](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FDT/UserGuide#BEDPOSTX).

**installation**

+ add `dti_conn.py` to your PATH.

**deterministic tractography details**

+ parameters determined from J. T. Duda, P. A. Cook, and J. C. Gee, “Reproducibility of graph metrics of human brain structural networks.,” Front Neuroinform, vol. 8, p. 46, 2014.
+ `wdfit` is used -- a wrapper for modelfit, using RK4 tracking (Fourth-order Runge-Kutta method). see Peter J. Basser, Sinisa Pajevic, Carlo Pierpaoli, Jeffrey Duda, Akram Aldrou. In vivo fiber tractography using DT-MRI data. Magnetic Resonance in Medicine, 2000., and Mariana Lazar, David M. Weinstein, Jay S. Tsuruda, Khader M. Hasan, Konstantinos Arfanakis, M. Elizabeth Meyerand, Benham Badie, Howard A. Rowley, Victor Haughton, Aaron Field, Andrew L. Alexander. White matter tractography using diffusion tensor deflection. Human Brain Mapping, 2003.

**probabilisitc tractography details**
+ uses probability look-up table (generated by `dtlutgen`) followed by PICo tractography, see Parker GJM, Haroon HA and Wheeler-Kingshott CAM, A Framework for a Streamline-Based Probabilistic Index of Connectivity (PICo) using a Structural Interpretation of MRI Diffusion Measurements, Journal of Magnetic Resonance Imaging, 18, 242-254, 2003., Cook PA, Alexander DC, Parker GJM, Modelling noise-induced fibre-orientation error in diffusion-tensor MRI, IEEE International Symposium on Biomedical Imaging, 332-335, 2004., Cook PA, Alexander DC,Modelling uncertainty in two fibre-orientation estimates within a voxel, Proc. Intl. Soc. Mag. Reson. Med. 14 (2006) pp 1629.
+ estimates SNR
+ this stage can take many hours to complete.

**bedpostx details**
+ uses camino's method for interpreting BEDPOSTX as either a deterministic or probabilisitic tractography method, and produces a connectivity matrix.


**connectivity matrix**

+ produced by `conmat`, major output are the `*${method}_sc.csv` files which represent fiber counts between pairs of ROIs.
+ for deterministic tractography, also outputs the max/mean/min FA and MD between each pair of ROIs.

**visualisation**
+ try `vtkstreamlines` or `pdview`.

