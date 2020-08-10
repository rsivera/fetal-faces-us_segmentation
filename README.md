# Segmentation algorithm for Fetal-faces Ultrasound
Atlas-based segmentation algorithm for ultrasound images of fetal faces 

We propose an implementation of the method described by Zuluaga et al. 2013<sup>1</sup>. 

The method requires a task-specific atlas of segmented images that is not included. 

### Dependencies
The script relies on Nifty tools for registration and segmentation.

nifty-reg: https://github.com/KCL-BMEIS/niftyreg

nifty-seg: https://github.com/KCL-BMEIS/NiftySeg


### References
<sup>1</sup> Zuluaga MA, Cardoso MJ, Modat M, Ourselin S. Multi-atlas propagation whole heart
segmentation from MRI and CTA using a local normalised correlation coefficient
criterion. Lect Notes Comput Sci. 2013; 7945 LNCS:174–81.
