# get CARBS !

![CARBS temp logo](/images/carbs.jpg)
### About
CARBS is a python wrapper for running **CA**dnano **R**igid **B**ody **S**imulations with HOOMD-blue. The first working version is written by Pablo Damasceno. 

Contributions by Tural Aksel (2017)

1. Code is organized into more classes, 
2. New nucleotide types are added for multiscale modeling.
3. Rigid and soft regions in DNA origami design are clustered using depth-first-search graph traversal algorithm.
4. Different spring constats are added for long- and short-range soft connections.
5. Skips in cadnano design files are handled properly. 


### Authors
[Pablo F. Damasceno](http://pablodamasceno.com/). San Francisco, CA (2017) - [Tural Aksel](https://www.linkedin.com/in/tural-aksel-28b97a15/). Douglas Lab, UCSF, 2017.

### Dependencies
- cadnano (https://github.com/cadnano/cadnano2.5)
- hoomd-blue (http://glotzerlab.engin.umich.edu/hoomd-blue/download.html)
- numpy quaternion (pip install numpy-quaternion)
