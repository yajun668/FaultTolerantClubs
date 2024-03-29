# Code for Finding Fault-tolerant s-Clubs


This code accompanies the paper "[On fault-tolerant low-diameter clusters in graphs](https://pubsonline.informs.org/doi/10.1287/ijoc.2022.1231)" and is written in C++. If you wish to use or cite this code, please cite the paper: 

    @article{LSBB2022fault-tolerant, 
		title = {On fault-tolerant low-diameter clusters in graphs}, 
		author = {Lu, Yajun and Salemi, Hosseinali and Balasundaram, Balabhaskar and Buchanan, Austin}, 
		journal = {INFORMS Journal on Computing},
		pages = {3181-3199},
		volume = {34}, 
		number = {6},
		year = {2022},
		doi = {https://doi.org/10.1287/ijoc.2022.1231},
	}

This repository includes:
1. r_robust_s_club: used to find r-robust s-clubs
2. h_hereditary_s_club: used to find h-hereditary s-clubs

## Datasets

In this paper, we use the below graph instances as our testbed. To keep it self-contained, we include all datasets used in our study in the "Dataset" folder.  The copyrights are credited to the original authors of datasets.

1. Gendreau instances used in [Veremyev and Boginski 2012](https://www.sciencedirect.com/science/article/pii/S0377221711009477).
2. Watts-Strogatz (WS) instances used in [Veremyev et al. 2022](https://www.sciencedirect.com/science/article/pii/S0377221721004227).
3. DIMACS-10 instances in [Bader et al. 2013](http://www.ams.org/books/conm/588/) (also publicly available [here](https://www.cc.gatech.edu/dimacs10/archive/clustering.shtml)).
4. Other Real-Life (ORL) instances used in [Veremyev et al. 2022](https://www.sciencedirect.com/science/article/pii/S0377221721004227) (also publicly available [here](https://sites.pitt.edu/~droleg/files/2-clubs.html)).

## Compiling the code
The following steps show how to compile and run the code to find r-robust s-clubs in a Linux environment using a makefile (you can also run the code in Mac or Windows environment by configuring your IDE appropriately). In the folder of r_robust_s_club, parameter.txt is used to configure parameters r and s, and each graph instance's folder includes InputFile.txt (used to determine which instances you want to test). It is similar to run the code to find h-hereditary s-clubs.


### Steps to run the code to find r-robust s-clubs in Linux environment:
1. Download or clone the repository to your machine.
2. Go to the folder "r_robust_s_club", open the "Makefile" and set GUROBI_HOME to the directory of your Gurobi installation, e.g.: /opt/gurobi/9.0.1/linux64.
3. From the terminal, go to the folder of "r_robust_s_club".
4. Type "make" and hit enter to compile. 
5. Type "./main" to run the code.


## Acknowledgments
We would like to thank [Alexander Veremyev](https://www.cecs.ucf.edu/faculty/alexander-veremyev) for sharing his test instances, codes for graph visualization, and solvers for the maximum r-robust 2-club problem.


## Terms and Use:

MIT License

Copyright (c) 2021 Yajun Lu, Hosseinali Salemi, Balabhaskar Balasundaram, and Austin Buchanan.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
