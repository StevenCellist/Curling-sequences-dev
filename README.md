# Curling-sequences

On curling sequences, trying to find a sequence with infinite tail length. Let's crunch some sequences!

Largest generator length achieved: *n* = 345.

Record values:
| *n* |&Omega;(*n*)| | *n* |&Omega;(*n*)| | *n* |&Omega;(*n*)|
|:---:|:---:|---|:---:|:---:|---|:---:|:---:|
|2    |2    |   |48   |131  |   |133  |342  |
|4    |4    |   |68   |132  |   |149  |343  |
|6    |8    |   |73   |133  |   |154  |356  |
|8    |58   |   |77   |173  |   |176  |406  |
|9    |59   |   |85   |178  |   |197  |1668 |
|10   |60   |   |115  |215  |   |199  |1669 |
|11   |112  |   |116  |228  |   |200  |1670 |
|14   |118  |   |118  |229  |   |208  |1708 |
|19   |119  |   |128  |332  |   |217  |1836 |
|22   |120  |   |132  |340  |   |290  |3382 |

Inspiration: http://neilsloane.com/doc/CNC.pdf

Disproval of P3: &Omega;(*115*) = 215  
Disproval of P4: &Omega;(*73*) = 133

**File instructions:**  
The default file is Sequences-negative.cpp. This can be run either on Windows (Visual Studio recommended) or Linux (using GCC). Its implementation, however, is limited to one CPU only. For extensive calculations, there is MPI_Sequences.cpp. This file has identical algorithms, but uses the MPI standard in order to enable the use of an HPC or Beowolf cluster. Below, you can find the flags that currently show the best performance on the GCC compiler. Note that MPICXX is just a wrapper of GCC that enables the use of the MPI standard. These flags make use of Profile-Guided Optimization. This requires two compiles: the first compile can use a shorter length (e.g. 200, depending on the cluster), while the second compile should be kept as identical as possible, but maybe a larger length or higher limit.

_GCC compile instructions for Sequences-negative.cpp:_  
First compile: `g++ -Ofast -march=native -fprofile-generate -pthread Sequences-negative.cpp -o final.out`  
Second compile: `g++ -Ofast -march=native -fprofile-use -fprofile-correction -flto -pthread Sequences-negative.cpp -o final.out`  

_MPI(CH) compile instructions for MPI_Sequences.cpp:_  
First compile: `mpicxx -Ofast -march=native -fprofile-generate MPI_Sequences.cpp -o final.out`  
Second compile: `mpicxx -Ofast -march=native -fprofile-use -fprofile-correction -flto MPI_Sequences.cpp -o final.out`

Contributors: Steven Boonstoppel, Vladimir Feinstein, Levi van de Pol
