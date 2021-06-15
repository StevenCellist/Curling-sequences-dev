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

GCC compile instructions for Sequences-negative.cpp (making use of PGO):  
First compile: `g++ -Ofast -march=native -fprofile-generate -pthread Sequences-negative.cpp -o final.out`  
Second compile: `g++ -Ofast -march=native -fprofile-use -fprofile-correction -flto -pthread Sequences-negative.cpp -o final.out`  

MPICXX compile instructions for MPI_Sequences.cpp (making use of PGO):  
First compile: `mpicxx -Ofast -march=native -fprofile-generate MPI_Sequences.cpp -o final.out`  
Second compile: `mpicxx -Ofast -march=native -fprofile-use -fprofile-correction -flto MPI_Sequences.cpp -o final.out`

Made by: Steven Boonstoppel, Vladimir Feinstein and Levi van de Pol
