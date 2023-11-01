# Smith-Waterman solver for local alignment
### Created with  
![Static Badge](https://img.shields.io/badge/Python-yellow?logo=python&logoColor=white&labelColor=%233776AB&link=https%3A%2F%2Fwww.python.org%2F)

## About Smith-Waterman algorithm
The Smith–Waterman algorithm performs local sequence alignment to determine similar regions between two strings of nucleic acid sequences or protein sequences. When we consider the lengths of the two strings to be n and m, the naive approach has a runtime of O(n3m3), while the Smith-Waterman algorithm has a runtime O(nm).

<br>

## Preparation
To start using our program, download the zip code or clone repository locally with
```sh
git clone https://github.com/nikakoga/Bioinformatic.git
``` 
<br>

# User guide
## Input file
The program needs a file in FASTA format to work. The file should contain two sequences. When there are too few sequences, the program will terminate. If there are more, it will take the first two. Remember to provide the correct path to the file.

## Output file
Here, enter the name of the file in which you would like to save the program result. The output file will be in the same place as the program unless you provide another path.

## Configuration
You can adjust parameters such as mismatch penalty, gap penalty and match reward. Depending on the configuration of these parameters, you can obtain results that are more adequate to your needs. If these parameters are not provided, the algorithm will assume the default values: 

`missmatch: -1 `

 `gap: -1` 
 
 `match: 1`
