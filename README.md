main.cpp has the code for Memory Mapping, it reads 2 text files with information about logical memory 
one file has the circuit# id depth width mode of the logical memory as this:
circuit#	ID		mode			depth	width
64			581		SimpleDualPort	128		32
modes can be SimpleDualPort, TrueDualPort,SinglePort, ROM
and the other file contains the number of logic blocks needed, (STRATIX IV including 10 6LUT in each BLOCK) as follow
circuit#	logic block count
6			1853	
having these 2 file after make by running the following the output file will produce 

 ./mapper 	-l 1 1 -b 8192 32 10 1 -b 131072 64 200 1 logical_rams.txt logic_block_count.txt physical_rams.txt 
./mapper is the exe file 

-l 1 1 means LUTRAM with ratio 1:1
-b 8192 32 10 1 means physical RAM with max size 8192 and max width 32 (not for TrueDualPort)
3 file names, 2 for inputs and last one for output file

code now can support max 3 types of RAM but by changing phRAMs global variable to vector more types can be accepted