#!/bin/bash
make mapper
for (( s=1*1024; s<=128*1024; s= s*2 ))
do  
	for (( w=4; w<=256; w= w*2 ))
	do  
		for (( r=5; r<=90; r= r+5 ))
		do
		./mapper -b $s $w $r 1 logical_rams.txt logic_block_count.txt physical_rams.txt 	
		echo "config size: $s width: $w ratio: $r" 
		./checker -b $s $w $r 1 logical_rams.txt logic_block_count.txt physical_rams.txt
		done
	done
done


