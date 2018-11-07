### Uracilation Index ###

* First extract the nucleotide count using bam-readcount:

bam-readcount -w0 -f <ref.fa> <bam_file> | awk -F ":|\t|=" 'BEGIN {OFS = "\t"}; {print $1, $2, $3 , $4, $21 , $35, $49 , $63}' > out.txt

# output columns :
("Chr.","POS","REF","COV","A","C","G","T")

* Run the U-idx Rscript 
