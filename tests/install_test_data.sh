ACCESSIONs=(NC_005088 NC_007337 NC_008459 NC_001735)

for ACCESSION in ${ACCESSIONs[@]}; do
    curl -L "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$ACCESSION&rettype=fasta&retmode=text" > ./data/fna/$ACCESSION.fasta
done

cat *.fasta > test_data.fna
rm *.fasta
