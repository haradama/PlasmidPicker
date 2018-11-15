ACCESSIONs=(NZ_CP007509 NZ_CP007510 NZ_CP008957 NZ_CP008958)

for ACCESSION in ${ACCESSIONs[@]}; do
    curl -L "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=$ACCESSION&rettype=fasta&retmode=text" > ./data/fna/$ACCESSION.fasta
done

cat *.fasta > test_data.fna
rm *.fasta
