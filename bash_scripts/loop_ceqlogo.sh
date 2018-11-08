BLAH=1
while ./ceqlogo -i$BLAH meme_format2.motifs -o logos/logo$BLAH.png -f PNG
do 
    BLAH=$((BLAH+1))
done
