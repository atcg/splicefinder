Splicefinder
======

This is the walkthrough and code that accompanies "Conservative prediction of intron splice sites for the design of exon-capture arrays."



###Step 1: Download mouse data from Ensembl

_________________

Connect to the server, navigate to the proper folder, and download data in bash:

```bash
ssh rover
cd /mnt/Data3/arrayDesignPaper
mkdir ensembl
cd ensembl
wget ftp://ftp.ensembl.org/pub/release-80/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz

wget ftp://ftp.ensembl.org/pub/release-80/fasta/mus_musculus/cdna/README

gunzip Mus_musculus.GRCm38.cdna.all.fa.gz
```

After downloading the data, `get_fasta_lengths.py --input Mus_musculus.GRCm38.cdna.all.fa` gives the following output:

```text
Reads:		90,891
Bp:		168,812,502
Avg. len:	1,857.30712612
STDERR len:	6.5011432467
Min. len:	9
Max. len:	106,824
Median len:	1,078.0
Contigs > 1kb:	47,261
```


These represent the latest set of mouse cDNAs from ENSEMBL. We are first going to randomly subset these into fewer than 90,891, because we are going to run many analyses on these sets, mapping to multiple genomes (which can be slow). Let's see what they look like:

```text
head --lines 6 Mus_musculus.GRCm38.cdna.all.fa
>ENSMUST00000196221 havana_ig_gene:known chromosome:GRCm38:14:54113468:54113476:1 gene:ENSMUSG00000096749 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene
ATGGCATAT
>ENSMUST00000179664 ensembl_ig_gene:known chromosome:GRCm38:14:54113468:54113478:1 gene:ENSMUSG00000096749 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene
ATGGCATATCA
>ENSMUST00000177564 ensembl_havana_ig_gene:known chromosome:GRCm38:14:54122226:54122241:1 gene:ENSMUSG00000096176 gene_biotype:TR_D_gene transcript_biotype:TR_D_gene
ATCGGAGGGATACGAG
```

Looking at this output tells us a few things. First of all, some of the transcripts are extremely short. Since the purpose of this method is to infer contiguous genomic sequence chunks long enough to use target enrichment array design, we know we'll want to filter out all cDNAs that are too short to actually yield a target. We'll arbitrarily set this limit at 200bp for now. The following script will output all sequences longer than 199bp:

```perl
#!/usr/bin/perl

#filterUnder200.pl

use strict;
use warnings;
use Bio::SeqIO;

my $seqIn = Bio::SeqIO->new(-file => "Mus_musculus.GRCm38.cdna.all.fa",
                            -format => 'fasta');
my $seqOut = Bio::SeqIO->new(-file => ">mmGRCm38.cdna.200bpPlus.fa",
                             -format => "fasta");

while (my $seq = $seqIn->next_seq()) {
    if ($seq->length() > 199) {
        $seqOut->write_seq($seq);
    }
}
```

After running, we have:

```text
get_fasta_lengths.py --input mmGRCm38.cdna.200bpPlus.fa
Reads:          90,060
Bp:             168,702,266
Avg. len:       1,873.22080835
STDERR len:     6.53762906167
Min. len:       200
Max. len:       106,824
Median len:     1,098.0
Contigs > 1kb:  47,261
```

That didn't very many, so we'll have to reduce further. You'll notice in the fasta headers that there is a gene name in the fourth column. Let's find out how many unique gene names there are total:

```text
grep "gene:" mmGRCm38.cdna.200bpPlus.fa | awk '{print $4}' | uniq -c | wc -l
30722
```

Randomly choosing one transcript from each gene name seems like one good way to reduce the dataset without losing genomic coverage. But choosing the longest cDNA from each gene name might give us better success at mapping to other genomes. Let's do that with the following script:

```perl
#!/usr/bin/perl

# longestFromGenes.pl

use strict;
use warnings;
use Bio::SeqIO;


my %seqHash;
my $seqIn1 = Bio::SeqIO->new(-file => "mmGRCm38.cdna.200bpPlus.fa",
                             -format => "fasta");
my $seqOut = Bio::SeqIO->new(-file => ">mmGRCm38.cdna.longestFromGenes.fa",
                             -format => "fasta");

while (my $seq = $seqIn1->next_seq()) {
    if ($seq->desc() =~ /.*gene\:(\S+)\s/) {
        # Each transcript is associated with a gene name, and has a
        # hash value equal to its length
        my $geneName = $1;
        $seqHash{$geneName}{$seq->display_id()} = $seq;
    }
}

foreach my $gene (sort keys %seqHash) {
    my $longestID;
    my $longestLength = 0;
    foreach my $transcript (sort keys %{$seqHash{$gene}}) {
        if ($seqHash{$gene}{$transcript}->length() > $longestLength) {
            $longestID = $transcript;
            $longestLength = $seqHash{$gene}{$transcript}->length();
        }
    }
    $seqOut->write_seq($seqHash{$gene}{$longestID})
}
```


This should get the list down to 30,722, and the average length should go up (it was 1,873 bp before).

```text
get_fasta_lengths.py --input mmGRCm38.cdna.longestFromGenes.fa
Reads:          30,722
Bp:             73,817,860
Avg. len:       2,402.76869995
STDERR len:     13.0889977284
Min. len:       200
Max. len:       106,824
Median len:     1,734.0
Contigs > 1kb:  20,279
```

Looking good! I think 30,722 transcripts is still a little high though, so I will subset the dataset down to 10,000 transcripts randomly.

```perl
use List::Util qw/shuffle/;

#random10k.pl

my $seqIn = Bio::SeqIO->new(-file => "mmGRCm38.cdna.longestFromGenes.fa",
                            -format => "fasta");
my $seqOut = Bio::SeqIO->new(-file => ">mmGRCm38.cdna.rand10kLongest.fa",
                            -format => "fasta");

my @seqArray;
while (my $seq = $seqIn->next_seq()) {
    push(@seqArray, $seq);
}

my @shuffledSeqs = shuffle(@seqArray);

my $counter = 0;
while ($counter < 10000) {
    $seqOut->write_seq($shuffledSeqs[$counter]);
    $counter++;
}
```

That should give us 10,000 transcripts from 10,000 distinct genes scattered throughout the mouse genome.

```text
get_fasta_lengths.py --input mmGRCm38.cdna.rand10kLongest.fa
Reads:		10,000
Bp:		23,940,240
Avg. len:	2,394.024
STDERR len:	21.8914472689
Min. len:	200
Max. len:	22,489
Median len:	1,732.0
Contigs > 1kb:	6,651
```

Excellent. Now we're ready to set up some mapping. In this experiment, we're going to want to map to approximately 10 genomes of varying sequence divergence from the mouse genome. For now, let's try the following genomes:
  1. Rat (_Rattus norvegicus_)
  2. Human (_Homo sapiens_)
  3. Pig (_Sus scrofa_)
  4. Elephant (_Loxodonta_africana_)
  5. Playtpus (_Ornithorynchus anatinus_)
  6. Chicken (_Gallus gallus_)
  7. _Xenopus tropicalis_
  8. _Anolis carolinensis_
  9. Coelacanth (_Latimeria chalumnae_)
  10. Puffer fish (_Takifugu rubripes_)
  11. Sea lamprey (*Petromyzon marinus*)
  12. Sea urchin (*Strongylocentrotus purpuratus*)


Let's download those genomes, getting the unmasked primary assembly if possible, otherwise taking the unmasked top-level assembly:

```bash
mkdir genomes
cd genomes
wget ftp://ftp.ensembl.org/pub/release-80/fasta/rattus_norvegicus/dna/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa10.2.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/loxodonta_africana/dna/Loxodonta_africana.loxAfr3.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/ornithorhynchus_anatinus/dna/Ornithorhynchus_anatinus.OANA5.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/gallus_gallus/dna/Gallus_gallus.Galgal4.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/xenopus_tropicalis/dna/Xenopus_tropicalis.JGI_4.2.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/anolis_carolinensis/dna/Anolis_carolinensis.AnoCar2.0.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/latimeria_chalumnae/dna/Latimeria_chalumnae.LatCha1.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/takifugu_rubripes/dna/Takifugu_rubripes.FUGU4.dna.toplevel.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/petromyzon_marinus/dna/Petromyzon_marinus.Pmarinus_7.0.dna.toplevel.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-26/fasta/strongylocentrotus_purpuratus/dna/Strongylocentrotus_purpuratus.GCA_000002235.2.26.dna.toplevel.fa.gz

gunzip *.gz*
```
```

We'll also need the full collection of protein sequences from these organisms:
```bash
cd /mnt/Data3/arrayDesignPaper/ensembl/genomes

#!/bin/bash
#downloadProteins.bash
mkdir proteins
cd proteins

wget ftp://ftp.ensembl.org/pub/release-80/fasta/rattus_norvegicus/pep/Rattus_norvegicus.Rnor_6.0.pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/sus_scrofa/pep/Sus_scrofa.Sscrofa10.2.pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/loxodonta_africana/pep/Loxodonta_africana.loxAfr3.pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/ornithorhynchus_anatinus/pep/Ornithorhynchus_anatinus.OANA5.pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/gallus_gallus/pep/Gallus_gallus.Galgal4.pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/xenopus_tropicalis/pep/Xenopus_tropicalis.JGI_4.2.pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/anolis_carolinensis/pep/Anolis_carolinensis.AnoCar2.0.pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/latimeria_chalumnae/pep/Latimeria_chalumnae.LatCha1.pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/takifugu_rubripes/pep/Takifugu_rubripes.FUGU4.pep.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/petromyzon_marinus/pep/Petromyzon_marinus.Pmarinus_7.0.pep.all.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-26/fasta/strongylocentrotus_purpuratus/pep/Strongylocentrotus_purpuratus.GCA_000002235.2.26.pep.all.fa.gz

gunzip *
```

We're going to use exonerate protein2genome to map all of these proteins to their respective genomes.
```bash
cd /mnt/Data3/arrayDesignPaper/ensembl/genomes

#!/bin/bash
#prepareExonerate.bash
cd ../

fasta2esd --softmask no Loxodonta_africana.loxAfr3.dna.toplevel.fa Loxodonta_africana.loxAfr3.dna.toplevel.esd
fasta2esd --softmask no Ornithorhynchus_anatinus.OANA5.dna.toplevel.fa Ornithorhynchus_anatinus.OANA5.dna.toplevel.esd
fasta2esd --softmask no Strongylocentrotus_purpuratus.GCA_000002235.2.26.dna.toplevel.fa Strongylocentrotus_purpuratus.GCA_000002235.2.26.dna.toplevel.esd
fasta2esd --softmask no Petromyzon_marinus.Pmarinus_7.0.dna.toplevel.fa Petromyzon_marinus.Pmarinus_7.0.dna.toplevel.esd
fasta2esd --softmask no Takifugu_rubripes.FUGU4.dna.toplevel.fa Takifugu_rubripes.FUGU4.dna.toplevel.esd
fasta2esd --softmask no Latimeria_chalumnae.LatCha1.dna.toplevel.fa Latimeria_chalumnae.LatCha1.dna.toplevel.esd
fasta2esd --softmask no Anolis_carolinensis.AnoCar2.0.dna.toplevel.fa Anolis_carolinensis.AnoCar2.0.dna.toplevel.esd
fasta2esd --softmask no Xenopus_tropicalis.JGI_4.2.dna.toplevel.fa Xenopus_tropicalis.JGI_4.2.dna.toplevel.esd
fasta2esd --softmask no Gallus_gallus.Galgal4.dna.toplevel.fa Gallus_gallus.Galgal4.dna.toplevel.esd
fasta2esd --softmask no Sus_scrofa.Sscrofa10.2.dna.toplevel.fa Sus_scrofa.Sscrofa10.2.dna.toplevel.esd
fasta2esd --softmask no Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.dna.primary_assembly.esd
fasta2esd --softmask no Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa Rattus_norvegicus.Rnor_6.0.dna.toplevel.esd

esd2esi Loxodonta_africana.loxAfr3.dna.toplevel.esd Loxodonta_africana.loxAfr3.dna.toplevel.trans.esi --translate yes --memorylimit 25600
esd2esi Ornithorhynchus_anatinus.OANA5.dna.toplevel.esd Ornithorhynchus_anatinus.OANA5.dna.toplevel.trans.esi --translate yes --memorylimit 25600
esd2esi Strongylocentrotus_purpuratus.GCA_000002235.2.26.dna.toplevel.esd Strongylocentrotus_purpuratus.GCA_000002235.2.26.dna.toplevel.trans.esi --translate yes --memorylimit 25600
esd2esi Petromyzon_marinus.Pmarinus_7.0.dna.toplevel.esd Petromyzon_marinus.Pmarinus_7.0.dna.toplevel.trans.esi --translate yes --memorylimit 25600
esd2esi Takifugu_rubripes.FUGU4.dna.toplevel.esd Takifugu_rubripes.FUGU4.dna.toplevel.trans.esi --translate yes --memorylimit 25600
esd2esi Latimeria_chalumnae.LatCha1.dna.toplevel.esd Latimeria_chalumnae.LatCha1.dna.toplevel.trans.esi --translate yes --memorylimit 25600
esd2esi Anolis_carolinensis.AnoCar2.0.dna.toplevel.esd Anolis_carolinensis.AnoCar2.0.dna.toplevel.trans.esi --translate yes --memorylimit 25600
esd2esi Xenopus_tropicalis.JGI_4.2.dna.toplevel.esd Xenopus_tropicalis.JGI_4.2.dna.toplevel.trans.esi --translate yes --memorylimit 25600
esd2esi Gallus_gallus.Galgal4.dna.toplevel.esd Gallus_gallus.Galgal4.dna.toplevel.trans.esi --translate yes --memorylimit 25600
esd2esi Sus_scrofa.Sscrofa10.2.dna.toplevel.esd Sus_scrofa.Sscrofa10.2.dna.toplevel.trans.esi --translate yes --memorylimit 25600
esd2esi Homo_sapiens.GRCh38.dna.primary_assembly.esd Homo_sapiens.GRCh38.dna.primary_assembly.trans.esi --translate yes --memorylimit 25600
esd2esi Rattus_norvegicus.Rnor_6.0.dna.toplevel.esd Rattus_norvegicus.Rnor_6.0.dna.toplevel.trans.esi --translate yes --memorylimit 25600
```

The above will take several hours. We want to map the orthologous proteins from each species to their own reference genomes. So, we will need to find the orthologous proteins from each species. First we'll make blast databases from each protein set:
```bash
cd proteins
for i in *.pep.all.fa; do makeblastdb -in $i -dbtype prot; done
```

Now we'll iterate through each of the protein sets and blastx the mouse cDNAs to each, only outputting the best match:

```bash
for i in *.pep.all.fa; do OUTFILE=$i".proteinMatches"; blastx -db $i -query ../../mmGRCm38.cdna.rand10kLongest.fa -outfmt 6 -max_target_seqs 1 -num_threads 10 -out $OUTFILE ; done
```

The blastx'ing will also take a few hours.

After the blastx is done, we'll parse the blastx results to pull out all the protein sequences that matched and store them into a new file for exonerate protein2genome mapping:

```perl
#!/usr/bin/perl

# pullOutMatchingProteins.pl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::SearchIO;

my @speciesArray = ("Strongylocentrotus_purpuratus.GCA_000002235.2.26", "Petromyzon_marinus.Pmarinus_7.0", "Takifugu_rubripes.FUGU4", "Latimeria_chalumnae.LatCha1", "Anolis_carolinensis.AnoCar2.0", "Xenopus_tropicalis.JGI_4.2", "Gallus_gallus.Galgal4", "Ornithorhynchus_anatinus.OANA5", "Loxodonta_africana.loxAfr3", "Sus_scrofa.Sscrofa10.2", "Homo_sapiens.GRCh38", "Rattus_norvegicus.Rnor_6.0");

foreach my $species (@speciesArray) {
    my $proteinFasta = $species . ".pep.all.fa";

    my %protHash;
    my $seqIn = Bio::SeqIO->new(-file => $proteinFasta,
                                -format => 'fasta');
    while (my $seq = $seqIn->next_seq()) {
        $protHash{$seq->display_id()} = $seq;
    }

    my $proteinOut = $species . ".pep.matching.fa";
    my $seqOut = Bio::SeqIO->new(-file => ">$proteinOut",
                                 -format => 'fasta');

    my $blastResults = $species . ".pep.all.fa.proteinMatches";
    open(my $blastFH, "<", $blastResults) or die "Couldn't open $blastResults for reading: $!\n";
    while (my $line = <$blastFH>) {
        my @fields = split(/\t/, $line);
        $seqOut->write_seq($protHash{$fields[1]});
    }
}
```

Now we'll set these aside:

```bash
mkdir matchingProteins
mv "*.pep.matching.fa" matchingProteins/
```

Exonerate protein2genome is very slow unless you change some settings around, and it seems to run faster when you separate each sequence into its own file. Using only one sequence per run also has the benefit of showing you where you left off if exonerate segfaults.

So we're first gonna need to separate the protein fasta files into new files, one sequence per file:

```perl
#!/usr/bin/perl

# separateProteinsIntoIndividuals.pl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::SearchIO;

my @speciesArray = ("Strongylocentrotus_purpuratus.GCA_000002235.2.26", "Petromyzon_marinus.Pmarinus_7.0", "Takifugu_rubripes.FUGU4", "Latimeria_chalumnae.LatCha1", "Anolis_carolinensis.AnoCar2.0", "Xenopus_tropicalis.JGI_4.2", "Gallus_gallus.Galgal4", "Ornithorhynchus_anatinus.OANA5", "Loxodonta_africana.loxAfr3", "Sus_scrofa.Sscrofa10.2", "Homo_sapiens.GRCh38", "Rattus_norvegicus.Rnor_6.0");

foreach my $species (@speciesArray) {
    # We'll have a different folder for each species
    unless (-d $species) {
        mkdir $species;
    }
    chdir $species;

    my $seqsFile = "../" . $species . ".pep.matching.fa";

    my $seqIn = Bio::SeqIO->new(-file => $seqsFile,
                                -format => 'fasta');

    my $counter = 0;
    while (my $seq = $seqIn->next_seq()) {
        $counter++;
        my $seqOutName = $seq->display_id() . ".pep.fasta";
        my $seqOut = Bio::SeqIO->new(-file => ">$seqOutName",
                                     -format => 'fasta');

        $seqOut->write_seq($seq);

    }
    print $counter . " total records processed\n";
    chdir "..";
}
```

Note that some mouse transcripts might have selected the same genes from the reference set as others for the best blastx hit. This means that there should be fewer actual peptide fastas than query sequences that had positive matches in the blastx search.

Now let's loop through all those sequences and align them to their proper genomes. The sending multiple queries to a single exonerate-server tends to make things segfault, so instead we'll instantiate exonerate-servers so that each thread gets its own.


```perl
#!/usr/bin/perl

# exonerateP2G.pl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::SearchIO;
use Parallel::ForkManager;

my @speciesArray = ("Strongylocentrotus_purpuratus.GCA_000002235.2.26", "Petromyzon_marinus.Pmarinus_7.0", "Takifugu_rubripes.FUGU4", "Latimeria_chalumnae.LatCha1", "Anolis_carolinensis.AnoCar2.0", "Xenopus_tropicalis.JGI_4.2", "Gallus_gallus.Galgal4", "Ornithorhynchus_anatinus.OANA5", "Loxodonta_africana.loxAfr3", "Sus_scrofa.Sscrofa10.2", "Homo_sapiens.GRCh38", "Rattus_norvegicus.Rnor_6.0");


my $portCounter = 12886; # We'll start on port 12887 and go up 29 more

foreach my $species (@speciesArray) {
    my $genomeFile;
    if ($species eq "Homo_sapiens.GRCh38") {
        $genomeFile = "Homo_sapiens.GRCh38.dna.primary_assembly.trans.esi";
    } else {
        my $genomeFile = $species . ".dna.toplevel.trans.esi";
    }

    # We'll use a total of 4 possible threads for the forkmanager, because
    # each thread will initiate two processes--the exonerate server and the
    # exonerate clients. We want to use about 8 CPUs total, so set this
    # to 4.
    my $forkManager = Parallel::ForkManager->new(4);

    foreach my $thread (0..3) {
        $portCounter++;
        $forkManager->start and next;
        #chdir("/home/evan/spliceFinder/genomes/");
        system("exonerate-server --port $portCounter $genomeFile &");
        sleep(30);
        chdir("proteins/matchingProteins/");
        chdir $species;
        unless(-d "exonerate") {
            mkdir "exonerate";
        }

        opendir(my $dirFH, "./");
        my @files = readdir($dirFH);
        closedir($dirFH);

        my @fastaFiles;
        foreach my $file (@files) {
            if ($file =~ /.fasta/) {
                push(@fastaFiles, $file);
            }
        }

        my $proteinsPerThread = scalar(@fastaFiles) / 4; # Decide how many proteins to provide to each thread

        my $beginning = $thread * $proteinsPerThread;
        my $end = ($thread+1) * $proteinsPerThread;
        for (my $fileIndex=$beginning; $fileIndex < $end; $fileIndex++) {
            my $outFileName = "exonerate/" . $fastaFiles[$fileIndex] . ".p2g.exonerate";
            system("exonerate $fastaFiles[$fileIndex] localhost:$portCounter --querytype protein --targettype dna --model p2g --bestn 1 --dnawordlen 12 --fsmmemory 4096 --hspfilter 50 --geneseed 250 > $outFileName");
        }
        $forkManager->finish;
    }
    $forkManager->wait_all_children;
    sleep(120);
    system("pkill -f exonerate-server");
    sleep(120);
}
```

This takes several hours to a day to run all these, assuming all goes well and you have > 30 cores to work with.

After this comes one of the trickier parts. We'll process all of these exonerate files to harvest the putative contiguous genomic region sequences and align to the proteins. For instance, rat protein ENSRNOP00000046335 was mapped to the rat chromosome 7 between base pairs 143497108 and 143489162. Additionally, exonerate uncovered 8 introns spread throughout the alignment. The first intron starts at about bp 143496581 and ends at bp 143495791. So, we'll designate the first intron as running from 143496581-143497108. We want to pull the rat genomic DNA sequence from 143496581-143497108 and align it to the original mouse EST sequence.

Let's say this exon sequence aligns at base pairs 10-537 in the mouse EST. We would then code this as rat-inferred splice sites before base 10 and after 537 in the mouse EST. We'll gather all such inferred splice sites for each EST, keeping track of which ones come from which species.

I'll present this all in heavily-commented code below.

```perl
#!/usr/bin/perl

# harvestSplicePositions.pl

use strict;
use warnings;
use Bio::SearchIO; # <= We'll use this to parse the exonerate reports
use Bio::SeqIO;


# We'll deal with all the different reference genomes one-by-one in a loop
my @speciesArray = ("Strongylocentrotus_purpuratus.GCA_000002235.2.26", "Petromyzon_marinus.Pmarinus_7.0", "Takifugu_rubripes.FUGU4", "Latimeria_chalumnae.LatCha1", "Anolis_carolinensis.AnoCar2.0", "Xenopus_tropicalis.JGI_4.2", "Gallus_gallus.Galgal4", "Ornithorhynchus_anatinus.OANA5", "Loxodonta_africana.loxAfr3", "Sus_scrofa.Sscrofa10.2", "Homo_sapiens.GRCh38", "Rattus_norvegicus.Rnor_6.0");

# All of the data about the splice positions will be printed into an intermediate
# output file for future use and inspection
open(my $intronOutFile, ">", "putativeExons.txt");
# The next line prints a header to the output file
print $intronOutFile "Species\tProtein\tHitScaffold\tPutativeGenomicExons\n";


foreach my $species (@speciesArray) {


    # First we'll put the genome scaffolds into a hash:
    my %genomeHash;
    my $genomeFile;
    if ($species eq "Homo_sapiens.GRCh38") {
        $genomeFile = "Homo_sapiens.GRCh38.dna.primary_assembly.fa";
    } else {
        my $genomeFile = $species . ".dna.toplevel.fa";
    }
    print $genomeFile . "\n";
    my $genomeIn = Bio::SeqIO->new(-file   => $genomeFile,
                                   -format => 'fasta');
    while (my $seq = $genomeIn->next_seq()) {
        $genomeHash{$seq->display_id()} = $seq;
    }
    
    # Now we'll go through the exonerate files.
    # Each reference genome has a collection of proteins associated with it, and each of these proteins have
    # been mapped to its respective reference genome. We'll iterate through each one of these exonerate output
    # files:
    my $speciesFolder = "proteins/matchingProteins/$species/exonerate";
    chdir($speciesFolder);
    opendir(my $dirFH, "./");
    my @exonerateFiles = readdir($dirFH);
    closedir($dirFH);
    
    foreach my $exonerateFile (@exonerateFiles) {
        my $exonIO = Bio::SearchIO->new(-file => $exonerateFile,
                                         -format => 'exonerate');
        while (my $result = $exonIO->next_result()) {
            print $intronOutFile $species . "\t" . $result->query_name();
            while (my $hit = $result->next_hit()) {
                print $intronOutFile "\t" . $hit->name() . "\t";
                my $numHSPS = $hit->num_hsps();
                my $hspCounter = 0;
                while (my $hsp = $hit->next_hsp()) {
                    $hspCounter++;
                    #print $intronOutFile $hsp->start('subject') . ":";
                    #print $intronOutFile $hsp->end('subject');
                    
                    # Note that filtering this way, using the BioPerl-parsed HSP start and stop
                    # base locations from an exonerate file, is a little conservative. If the
                    # bases coding for a single amino acid are split by an intron, then it
                    # chops these bases off and defines the start site as the terminal end of
                    # the last amino acid
                    print $intronOutFile $genomeHash{$hit->name()}->subseq($hsp->start('subject'), $hsp->end('subject'));
                    unless ($hspCounter == $numHSPS) {
                        print $intronOutFile ",";
                    }
                }  
            }
            print $intronOutFile "\n";
        }
    }
    chdir("../../../..");
}
```

This gets us pretty close to where we want to be. Here's an excerpt from what the file we created ("intronPositions.txt") looks like:
```text
Species	Protein	HitScaffold	PutativeGenomicExons
Strongylocentrotus_purpuratus.GCA_000002235.2.26	SPU_005493-tr	Scaffold459	CATCTCTTTGCGGAGCAAGGCCAAGGCAGCATCTAGATTTCCTCTTGGTATTTCAATGCGTTTCCTTTTTTTCCTTCCTCCACTCTTATGGACATGGCTTGGAGAAAATCCAGCAGTAAGACTGTCTTGTATATTCGTTTGGTGGGCTAAAATTCTACCTGTTTCGCGCCAAGTTCCACTGTTAGCATTTAGAAGTCCACTCAAACTTTTTGGTAACGCTGGCAGTCCTGGGAGTTTCATGTTTGAACTGGAATCCAT,GTCTTCTAGTAGTTCTTCATCAAGGGCTGATAGATGGCTGAGCTGAGAATCGTCGGCGTCCAGATCTAAAGAGCTATCAACTACCTTCTCCATATCAGACTTATAATCTTGGATGGATTCATGAAGAGAATAAAGCTCACTTAGAAGAGACATATCCAGCTGACGTAGACCAAC
Strongylocentrotus_purpuratus.GCA_000002235.2.26	SPU_003577-tr	Scaffold292	GTCATTATGTATGAAACGCTTTAGCGACCTTGGTACTAACGTGTCTCCATTTTCAGTTCCCTCTAATGTTGAGTACTCTTCACGGTTGCACACTGTCAGTCTGATATCGTTGAGGATTATGGAAGCAGCCAT,TGCCAATGCACACACACAATGGCGTGCTATACCTTTTATGGACAGGGTAAACGTATCCATTGGTTTAATTTCAGTTATGAGTTTATCTCCGTTCATATTGATGATGGGACGGCACGGGTTACTCTCCTTTTGTAAAACACCACCAGCCTTCCAGTATCCGAATACGTC,TCTGTGACACCTCGACCCTGAAAATCCGTCTGCACATGAACAAGAAGACTCCCCCGATAAGGGTACACACGTTCCGTTGTTTTCACAAGGGCTGGGGTTGCAAACTGTATC,TGTGCAGCACTTTACAGCGATGTGATATATGATACTCGTGTTTCGAAGCCTGAGTGTAGAGTACCAGCATTCCGTGAATGTGGTGTAGTGAAGCGGGCAATGGATTTCTCCAATCTCCACCAAAGCTTCCGTAGTTGGGAACATGTCACCACTCATCGTAGCAATGGCGCCACCGAAGCCCAAATCTTGACAAACCATGTGAGACATGGACATTGACCACGAGTGGCTCATGTACGCTATCGTATATGGATCATAATTGTTGCGAAACCCCACTCGTCCTTCACCGCTATGGTTGCCACCGTGTAGACTCACACTCAACAGGTAACG,CATTTCAGTTCCTTCGTCCACGTCAGTATAGATCGTCAAATCAGTTTCCTCAACACCATACGAGATCTTAAAAGAAGTCACCCACTTATCTGATCCGCCCTGTGTGACGATAGCAGTGATCATGTGATACGCCCTTAGATGGATCTGGATATAGGGTTCCGTCTGATCGTTCAAGGCAGCCATCCATGACGTCATGGCATTGAGACGGGCCTCGTGACCACACGGTGACGTAGTCGCACATGATGACGTAGTGATATCGGCGTTTAAAATATCGCCGGACTCCATTCCTAGTGGACGTCCGGATGGAGCACATTGTGA,CACGCGATCAGGTAGTGGTCCATAGCCTATGAGCTCAAATCGCATGCTGACAGTGCTGTTACTTGACTTGGGTCTTATACTAATGTACTTGGCCAGGATGTATGGAGTCAAAGATGTCGTCACGGATGTAGTATTATCGTAATTACCGGGAAAAAC,CTGCAGAAAAGGCGTGGAGTCTGTATTAAGTGGTATCCAGCCACCACCGGTAACACTGTTCAATCTAGCTGTATGAGATGGGTCACTGGGTTCGCTGGTGTGCGCTGTCAGTGATTCGTCGCCAATGTCACCATTCTCCACCCCTAGCGGAATCCCCTTATCCAAGCACGTTCCGTCTCGCTTGTGTATGTC,CTTTCTGGATCCGCATTCGTCCTCGTATAAATAAATCCATCCCTGATCAGCTTCATATGCTAGCGTGATCGACGTCGTCCATTGGTCAAGAGTACGATGACCTTGAGTTATCACACCGGTGACCAGATGTATTTCTGTAAGATCCAC,ACCTACAACATGACCTAAAATTTCACATCGAAGTGTGATCCTGTTCATCCATGATTGTGGAGAAAAGAGAACCTTTCTTGCCATCACAGGCTCAGCAAGGCGAATAGTTACAGGTGTGTTGTTGTCGAAGTTAGTTGGATAAAT
```

As you can see, we get an output file with four columns: the species/reference genome, the identifier of the protein we're mapping, the scaffold of the genome where we found a match, and a comma-separated list of genomic sequences of the putative exons which have been pulled directly from the reference genome using the coordinates given by exonerate protein2genome.

As noted in the script above, the boundaries of these putative exons are slightly conservative: if an intron splice site is inferred in the middle of an amino acid exonerate defines the edges of the exon as the first (or last) base of the next exon. For instance, if the genomic string of a match looks like this (where intron sequence is denoted with ...'s): CTT{AG}gt...ag{A}TTT, where the CTT codes for Leu, and {AG}+{A} codes for Arg, and the TTT codes for Phe, then the CTT in the first chunk forms the terminal end of one exon and the TTT in the second chunk forms the beginning end of the next exon. The end result of this is that some exons have one or two bases clipped from one or both ends, which is conservative if the objective is to designate contiguous stretches of sequence.







