Splicefinder
======

This is the walkthrough and code that accompanies "Conservative prediction of intron splice sites for the design of exon-capture arrays."



###Step 1: Download mouse data from Ensembl

_________________

Connect to the server, navigate to the proper folder, and download data in bash:

```bash
ssh rover
cd /mnt/Data3/arrayDesignPaper2
mkdir ensembl
cd ensembl
wget ftp://ftp.ensembl.org/pub/release-80/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz

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

That didn't nuke very many, so we'll have to reduce further. You'll notice in the fasta headers that there is a gene name in the fourth column. Let's find out how many unique gene names there are total:

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
#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
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
#!/usr/bin/bash
#downloadGenomes.bash
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

We'll also need the full collection of cdna sequences from these organisms:
```bash
cd /mnt/Data3/arrayDesignPaper/ensembl/genomes

#!/bin/bash
#downloadCDNAs.bash
mkdir cdnas
cd cdnas

wget ftp://ftp.ensembl.org/pub/release-80/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/sus_scrofa/cdna/Sus_scrofa.Sscrofa10.2.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/loxodonta_africana/cdna/Loxodonta_africana.loxAfr3.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/ornithorhynchus_anatinus/cdna/Ornithorhynchus_anatinus.OANA5.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/gallus_gallus/cdna/Gallus_gallus.Galgal4.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/xenopus_tropicalis/cdna/Xenopus_tropicalis.JGI_4.2.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/anolis_carolinensis/cdna/Anolis_carolinensis.AnoCar2.0.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/latimeria_chalumnae/cdna/Latimeria_chalumnae.LatCha1.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/takifugu_rubripes/cdna/Takifugu_rubripes.FUGU4.cdna.all.fa.gz
wget ftp://ftp.ensembl.org/pub/release-80/fasta/petromyzon_marinus/cdna/Petromyzon_marinus.Pmarinus_7.0.cdna.all.fa.gz
wget ftp://ftp.ensemblgenomes.org/pub/metazoa/release-26/fasta/strongylocentrotus_purpuratus/cdna/Strongylocentrotus_purpuratus.GCA_000002235.2.26.cdna.all.fa.gz

gunzip *
```

We're going to use exonerate est2genome to map all of these proteins to their respective genomes.
```bash
cd /mnt/Data3/arrayDesignPaper/ensembl/genomes
```
```perl
#!/usr/bin/perl
#prepareExonerate.pl

use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;


my @fasta2esdCommands = ("fasta2esd --softmask no Loxodonta_africana.loxAfr3.dna.toplevel.fa Loxodonta_africana.loxAfr3.dna.toplevel.esd",
                         "fasta2esd --softmask no Ornithorhynchus_anatinus.OANA5.dna.toplevel.fa Ornithorhynchus_anatinus.OANA5.dna.toplevel.esd",
                         "fasta2esd --softmask no Strongylocentrotus_purpuratus.GCA_000002235.2.26.dna.toplevel.fa Strongylocentrotus_purpuratus.GCA_000002235.2.26.dna.toplevel.esd",
                         "fasta2esd --softmask no Petromyzon_marinus.Pmarinus_7.0.dna.toplevel.fa Petromyzon_marinus.Pmarinus_7.0.dna.toplevel.esd",
                         "fasta2esd --softmask no Takifugu_rubripes.FUGU4.dna.toplevel.fa Takifugu_rubripes.FUGU4.dna.toplevel.esd",
                         "fasta2esd --softmask no Latimeria_chalumnae.LatCha1.dna.toplevel.fa Latimeria_chalumnae.LatCha1.dna.toplevel.esd",
                         "fasta2esd --softmask no Anolis_carolinensis.AnoCar2.0.dna.toplevel.fa Anolis_carolinensis.AnoCar2.0.dna.toplevel.esd",
                         "fasta2esd --softmask no Xenopus_tropicalis.JGI_4.2.dna.toplevel.fa Xenopus_tropicalis.JGI_4.2.dna.toplevel.esd",
                         "fasta2esd --softmask no Gallus_gallus.Galgal4.dna.toplevel.fa Gallus_gallus.Galgal4.dna.toplevel.esd",
                         "fasta2esd --softmask no Sus_scrofa.Sscrofa10.2.dna.toplevel.fa Sus_scrofa.Sscrofa10.2.dna.toplevel.esd",
                         "fasta2esd --softmask no Homo_sapiens.GRCh38.dna.primary_assembly.fa Homo_sapiens.GRCh38.dna.primary_assembly.esd",
                         "fasta2esd --softmask no Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa Rattus_norvegicus.Rnor_6.0.dna.toplevel.esd");

my @esd2esiCommands = ("esd2esi Loxodonta_africana.loxAfr3.dna.toplevel.esd Loxodonta_africana.loxAfr3.dna.toplevel.esi --memorylimit 25600",
                       "esd2esi Ornithorhynchus_anatinus.OANA5.dna.toplevel.esd Ornithorhynchus_anatinus.OANA5.dna.toplevel.esi --memorylimit 25600",
                       "esd2esi Strongylocentrotus_purpuratus.GCA_000002235.2.26.dna.toplevel.esd Strongylocentrotus_purpuratus.GCA_000002235.2.26.dna.toplevel.esi --memorylimit 25600",
                       "esd2esi Petromyzon_marinus.Pmarinus_7.0.dna.toplevel.esd Petromyzon_marinus.Pmarinus_7.0.dna.toplevel.esi --memorylimit 25600",
                       "esd2esi Takifugu_rubripes.FUGU4.dna.toplevel.esd Takifugu_rubripes.FUGU4.dna.toplevel.esi --memorylimit 25600",
                       "esd2esi Latimeria_chalumnae.LatCha1.dna.toplevel.esd Latimeria_chalumnae.LatCha1.dna.toplevel.esi --memorylimit 25600",
                       "esd2esi Anolis_carolinensis.AnoCar2.0.dna.toplevel.esd Anolis_carolinensis.AnoCar2.0.dna.toplevel.esi --memorylimit 25600",
                       "esd2esi Xenopus_tropicalis.JGI_4.2.dna.toplevel.esd Xenopus_tropicalis.JGI_4.2.dna.toplevel.esi --memorylimit 25600",
                       "esd2esi Gallus_gallus.Galgal4.dna.toplevel.esd Gallus_gallus.Galgal4.dna.toplevel.esi --memorylimit 25600",
                       "esd2esi Sus_scrofa.Sscrofa10.2.dna.toplevel.esd Sus_scrofa.Sscrofa10.2.dna.toplevel.esi --memorylimit 25600",
                       "esd2esi Homo_sapiens.GRCh38.dna.primary_assembly.esd Homo_sapiens.GRCh38.dna.primary_assembly.esi --memorylimit 25600",
                       "esd2esi Rattus_norvegicus.Rnor_6.0.dna.toplevel.esd Rattus_norvegicus.Rnor_6.0.dna.toplevel.esi --memorylimit 25600");

my $forkManager = Parallel::ForkManager->new(8);
foreach my $command (@fasta2esdCommands) {
    $forkManager->start and next;
    system($command);
    $forkManager->finish;
}
$forkManager->wait_all_children;

my $forkManager 2= Parallel::ForkManager->new(8);
foreach my $command (@esd2esiCommands) {
    $forkManager2->start and next;
    system($command);
    $forkManager2->finish;
}
$forkManager2->wait_all_children;
```

The above will take several hours. We want to map the orthologous proteins from each species to their own reference genomes. So, we will need to find the orthologous cDNAs from each species. First we'll make blast databases from each cDNA set:
```bash
cd cdnas
for i in *.cdna.all.fa; do makeblastdb -in $i -dbtype nucl; done
```

Now we'll iterate through each of the cDNA sets and tblastx the mouse cDNAs to each, only outputting the best match:


```perl
#!/usr/bin/perl

#runTBlastX.pl
use strict;
use warnings;

my @speciesArray = ("Strongylocentrotus_purpuratus.GCA_000002235.2.26", "Petromyzon_marinus.Pmarinus_7.0", "Takifugu_rubripes.FUGU4", "Latimeria_chalumnae.LatCha1", "Anolis_carolinensis.AnoCar2.0", "Xenopus_tropicalis.JGI_4.2", "Gallus_gallus.Galgal4", "Ornithorhynchus_anatinus.OANA5", "Loxodonta_africana.loxAfr3", "Sus_scrofa.Sscrofa10.2", "Homo_sapiens.GRCh38", "Rattus_norvegicus.Rnor_6.0");

foreach my $species (@speciesArray) {
    my $speciesCDNAfile = $species . ".cdna.all.fa";
    my $outFile = $species . ".cdna.all.fa.cdnaMatches";
    system("tblastx -db $speciesCDNAfile -query ../../mmGRCm38.cdna.rand10kLongest.fa -outfmt 6 -max_target_seqs 1 -num_threads 32 -out $outFile);
}
```


The blastx'ing will also take a few hours.

After the blastx is done, we'll parse the blastx results to pull out all the cDNA sequences that matched and store them into a new file for exonerate est2genome mapping:

```perl
#!/usr/bin/perl

# pullOutMatchingCDNAs.pl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::SearchIO;

my @speciesArray = ("Strongylocentrotus_purpuratus.GCA_000002235.2.26", "Petromyzon_marinus.Pmarinus_7.0", "Takifugu_rubripes.FUGU4", "Latimeria_chalumnae.LatCha1", "Anolis_carolinensis.AnoCar2.0", "Xenopus_tropicalis.JGI_4.2", "Gallus_gallus.Galgal4", "Ornithorhynchus_anatinus.OANA5", "Loxodonta_africana.loxAfr3", "Sus_scrofa.Sscrofa10.2", "Homo_sapiens.GRCh38", "Rattus_norvegicus.Rnor_6.0");

foreach my $species (@speciesArray) {
    my $cdnaFasta = $species . ".cdna.all.fa";

    my %cdnaHash;
    my $seqIn = Bio::SeqIO->new(-file => $cdnaFasta,
                                -format => 'fasta');
    while (my $seq = $seqIn->next_seq()) {
        $cdnaHash{$seq->display_id()} = $seq;
    }

    my $cdnaOut = $species . ".cdna.matching.fa";
    my $seqOut = Bio::SeqIO->new(-file => ">$cdnaOut",
                                 -format => 'fasta');

    my $blastResults = $species . ".cdna.all.fa.cdnaMatches";
    my %keeperHash;
    open(my $blastFH, "<", $blastResults) or die "Couldn't open $blastResults for reading: $!\n";
    while (my $line = <$blastFH>) {
        my @fields = split(/\t/, $line);
        $seqOut->write_seq($cdnaHash{$fields[1]});
    }
}
```

Now we'll set these aside:

```bash
mkdir matchingCDNAs
mv "*.cdna.matching.fa" matchingCDNAs/
```

We also need to create a map of which cDNAs match the query transcripts from
the different reference genomes. We can do that by parsing the blast results with
this script, which outputs a file called "blastMap.txt" with that info:

```perl
#!/usr/bin/perl

# makeBlastMap.pl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::SearchIO;

my @speciesArray = ("Strongylocentrotus_purpuratus.GCA_000002235.2.26", "Petromyzon_marinus.Pmarinus_7.0", "Takifugu_rubripes.FUGU4", "Latimeria_chalumnae.LatCha1", "Anolis_carolinensis.AnoCar2.0", "Xenopus_tropicalis.JGI_4.2", "Gallus_gallus.Galgal4", "Ornithorhynchus_anatinus.OANA5", "Loxodonta_africana.loxAfr3", "Sus_scrofa.Sscrofa10.2", "Homo_sapiens.GRCh38", "Rattus_norvegicus.Rnor_6.0");

open(my $outFH, ">", "blastMap.txt");
my %resultsHash;
foreach my $species (@speciesArray) {
    my $blastResults = $species . ".cdna.all.fa.cdnaMatches";
    open(my $blastFH, "<", $blastResults) or die "Couldn't open $blastResults for reading: $!\n";
    my $result = " ";
    while (my $line = <$blastFH>) {
        my @fields = split(/\t/, $line);
        unless ($result eq "$species\t$fields[0]\t$fields[1]\n") {
    	    print $outFH $species . "\t" . $fields[0] . "\t" . $fields[1] . "\n";
    	    $result = "$species\t$fields[0]\t$fields[1]\n";
    	}
    }
}
```





Exonerate est2genome is very slow unless you change some settings around, and it seems to run faster when you separate each sequence into its own file. Using only one sequence per run also has the benefit of showing you where you left off if exonerate segfaults.

So we're first gonna need to separate the protein fasta files into new files, one sequence per file:

```perl
#!/usr/bin/perl

# separateCDNAsIntoIndividuals.pl

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

    my $seqsFile = "../" . $species . ".cdna.matching.fa";

    my $seqIn = Bio::SeqIO->new(-file => $seqsFile,
                                -format => 'fasta');

    my $counter = 0;
    while (my $seq = $seqIn->next_seq()) {
        $counter++;
        my $seqOutName = $seq->display_id() . ".cdna.fasta";
        my $seqOut = Bio::SeqIO->new(-file => ">$seqOutName",
                                     -format => 'fasta');

        $seqOut->write_seq($seq);

    }
    print $counter . " total records processed\n";
    chdir "..";
}
```

Note that some mouse transcripts might have selected the same genes from the reference set as others for the best blastx hit. This means that there should be fewer actual cDNA fastas than query sequences that had positive matches in the blastx search.


#### Finding splice sites in reference cDNAs ####
Now let's loop through all those sequences and align them to their proper genomes. The sending multiple queries to a single exonerate-server tends to make things segfault, so instead we'll instantiate exonerate-servers so that each thread gets its own.


```perl
#!/usr/bin/perl

# exonerateE2G.pl

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
        $genomeFile = "Homo_sapiens.GRCh38.dna.primary_assembly.esi";
    } else {
        $genomeFile = $species . ".dna.toplevel.esi";
    }


    my $forkManager = Parallel::ForkManager->new(8);

    foreach my $thread (0..7) {
        $portCounter++;
        $forkManager->start and next;
        system("exonerate-server --port $portCounter $genomeFile &");
        sleep(90);
        chdir("cdnas/matchingCDNAs/");
        chdir $species;
        unless(-d "exonerate") {
            mkdir "exonerate";
        }

        # Read in all of the cDNA file names
        opendir(my $dirFH, "./");
        my @files = readdir($dirFH);
        closedir($dirFH);

        my @fastaFiles;
        foreach my $file (@files) {
            if ($file =~ /.fasta/) {
                push(@fastaFiles, $file);
            }
        }

        # TODO: Change the hard coded number for threads below to a variable for number of threads
        my $cDNAsPerThread = scalar(@fastaFiles) / 8; # Decide how many cDNAs to provide to each thread

        my $beginning = $thread * $cDNAsPerThread;
        my $end = ($thread+1) * $cDNAsPerThread;
        for (my $fileIndex=$beginning; $fileIndex < $end; $fileIndex++) {
            my $outFileName = "exonerate/" . $fastaFiles[$fileIndex] . ".e2g.exonerate";
            system("exonerate $fastaFiles[$fileIndex] localhost:$portCounter --querytype dna --targettype dna --model est2genome --bestn 1 --dnawordlen 12 --fsmmemory 4096 --hspfilter 50 --geneseed 250 > $outFileName");
        }
        $forkManager->finish;
    }
    $forkManager->wait_all_children;
    sleep(120);
    system("pkill -f exonerate-server");
    sleep(120);
}
```

This takes several hours to a day to run all these, assuming all goes well and you have > 8 cores or so to work with.





#### Harvesting splice sites and mapping to ESTs ####

After this comes one of the trickier parts. We'll process all of these exonerate files to harvest the
putative contiguous genomic region sequences. For instance, rat cDNA ENSRNOT00000018050 was mapped to
rate chromosome 4 between base pairs 148761547 and 148757362. Additionally, exonerate found 1 intron
between bases 206 and 207 of the cDNA. We want to pull the rat DNA sequence from 1-206 in the cDNA (and also the one from 207-2378) and align it to the original mouse EST sequence.

Let's say this exon sequence aligns at base pairs 10-215 in the mouse EST. We would then code this as rat-inferred splice sites before base 10 and after 215 in the mouse EST. We'll gather all such inferred splice sites for each EST, keeping track of which ones come from which species.

I'll present this all in heavily-commented code below.

```perl
#!/usr/bin/perl

# harvestSplicePositions.pl

use strict;
use warnings;
use Bio::SearchIO; # <= We'll use this to parse the exonerate reports
use Bio::SeqIO;


# We have information on where a species' cDNA mapped to its own genome, but we also
# need to know which original species cDNA that cDNA corresponded to. Here we'll build both
# a hash of all of the original cDNA sequences and also an index of all the original mRNA
# query sequences:
my %mRNAhash;
###TODO: Change the filename in the next line to a variable set in the config file
my $mRNAin = Bio::SeqIO -> new(-file => "../mmGRCm38.cdna.rand10kLongest.fa",
                               -format => 'fasta');
while (my $seq = $mRNAin->next_seq()) {
    $mRNAhash{$seq->display_id()} = $seq;
}

# Here's the map of which cDNAs from the different reference genomes
# correspond to the original mRNA query sequences. We created this earlier in a file called
# blastMap.txt, which is in the cdnas folder:
my %blastMap;
open(my $blastMapFH, "<", "cdnas/blastMap.txt") or die "Couldn't open blastMap.txt: $!\n";
while(my $line = <$blastMapFH>) {
    chomp($line);
    my @fields = split(/\t/, $line);
    $blastMap{$fields[0]}{$fields[2]} = $fields[1]; # This equates to $blastMap{species}{speciesCDNAName} = queryRNAname
}



# We'll deal with all the different reference genomes one-by-one in a loop
my @speciesArray = ("Strongylocentrotus_purpuratus.GCA_000002235.2.26", "Petromyzon_marinus.Pmarinus_7.0", "Takifugu_rubripes.FUGU4", "Latimeria_chalumnae.LatCha1", "Anolis_carolinensis.AnoCar2.0", "Xenopus_tropicalis.JGI_4.2", "Gallus_gallus.Galgal4", "Ornithorhynchus_anatinus.OANA5", "Loxodonta_africana.loxAfr3", "Sus_scrofa.Sscrofa10.2", "Homo_sapiens.GRCh38", "Rattus_norvegicus.Rnor_6.0");

# All of the data about the splice positions will be printed into an intermediate
# output file for future use and inspection
open(my $intronOutFile, ">", "putativeSplicePositions.txt");
# The next line prints a header to the output file
print $intronOutFile "targetCDNA\tSpecies\tPutativeSplicePositions\n";

foreach my $species (@speciesArray) {
    # First we'll put the genome scaffolds into a hash:
    my %genomeHash;
    my $genomeFile;
    if ($species eq "Homo_sapiens.GRCh38") {
        $genomeFile = "Homo_sapiens.GRCh38.dna.primary_assembly.fa";
    } else {
        $genomeFile = $species . ".dna.toplevel.fa";
    }
#    print $genomeFile . "\n";
    my $genomeIn = Bio::SeqIO->new(-file   => $genomeFile,
                                   -format => 'fasta');
    while (my $seq = $genomeIn->next_seq()) {
        $genomeHash{$seq->display_id()} = $seq;
    }

    my $speciesFolder = "cdnas/matchingCDNAs/$species/exonerate";
    chdir($speciesFolder);
    opendir(my $dirFH, "./");
    my @exonerateFiles = readdir($dirFH);
    closedir($dirFH);
    
    foreach my $exonerateFile (@exonerateFiles) {
        unless ($exonerateFile =~ /\.exonerate$/) {
	    next;
	}
 #       print $exonerateFile . "\n";
        my $exonIO = Bio::SearchIO->new(-file => $exonerateFile,
                                         -format => 'exonerate');
        while (my $result = $exonIO->next_result()) {
            while (my $hit = $result->next_hit()) {
                my $numHSPS = $hit->num_hsps();
                my $hspCounter = 0;
                while (my $hsp = $hit->next_hsp()) {

                    # Pull out the genomic sequence from the original genome file that has the same bounds as the HSP
                    # Then we'll write that HSP sequence to a new file
                    my $sequence = $genomeHash{$hit->name()}->subseq($hsp->start('subject'), $hsp->end('subject'));
                    print $sequence . "\n";
                    my $seqFileName = $hit->name() . "." . $hsp->start . "-" . $hsp->end . ".fasta";
                    open(my $exonFH, ">", $seqFileName);
                    print $exonFH ">$seqFileName\n$sequence\n";
                    close($exonFH);
                    
                    
                    # We need to pull out the original cDNA sequence from the original possible target set,
                    # and write it to a temporary file
                    my $originalCDNAname = $mRNAhash{$blastMap{$species}{$result->query_name}}->display_id();
                    print "originalCDNAname = $originalCDNAname\n";
                    my $originalCDNAfileName = "originalSeq." . $originalCDNAname . ".fasta";
                    open(my $tempSeqFH, ">", $originalCDNAfileName);
                    print $tempSeqFH ">" . $mRNAhash{$blastMap{$species}{$result->query_name}}->display_id() . "\n" . $mRNAhash{$blastMap{$species}{$result->query_name}}->seq() . "\n";

                    

                    my $mapExonFile = $result->query_name . "." . $hsp->start . "-" . $hsp->end . ".map2." . $originalCDNAname . ".exonerate";

                    system("exonerate --model coding2coding --bestn 1 $seqFileName $originalCDNAfileName > $mapExonFile");
                    #system("tblastx -query $seqFileName -db $originalCDNAfileName > $mapExonFile");
                    
                    
                    # Now we need to parse that exonerate output
                    
                    
                    
                    # Now we delete the exonerate output
                    
                    


                    # Note that filtering this way, using the BioPerl-parsed HSP start and stop
                    # base locations from an exonerate file, is a little conservative. If the
                    # bases coding for a single amino acid are split by an intron, then it
                    # chops these bases off and defines the start site as the terminal end of
                    # the last amino acid
                    unless ($hspCounter == $numHSPS) {
                        #print $intronOutFile ",";
                    }
                }  
            }
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

As you can see, we get an output file with four columns: the species/reference genome, the
identifier of the protein we're mapping, the scaffold of the genome where we found a match, and a
comma-separated list of genomic sequences of the putative exons which have been pulled directly
from the reference genome using the coordinates given by exonerate protein2genome.

As noted in the script above, the boundaries of these putative exons are slightly conservative: if
an intron splice site is inferred in the middle of an amino acid exonerate defines the edges of
the exon as the first (or last) base of the next exon. For instance, if the genomic string of a
match looks like this (where intron sequence is denoted with ...'s): CTT{AG}gt...ag{A}TTT, where
the CTT codes for Leu, and {AG}+{A} codes for Arg, and the TTT codes for Phe, then the CTT in the
first chunk forms the terminal end of one exon and the TTT in the second chunk forms the beginning
end of the next exon. The end result of this is that some exons have one or two bases clipped from
one or both ends, which is conservative if the objective is to designate contiguous stretches of
sequence.

Now we want to take those putative contiguous stretches of DNA in the fourth column of the output
file and align them back to the original mRNA sequences that we started with (the sequenced
transcriptome from our non-model species that we want to use to design capture targets, for
instance). We'll do that with the following script:

```perl
#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;
use Bio::SearchIO;
use Data::Dumper;

# We'll need to store the locations of all of our inferred splice sites:
open(my $spliceSitesFile, ">", "inferredSpliceSitesEST2Genome.txt");


# We also need to have an index of all the original mRNA query sequences:
my %mRNAhash;
my $mRNAin = Bio::SeqIO -> new(-file => "../mmGRCm38.cdna.rand10kLongest.fa",
                               -format => 'fasta');
while (my $seq = $mRNAin->next_seq()) {
    $mRNAhash{$seq->display_id()} = $seq;
}


# Finally, we need a map of which proteins from the different reference genomes
# correspond to the mRNA query sequence. We created this earlier in a file called
# blastMap.txt:
my %blastMap;
open(my $blastMapFH, "<", "/mnt/Data3/arrayDesignPaper/ensembl/genomes/proteins/blastMap.txt") or die "Couldn't open blastMap.txt: $!\n";
while(my $line = <$blastMapFH>) {
    chomp($line);
    my @fields = split(/\t/, $line);
    $blastMap{$fields[0]}{$fields[2]} = $fields[1]; # This equates to $blastMap{species}{speciesProteinName} = queryRNAname
}


# We'll be going line-by-line through the following file to process the putative exons
open(my $exonsFile, "<", "putativeExons.txt");

my $currentGenome; # We'll hold the current genome we're working on in a variable here
                   # so that we don't have to load it into memory on every iteration, just
                   # when we change genomes
my %genomeHash;
               
               


while (my $line = <$exonsFile>) {
    chomp($line);
    next if ($line =~ /^Species\tProtein\tHitScaffold\tProteinIntronSites/); # Skip the header line
    my @fields = split(/\t/, $line);
    # $fields[0] is the species name, $fields[1] is the Protein ID, $fields[2] is
    # the genomic sequence, and $fields[3] are the comma-delimited putative exons
    
    
    ### Make sure we have the right genome indexed
    if (!$currentGenome or $currentGenome ne $fields[0]) {
        %genomeHash = (); # empty out the hash
        my $genomeFile;
        if ($fields[0] eq "Homo_sapiens.GRCh38") {
           $genomeFile = "Homo_sapiens.GRCh38.dna.primary_assembly.fa";
        } else {
            $genomeFile = $fields[0] . ".dna.toplevel.fa";
        }
        my $genomeIn = Bio::SeqIO->new(-file   => $genomeFile,
                                       -format => 'fasta');
        while (my $seq = $genomeIn->next_seq()) {
            $genomeHash{$seq->display_id()} = $seq;
        }    
    }
    
    
    # Now let's pull out all of the strings that represent putative exons for that target
    my @putativeRefExons = split(/,/, $fields[3]);

    foreach my $putativeRefExon (@putativeRefExons) {
        my $proteinFileName = $fields[1] . ".fasta";
        open(my $tempRefFasta, ">", "$proteinFileName");
        print $tempRefFasta ">tempRefExon\n";
        print $tempRefFasta $putativeRefExon . "\n";
        close($tempRefFasta);
        
        # We now have the sequence that we want in tempRefExon.fasta, which
        # will be rewritten in every iteration of the loop. We want to align
        # that to a single target transcript, which we'll pull from the blast map
        my $name = $blastMap{$fields[0]}{$fields[1]} . ".fasta";
        open(my $tempMRNAfasta, ">", $name);
        print $tempMRNAfasta ">$blastMap{$fields[0]}{$fields[1]}\n" . $mRNAhash{$blastMap{$fields[0]}{$fields[1]}}->seq ."\n";
        close($tempMRNAfasta);
        
        # Now align those two sequences:
        #system("exonerate -m est2genome --bestn 1 -q $proteinFileName -t $name > exonerate.output");
        
        # Or maybe use the ungapped:trans? I used est2genome in the salamander array design
        # my $exonerateOutput = `exonerate -m ungapped:trans --bestn 1 -q $proteinFileName -t $name --showtargetgff 1`;
        # my @exonerateOutputLines = split(/\n/, $exonerateOutput);
        # 
        # my $lineCounter = 0;
        # foreach my $exonerateLine (@exonerateOutputLines) {
        #     $lineCounter++;
        #     if ($exonerateLine =~ /^# seqname source feature start end score strand frame attributes/) {
        #         my $GFFline = $exonerateOutputLines[$lineCounter+2]; # This should be the actual GFF line
        #         my @GFFfields = split(/\t/, $GFFline);
        #         print $spliceSitesFile $blastMap{$fields[0]}{$fields[1]} . "\t" . $fields[0] . "\t" . $fields[1] . "\t" . $GFFfields[3] . "," . $GFFfields[4] . "\n";
        #         
        #     }
        # }
        
        
        # ORRRRRR should I maybe use est2genome, but switch the query and target, like so:
                my $exonerateOutput = `exonerate -m est2genome --bestn 1 -q $name -t $proteinFileName --showtargetgff 1`;
        my @exonerateOutputLines = split(/\n/, $exonerateOutput);
        
        my $lineCounter = 0;
        foreach my $exonerateLine (@exonerateOutputLines) {
            $lineCounter++;
            if ($exonerateLine =~ /^# seqname source feature start end score strand frame attributes/) {
                my $GFFline = $exonerateOutputLines[$lineCounter+2]; # This should be the actual GFF line
                my @GFFfields = split(/\t/, $GFFline);
                print $spliceSitesFile $blastMap{$fields[0]}{$fields[1]} . "\t" . $fields[0] . "\t" . $fields[1] . "\t" . $GFFfields[3] . "," . $GFFfields[4] . "\n";
                
            }
        }





#        # And now we'll process the output to define the splice sites
#        my $exonIn = Bio::SearchIO->new(-file => "exonerate.output",
#                                        -format => 'exonerate');
#        while (my $result = $exonIn->next_result()) {
#            while (my $hit = $result->next_hit()) {
#                # If there are multiple HSPs in our file, then that means that
#                # there is a putative intron introduced. We don't expect that to
#                # happen, so we'll at least keep track of how many times it does...
#                
#                if ($hit->num_hsps() > 1) { print "More than one HSP for a putative exon in $fields[1]\n"; }
#                while (my $hsp = $hit->next_hsp()) {
#                    my $start = $hsp->start('subject');
#                    my $end = $hsp->end('subject');
#                    
#                    if ($start < $end) {
#                        $start = $start - 0.5;
#                        $end = $end + 0.5;
#                    } elsif ($start > $end) {
#                        $start = $start + 0.5;
#                        $end = $end - 0.5;
#                    } else {
#                        die "Contiguous region start site isn't larger or smaller than end site\n";
#                    }
#                    # We want to print the following to the output file: MouseProteinID  Species  InferredSpliceSites
#                    print $spliceSitesFile $blastMap{$fields[0]}{$fields[1]} . "\t" . $fields[0] . "\t" . $fields[1] . "\t" . $start . "," . $end . "\n";
#                }
#            }
#    
#        }                
#        # And finally we'll get rid of all those temporary files:
        unlink($name, $proteinFileName);






    }
}

```










### TODO ###
    --Make it easier to download the ENSEMBL genomes and cDNAs
    --Use a config file to specify locations of genome and cDNA files
    --Run first with 100 cDNAs, then change to 10,000
    --Make explicit that we want exonerate 2.4
    --Make sure the soft-masking can be set in the config file



### Config file should include ###
    1) How many threads
    2) Which genomes
    3) Where genomes are
    4) How many loci to include
    5) Minimum target length
        5a) Maximum target length--distribution (beta? uniform?)
    6) Path to tblastx
    7) Ports for exonerate-server
    8) Path to exonerate est2genome
    9) Output directory name
    10) Name of cDNA fasta for original targets
    11) Change the pullOutMatchingCDNAs.pl to not output duplicates
    12) Sort out logic of multiple cDNA queries matching the same target cDNAs
    13) Soft-masked genome files?







