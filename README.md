# Cblaster_stats
**Kindly if you find this repo useful for your work, cite & star this repo**

**What is this script?**

[Cblaster](https://github.com/gamcil/cblaster) is a great tool for finding co-located homologous sequences using BLAST. Sometimes in the remote mode, you can get the binary (~presnese absence matrix) file with thousands of hits (per isolate) depending on what you are working with. 
Here, you can break this file into a simple file per species which will count how many hits per species. 
Then, it will efetch the number of the assembled genomes per species (via Biopython [Enterez](https://biopython.org/docs/1.75/api/Bio.Entrez.html) based on the NCBI assembly database). Well, dividing both numbers will hint at the spread of your cluster among different species.
Finally, this script will draw a tree based on the pre-defined NCBI taxonomy among your species using [ete3 toolkit](http://etetoolkit.org/docs/latest/tutorial/index.html).

**What do you need?**

You shall have the binary file as easy as I get it like this.

```bash
cblaster search --query_file CP018841.1.faa --binary example_binary.csv -bde "," -bhh -bdc 6 -mi 50 -mc 50
```

 *What about dependencies?*

Pandas, Biopython, ete3, argparse

Well, for [ete3](http://etetoolkit.org/download/), I recommend installing it via conda env (even if it takes a lot of time), if the pip does not work properly.

So, just type this command effortlessly.

```bash
 python cblaster_stats.py -i  example_binary.csv
```
"-i /--input_dir"  is your path to the directory for your binary file  

**What do you get?**

Currently, there are three files.

1. database_percentage_*your_file*.csv (The main output where you can find for each species the count of this species in the binary file, the number of assembled genomes per NCBI assembly database, and the percentage of (count/assembly)*100. 
2. *your_file*_tree.png. This is just a basic tree that links your isolates together.
3. *your_file*_tree.nwk. If you like to take this tree to another tool (iTOL,FigTree,..)

I hope this helps.

Thanks
