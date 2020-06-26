# besolverd.py

As promised, I made a small python (2.7) script which encapsulates
mosdepth, bedtools and rtg tools for our benchmarking.

## STEPS OF THE SCRIPT

The tool first runs mosdepth to compute the global coverage per base.

 From this file, it will then use bedtools to find the regions covered
at a depth equal or higher than 0, 5, 10, 20, 30, 42, 50 and 100 and
intersect those with the "high confidence regions" that are given by the
GIAB for both reference samples. It then runs the comparaison on those
region with rtg-tool.

At the very end, it creates a small tab delimited file containing the
summary results of all iterations.

I made the the script (besolved.py) available via dockerhub. This image
(which is a bit heavy) does contain all the reference files for hg37 and
hg38. Provided you run it within the docker image, all you have to
provide are

- VCF file delivered by your pipeline
- BAM file that was created by your pipeline
- Genome fasta file used for the alignment

The chromosomes must be encoded with the chr prefix before the
chromosomes id. If it is a real problem, I can of course figure out a
tweak around it. Hence, let me know.

The procedure is quite slow as it deals with big files and there are
quite a few iterations (3h per benchmark). It would be of interest to
parallelize the script, but I haven't found the courage yet for that. I
also attach the source code of the script to this email.


## INSTALLATION
```
docker pull sbrohee/besolverd
```

## RUNNING (with docker)


### Get the options
```
docker run  -it besolverd python2.7 besolverd.py -h
```

### Run in your home directory being the same user (not as root) having the same mount points within the container
```
docker run -u $(id -u ${USER}):$(id -g ${USER}) -v
/my/home/absolute/path/:/my/home/absolute/path/:rw besolverd.py -h
```

### Example run on cell line N124385 with hg38 genome build

```
docker run -u $(id -u ${USER}):$(id -g ${USER}) -v
/data/home:/data/home:rw -it besolverd python2.7 besolverd.py -q
INPUTVCF  -b INPUTBAM  -g hg38 -o /outputdir/outputprefix -f INPUTFASTA
-u NA24385
```
## SOME REMARKS

It is of course not mandatory to use the tool, nor the procedure but at
least it can be used as a starting point strategy to evaluate our pipelines.

I have been the only developer of the script so it might contain error.
I won't be offended if you ask that I adapt, correct or add something to
the script. Of course, you can also modify the script yourself (even if
it is never a pleasure to dig into the awful python code of someone else).
