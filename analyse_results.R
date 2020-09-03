library(optparse)
library(ggplot2)
library(data.table)
library(dplyr) 

# Example command 
#Rscript analyse_results.R -i results_sub/ -o figures/ -c ULiege,IPG,KULeuven,UGent,VUB -l NA24385,NA12878 -d 42x,30x,25x -k Truseq,Kapa,NexteraFlex



#########################
## Define the arguments
option_list = list(
    make_option(c("-i", "--inputDir"), 
    	type="character",
    	default=NA, 
    	help="Directory containing the directories with all the final results files (prefix _results.tab)", 
    	metavar="character"),
    make_option(c("-o", "--outputDir"), 
    	type="character",
    	default=NA, 
    	help="Directory were the results will be written", 
    	metavar="character"),    	
    make_option(c("-c", "--centerMatch"),
    	type="character",
    	default=NA, 
    	help="Strings to be found in the file name to obtain the genetic center id. The different centers must be separated by a comma. e.g. -c ULiege,IPG,KULeuven,UGent,VUB", 
    	metavar="character"),
    make_option(c("-l", "--cellineMatch"),
    	type="character",
    	default=NA, 
    	help="Strings to be found in the file name to obtain the cell line. The different cell lines must be separated by a comma. e.g. -l NA24385,NA12878", 
    	metavar="character"),  
    make_option(c("-d", "--coverageMatch"),
    	type="character",
    	default=NA, 
    	help="Strings to be found in the file name to obtain the coverage. The strings must be comma separated. e.g. -d 42x,30x,25x ", 
    	metavar="character"),
    make_option(c("-k", "--kitMatch"),
    	type="character",
    	default=NA, 
    	help="Strings to be found in the file name to obtain the kit name The strings must be comma separated. e.g. -d Truseq,Kapa,NexteraFlex ", 
    	metavar="character")    	
    
    
); 



opt_parser = OptionParser(usage = "%prog [options] ", option_list=option_list);
arguments = parse_args(opt_parser)
opt = arguments
inputDir <- opt$inputDir
outputDir <- opt$outputDir
centerMatch <- opt$centerMatch
coverageMatch <- opt$coverageMatch
cellineMatch <- opt$cellineMatch
kitMatch <- opt$kitMatch
if (is.na(inputDir)) {
    stop("Error : missing argument -i")
}
if (is.na(centerMatch)) {
    stop("Error : missing argument -c")
}
if (is.na(coverageMatch)) {
    stop("Error : missing argument -d")
}
if (is.na(kitMatch)) {
    stop("Error : missing argument -k")
}
if (!file.exists(inputDir)) {
    stop(paste("Directory" , inputDir, "does not seem to exist. Please check"));
}
if (!file.exists(outputDir)) {
    stop(paste("Directory" , outputDir, "does not seem to exist. Please check"));
}

results.files <- list.files(path = opt$inputDir, pattern = "_results.tab$", full.names = T, recursive= T)
coverage.strings <- strsplit(coverageMatch, ",")[[1]]
center.strings <- strsplit(centerMatch, ",")[[1]]
cell.lines.strings <- strsplit(cellineMatch, ",")[[1]]
kits.strings <- strsplit(kitMatch, ",")[[1]]
table.list <- list()


if (length(results.files) == 0) {
    stop(paste("Directory",inputDir, "does not contain any result file produced by besolverd.py (suffix : _results.tab"))
}

for (result.file in results.files) {
    #parse name
    coverage <- NA
    center <- NA
    cell.line <- NA
    for (coverage.string in coverage.strings) {
        if (grepl(coverage.string, result.file)) {
            coverage <- coverage.string
        }
    }
    for (center.string in center.strings) {
        if (grepl(center.string, result.file)) {
            center <- center.string
        }
    }
    for (cell.lines.string in cell.lines.strings) {
        if (grepl(cell.lines.string, result.file)) {
            cell.line <- cell.lines.string
        }
    }
    for (kit.string in kits.strings) {
        if (grepl(kit.string, result.file)) {
            kit <- kit.string
        }
    }    
    if (is.na(center)) {
        stop(paste("Could not find center for file", result.file))
    } 
    if (is.na(coverage)) {
        stop(paste("Could not find coverage for file", result.file))
    }     
    results <- fread(paste("cat",result.file ,"| perl -pe 's/\t/ /g'", sep = " "))
    results$coverage <- coverage
    
    results$center <- center
    results$cell.line <- cell.line
    results$kit <- kit
    names(results)[1] <- "mincov"
    table.list[[length(table.list)+1]] <- results
}

table.dt <- rbindlist(table.list)
# print(table.dt)
table.dt$centerkit <- apply(table.dt[,.(center, kit)],1,paste, collapse = "_")



cell.lines <- unique(table.dt$cell.line)
mycols <- c("darkolivegreen1", "cadetblue1")

for (cell.linei in cell.lines) {

    table.cl <- table.dt[cell.line == cell.linei]
    mincovs <- unique(table.cl$mincov)
    for (mincovi in mincovs) {
        table.cl.mc <- table.cl[mincov == mincovi]
        coverages <- unique(table.cl.mc$coverage)
        file.name <- file.path(outputDir, paste0("besolverd_snprecision_", cell.linei, "_mincov", mincovi, ".pdf"))
        pdf( file.name, width = 6, height = 9)
        par(mfrow = c(3,1))
        for (coveragei in coverages) {
            table.cl.mc.cov <- table.cl.mc[coverage == coveragei]
            centers <- unique(table.cl.mc.cov$centerkit)

            mat <- t(table.cl.mc.cov[,.(Sensitivity, Precision)])
            colnames(mat) <- table.cl.mc.cov$centerkit
            barplot(mat, beside = T, ylim = c(0,1.2), col = mycols, cex.names = 0.7,panel.first = grid(col = "black",nx = NA, ny = NULL))
            legend ("topright",legend = c("Sensitivity", "Precision"), fill = mycols, cex = 0.55, bty = "n")
            title(paste("cellline", cell.linei, "subsampling" ,coveragei, "minimum coverage", mincovi))
            grid(col = grey)

        }
        dev.off()
    }
}
setnames(table.dt, 'False-neg', 'FN')
table.dt %>% 
  ggplot(aes(x=centerkit,y=FN, fill=centerkit)) +
  geom_boxplot() + geom_jitter(width=0.1,alpha=0.2) + facet_wrap(~cell.line,ncol = 2) + theme(axis.text.x = element_text(angle = 90))
ggsave(file.path(outputDir,"results_FN.pdf"))

table.dt %>% 
  ggplot(aes(x=centerkit,y=Sensitivity, fill=centerkit)) +
  geom_boxplot() + geom_jitter(width=0.1,alpha=0.2) + theme(axis.text.x = element_text(angle = 90)) + ylim(0.98,1)
ggsave(file.path(outputDir,"results_sensitivity.pdf"))

table.dt %>% 
  ggplot(aes(x=centerkit,y=Precision, fill=centerkit)) +
  geom_boxplot() + geom_jitter(width=0.1,alpha=0.2) + theme(axis.text.x = element_text(angle = 90)) + ylim(0.95,1) 
ggsave(file.path(outputDir,"results_precision.pdf"))

