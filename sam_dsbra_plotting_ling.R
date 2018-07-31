#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))

# read in arguments
option_list <- list(
    make_option('--input', help='input file'),
    make_option('--mut_type', action='store_true', default=FALSE,
                help="plot mutation_event_frequency_by_type.csv,\
                      plot Mutation Event Frequency by Type"),
    make_option("--del_len", action='store_true', default=FALSE,
                help="plot deletion_lens.csv,\
                plot Frequency of Deletions by Length"),
    make_option("--del_seq", action='store_true', default=FALSE,
                help="plot sequences_with_deletion_events.csv\
                plot Sequences with Deletion Event"),
    make_option("--min_count", default=50,
                help="determine how the minimum counts to include in\
                      Sequences with Deletion/Insertion/Mutation Event"),
    make_option("--ins_len", action='store_true', default=FALSE,
                help="plot insertion_lens.csv\
                plot Frequency of Insertions by Length"),
    make_option("--ins_seq", action='store_true', default=FALSE,
                help="plot sequences_with_insertion_events.csv\
                plot Sequences with Insertion Event"),
    make_option("--repair_seq", action='store_true', default=FALSE,
                help="plot sequences_with_mutation_events.csv\
                plot Sequences with Mutation Event")
    )
opt <- parse_args(OptionParser(option_list=option_list))

# read input dataframe
df <- read.csv(opt$input)
basename <- tools::file_path_sans_ext(opt$input)

# output name
outname <- sprintf('%s.png', basename)
print(sprintf('Writing to %s', outname))

min_count <- opt$min_count

# Mutation Event Frequency by Type, pie
#if (opt$mut_type){
#    png(outname)
#    dft <- t(data.frame(count=df$count, row.names = df$mutation_type))
#    lbls <- paste(colnames(dft), "\n", df$count)
#    pie(dft, labels = lbls,
#        main="Mutation Event Frequency by Type")
#    dev.off()
#}

# Mutation Event Frequency by Type, bar
if (opt$mut_type){
    df$mutation_type <- factor(df$mutation_type, levels = df[order(df$count), ]$mutation_type)
    p <- ggplot(data=df, aes(x=mutation_type, y=count)) +
         geom_bar(stat="identity", width=.5) +
         labs(x='Mutation Type', y='Count', title='Mutation Event Frequency by Type') +
         geom_text(aes(label=count), position=position_dodge(width=0.9), vjust=-0.3, size=3) +
         theme_bw() +
         theme(aspect.ratio = 0.6)
    ggsave(filename=outname, plot=p)
}

# Frequency of Deletions by Length
if (opt$del_len){
    p <- ggplot(df, aes(x=deletion_length)) +
        geom_histogram(colour="black") +
        labs(x='Deletion Length (bp)', y='Count', title='Frequency of Deletions by Length') +
        theme_bw() +
        theme(aspect.ratio = 0.6)
    ggsave(filename=outname, plot=p)
}

# Sequences with deletion events
if (opt$del_seq){
    df$sequence <- factor(df$sequence, levels = df[order(df$count),]$sequence)
    df <- subset(df, count > min_count)
    p <- ggplot(data=df[order(-df$count),], aes(x=sequence, y=count)) +
    geom_bar(stat="identity", width=.5) +
    labs(x='Sequences', y='Count', title='Sequences With Deletion Events') +
    geom_text(aes(label=count), position=position_dodge(width=0.9), hjust=-0.3, size=4*20/nrow(df)) +
    theme_bw() +
    theme(aspect.ratio = 1 * nrow(df) / 20,
          axis.text=element_text(size=16*20/nrow(df))) +
    coord_flip(ylim = c(0, max(df$count) * 1.2))
    ggsave(filename=outname, plot=p)
}


# Frequency of Insertion by Length
if (opt$ins_len){
    p <- ggplot(df, aes(x=insertion_length)) +
        geom_histogram(colour="black") +
        labs(x='Insertion Length (bp)', y='Count', title='Frequency of Insertions by Length') +
        theme_bw() +
        theme(aspect.ratio = 0.6)
    ggsave(filename=outname, plot=p)
}

# Sequences with insertion events
if (opt$ins_seq){
    df$sequence <- factor(df$sequence, levels = df[order(df$count),]$sequence)
    df <- subset(df, count > min_count)
    p <- ggplot(data=df[order(-df$count),], aes(x=sequence, y=count)) +
    geom_bar(stat="identity", width=.5) +
    labs(x='Sequences', y='Count', title='Sequences With Insertion Events') +
    geom_text(aes(label=count), position=position_dodge(width=0.9), hjust=-0.3, size=3) +
    theme_bw() +
    theme(aspect.ratio = 0.5 * nrow(df) / 20) +
    coord_flip(ylim = c(0, max(df$count) * 1.2))
    tryCatch(ggsave(filename=outname, plot=p), error= function(err){print('insetion seq fails')})
}

# Sequences with mutation events
if (opt$repair_seq){
    df$sequence <- factor(df$sequence, levels = df[order(df$count),]$sequence)
    df <- subset(df, count > min_count)
    p <- ggplot(data=df[order(-df$count),], aes(x=sequence, y=count)) +
    geom_bar(stat="identity", width=.5) +
    labs(x='Sequences', y='Count', title='Repair Patterns') +
    geom_text(aes(label=count), position=position_dodge(width=0.9), hjust=-0.3, size=4*20/nrow(df)) +
    theme_bw() +
    theme(aspect.ratio = 1 * nrow(df) / 20,
          axis.text=element_text(size=16*20/nrow(df))) +
    coord_flip(ylim = c(0, max(df$count) * 1.2))
    tryCatch(ggsave(filename=outname, plot=p), error= function(err){print('repair patterns fails')})
}
